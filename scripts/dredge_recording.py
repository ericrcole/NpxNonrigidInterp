"""Run DREDge on the LFP stream of a recording, over a chosen time range.

    python -u dredge_recording.py <rec_key> <t0_s> <t1_s> [config ...]

rec_key    'msit001'  Neuropix_2026-03-05_11-20-16  (Neuropix-PXI, 385 ch AP incl. sync)
           'msit004'  NPIX_MSIT_004                 (OneBox, 384 ch AP, no sync channel)

Times are FILE-RELATIVE seconds (0 = first sample), not the absolute Open Ephys clock.

Preprocessing matches oep_spikeinterface_pipeline.ipynb cell 2:
    bandpass fmin-250 -> phase_shift -> resample 250 Hz -> average_across_direction
with fmin = 0.5 Hz unless a config overrides it. average_across_direction collapses the
column pairs to unique depths, which dredge_lfp requires.

Window geometry: estimate_motion() defaults to win_step_um=200/win_scale_um=300, which
override dredge_online_lfp's own 800/850. On msit001 the 200 um setting gave physically
impossible tissue strain between adjacent windows (median 0.26). Prefer 800 um.

Variants for the cardiac-amplitude question
-------------------------------------------
On msit004 the baseline estimate gets cardiac TIMING right but amplitude ~2.5-3.5x too
large relative to the motion of nearby units. Two DREDge-side hypotheses:
  fmin       high-pass the LFP before registration. Pulsatile (blood-volume) LFP changes
             are low-frequency; removing that band may remove the inflation.
  thomas_kw  solver priors. lambda_t (default 1.0) penalises frame-to-frame change, so
             raising it shrinks fast displacement; lambda_s does the same across windows.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import spikeinterface.full as si
from scipy.signal import butter, filtfilt
from spikeinterface.sortingcomponents.motion import estimate_motion

import recinfo  # noqa: E402


def _rec(key):
    r = recinfo.get(key)
    return dict(base=r["lfp_base"], stream=r["lfp_stream"], outdir=r["motion"],
                max_depth_um=r.get("max_depth_um"))

W800 = dict(rigid=False, win_step_um=800.0, win_scale_um=850.0)
CONFIGS = {
    "nonrigid_800": dict(W800),
    "nonrigid_200": dict(rigid=False, win_step_um=200.0, win_scale_um=300.0),
    "nonrigid_400": dict(rigid=False, win_step_um=400.0, win_scale_um=600.0),
    "rigid": dict(rigid=True),
    "n800_hp10": dict(W800, fmin=10.0),
    "n800_hp20": dict(W800, fmin=20.0),
    "n800_lt10": dict(W800, thomas_kw=dict(lambda_t=10.0)),
    "n800_lt100": dict(W800, thomas_kw=dict(lambda_t=100.0)),
}
LFP_FS = 250
_CACHE = {}


def build(key, t0, t1, fmin):
    ck = (key, t0, t1, fmin)
    if ck in _CACHE:
        return _CACHE[ck]
    r = _rec(key)
    raw = si.read_openephys(r["base"], stream_name=r["stream"], load_sync_channel=False)
    rec = si.bandpass_filter(raw, freq_min=fmin, freq_max=250, dtype="float32",
                             ignore_low_freq_error=True)
    rec = si.phase_shift(rec)
    rec = si.resample(rec, resample_rate=LFP_FS, margin_ms=1000)
    rec = si.average_across_direction(rec)
    if r["max_depth_um"] is not None:
        y = rec.get_channel_locations()[:, 1]
        rec = rec.select_channels(rec.channel_ids[y <= r["max_depth_um"]])
    n = rec.get_num_frames()
    f0, f1 = int(round(t0 * LFP_FS)), min(int(round(t1 * LFP_FS)), n)
    rec = rec.frame_slice(start_frame=f0, end_frame=f1)
    y = rec.get_channel_locations()[:, 1]
    print("%s fmin=%.1f: %d ch, %.0f-%.0f um, t=%.0f-%.0f s (%d frames @ %d Hz)"
          % (key, fmin, rec.get_num_channels(), y.min(), y.max(), t0, f1 / LFP_FS,
             rec.get_num_frames(), LFP_FS), flush=True)
    _CACHE[ck] = (rec, r["outdir"])
    return _CACHE[ck]


def hf_rms(d):
    lp = filtfilt(*butter(3, 0.1 / (LFP_FS / 2), btype="low"), d, axis=0)
    return float(np.median((d - lp).std(0)))


def run(key, name, t0, t1):
    cfg = dict(CONFIGS[name])
    fmin = cfg.pop("fmin", 0.5)
    rec, outdir = build(key, t0, t1, fmin)
    kw = dict(method="dredge_lfp", progress_bar=False, **cfg)
    print("\n=== %s / %s === fmin=%.1f %s" % (key, name, fmin, cfg), flush=True)
    tic = time.time()
    m = estimate_motion(rec, **kw)
    if isinstance(m, tuple):
        m = m[0]
    d = np.asarray(m.displacement[0])
    c = np.atleast_1d(np.asarray(m.spatial_bins_um)).ravel()
    print("  %.0f s -> %s, %d window(s) at %s" % (time.time() - tic, d.shape, len(c),
                                                  np.round(c).astype(int)), flush=True)
    outdir.mkdir(parents=True, exist_ok=True)
    tag = "%s_%s_t%d-%d" % (key, name, int(t0), int(t0 + d.shape[0] / LFP_FS))
    np.savetxt(outdir / (tag + "_displacement.csv"), d.T, delimiter=",", fmt="%.6f")
    np.savetxt(outdir / (tag + "_window_centers.csv"), c, delimiter=",", fmt="%.4f")
    (outdir / (tag + "_t0.txt")).write_text("%.6f" % t0)
    (outdir / (tag + "_config.json")).write_text(json.dumps(dict(CONFIGS[name]), default=str))
    if d.shape[1] > 1:
        dz = np.median(np.diff(c))
        st = np.abs(np.diff(d.T, axis=0) / dz)
        print("  strain: median %.3f, %.1f%% > 0.5" % (np.median(st), 100 * np.mean(st > 0.5)))
    card = filtfilt(*butter(3, [0.8 / (LFP_FS / 2), 2.0 / (LFP_FS / 2)], btype="band"), d, axis=0)
    print("  SD %.1f um   RMS>0.1Hz %.1f um   cardiac-band RMS %.1f um"
          % (np.median(d.std(0)), hf_rms(d), float(np.median(card.std(0)))), flush=True)
    print("  saved %s" % tag, flush=True)


if __name__ == "__main__":
    key, t0, t1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    names = sys.argv[4:] or ["nonrigid_800"]
    for n in names:
        try:
            run(key, n, t0, t1)
        except Exception as exc:  # one failing variant must not kill the sweep
            print("  !! %s failed: %r" % (n, exc), flush=True)
    print("\nDONE", flush=True)
