"""Ground-truth injection test of the full MATLAB pipeline and of both trackers.

    python injection_test.py gen        # write the synthetic segment (msit004 AP folder / syn_card)
    matlab -batch run_interp_syn        # interpolate it with the real pipeline
    python injection_test.py analyze    # realised displacement per band, NO tracker involved
    python cardiac_calib.py measure syn_card 5 115      (NPX_REC=msit004)
    NPX_MIN_P2P_UM=0 python template_tracking.py syn_card 5 115 5.0
    python injection_test.py trackers   # each tracker's gain on a signal of known truth

Why
---
On msit004 the cardiac tracker reads only |g| ~ 0.35 of a known applied shift. That is
either the tracker under-reading, or the interpolation not realising cardiac-frequency
shifts (only slow shifts were verified so far: after-before slope -0.99 on msit001).
The two explanations lead to opposite conclusions about DREDge, so it is settled here with
ground truth.

Synthetic recording: 120 s, 384 channels in the msit004 OneBox layout (no sync channel),
Gaussian noise (sd 40 counts), per-channel DC offsets (sd 190 counts, like the real file),
and 12 STATIC units at depths 800-3000 um: biphasic spike, exp(-d/30 um) footprint, trough
400 counts, Poisson 20 Hz. Nothing in the input moves. The pipeline is then given
        h(t) = 40 sin(2 pi 0.02 t)  +  50 sin(2 pi 1.53 t)   um, identical at all depths
so every unit in the output should sit at y0 - h(t).
"""
import json
import os
import shutil
import sys
from pathlib import Path

os.environ["NPX_REC"] = "msit004"

import numpy as np  # noqa: E402
from scipy.signal import butter, filtfilt  # noqa: E402

import recinfo  # noqa: E402

SEG = "syn_card"
DUR = 120.0
FS = 30000
NCH = 384
MFS = 250.0
SLOW_AMP, SLOW_F = 40.0, 0.02
CARD_AMP, CARD_F = 50.0, 1.53
NOISE, OFFSET_SD, TROUGH, RATE, LAM = 40.0, 190.0, 400.0, 20.0, 30.0
UNIT_ROWS = list(range(40, 151, 10))
WIN_CENTERS = [710.0, 1510.0, 2310.0, 3110.0]
T_ANA = (5.0, 115.0)
HERE = recinfo.RESULTS


def h_slow(t):
    return SLOW_AMP * np.sin(2 * np.pi * SLOW_F * t)


def h_card(t):
    return CARD_AMP * np.sin(2 * np.pi * CARD_F * t)


def seg_dir():
    return recinfo.get()["ap"] / SEG


def gen():
    S = seg_dir()
    S.mkdir(exist_ok=True)
    yc_ch, xc_ch = recinfo.chanmap_yx()
    rng = np.random.default_rng(7)
    nfr = int(DUR * FS)
    chunk = FS
    offsets = rng.normal(0, OFFSET_SD, NCH).astype(np.float32)
    kk = np.arange(-20, 41)
    wf = ((-np.exp(-(kk / 4.0) ** 2) + 0.35 * np.exp(-((kk - 12) / 8.0) ** 2)) * TROUGH).astype(np.float32)
    units = []
    for k, row in enumerate(UNIT_ROWS):
        x0, y0 = [0.0, 32.0][k % 2], row * 20.0
        foot = np.exp(-np.hypot(xc_ch - x0, yc_ch - y0) / LAM).astype(np.float32)
        isi = rng.exponential(1.0 / RATE, int(DUR * RATE * 1.6)) + 0.002
        st = (np.cumsum(isi) * FS).astype(np.int64)
        st = st[st < nfr - 100]
        st = st[((st % chunk) >= 20) & ((st % chunk) < chunk - 41)]   # keep whole waveforms
        units.append(dict(row=row, x=x0, y=y0, foot=foot, st=st))
    out = np.memmap(S / "continuous.dat", dtype="int16", mode="w+", shape=(nfr, NCH))
    for c0 in range(0, nfr, chunk):
        c1 = min(nfr, c0 + chunk)
        blk = rng.normal(0, NOISE, (c1 - c0, NCH)).astype(np.float32) + offsets
        for u in units:
            for s in u["st"][(u["st"] >= c0) & (u["st"] < c1)] - c0:
                blk[s - 20:s + 41] += wf[:, None] * u["foot"][None, :]
        out[c0:c1] = np.clip(np.round(blk), -32768, 32767).astype(np.int16)
    out.flush()
    del out
    t = np.arange(int(DUR * MFS) + 8) / MFS
    h = h_slow(t) + h_card(t)
    np.savetxt(S / "nonrigid_displacement.csv", np.tile(h, (len(WIN_CENTERS), 1)), delimiter=",", fmt="%.6f")
    np.savetxt(S / "dredge_window_centers.csv", np.array(WIN_CENTERS), delimiter=",", fmt="%.4f")
    np.savetxt(S / "ref_displacement.csv", np.zeros(len(WIN_CENTERS)), delimiter=",", fmt="%.6f")
    (S / "segment_t0.txt").write_text("0.000000")
    shutil.copy(recinfo.get()["chanmap"], S / "chanMap.mat")
    np.savez(S / "truth.npz", rows=np.array(UNIT_ROWS), x=np.array([u["x"] for u in units]),
             y=np.array([u["y"] for u in units]),
             spikes=np.array([u["st"] for u in units], dtype=object),
             params=json.dumps(dict(DUR=DUR, SLOW_AMP=SLOW_AMP, SLOW_F=SLOW_F, CARD_AMP=CARD_AMP,
                                    CARD_F=CARD_F, NOISE=NOISE, TROUGH=TROUGH, RATE=RATE, LAM=LAM)))
    print("synthetic segment written: %s  (%d units, %d spikes, %.0f s)"
          % (S, len(units), sum(len(u["st"]) for u in units), DUR))


def centroid_positions(path, truth, yc, RMAT, win=10):
    nfr = path.stat().st_size // (NCH * 2)
    mm = np.memmap(path, dtype="int16", mode="r", shape=(nfr, NCH))
    bp = butter(3, [300 / (FS / 2), 6000 / (FS / 2)], btype="band")
    pos = {k: ([], []) for k in range(len(truth["rows"]))}
    for b0 in np.arange(T_ANA[0], T_ANA[1], 10.0):
        s0, s1 = int(b0 * FS), int(min(T_ANA[1], b0 + 10.0) * FS)
        x = np.asarray(mm[s0 - 200:s1 + 200]).astype(np.float32)
        x = filtfilt(*bp, x, axis=0)
        x = (x - np.median(x, axis=1, keepdims=True)) @ RMAT
        for k, row in enumerate(truth["rows"]):
            st = truth["spikes"][k]
            st = st[(st >= s0) & (st < s1)]
            if not len(st):
                continue
            rr = np.arange(row - win, row + win + 1)
            idx = (st - s0 + 200)[:, None] + np.arange(-2, 3)[None, :]
            trough = x[idx][:, :, rr].min(axis=1)                 # (nspk, 2*win+1)
            w = np.clip(-trough, 0, None) ** 2
            cen = (w * rr[None, :]).sum(1) / np.maximum(w.sum(1), 1e-9)
            pos[k][0].append(st / FS)
            pos[k][1].append(cen * float(np.median(np.diff(yc))))
    return {k: (np.concatenate(v[0]), np.concatenate(v[1])) for k, v in pos.items() if v[0]}


def fit_bands(t, y):
    X = np.c_[np.ones_like(t), h_slow(t), np.cos(2 * np.pi * CARD_F * t), np.sin(2 * np.pi * CARD_F * t)]
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    z = complex(c[2], -c[3])                          # y_card = Re(z * exp(i*2*pi*f*t))
    slow_gain = -c[1]                                 # expected +1 (position = y0 - h)
    card_gain = z / (1j * CARD_AMP)                   # -h_card = Re(50i * e^{i phi}) -> expected 1+0j
    return slow_gain, card_gain


def analyze():
    S = seg_dir()
    truth = np.load(S / "truth.npz", allow_pickle=True)
    yc, RMAT = recinfo.row_matrix()
    # the centroid estimator's own gain, on noiseless footprints of known shift
    yc_ch, xc_ch = recinfo.chanmap_yx()
    est = []
    for k, row in enumerate(truth["rows"]):
        for d in (-40.0, 40.0):
            f = -np.exp(-np.hypot(xc_ch - truth["x"][k], yc_ch - (truth["y"][k] + d)) / LAM) @ RMAT
            rr = np.arange(row - 10, row + 11)
            w = np.clip(-f[rr], 0, None) ** 2
            est.append(((w * rr).sum() / w.sum() * 20.0 - truth["y"][k]) / d)
    print("centroid estimator gain on noiseless shifted footprints: %.3f (1 = unbiased)" % np.median(est))
    res = {}
    for lab, path in [("before", S / "continuous.dat"), ("after", S / "interpsyn_card" / "interpsyn_card.dat")]:
        pos = centroid_positions(path, truth, yc, RMAT)
        sg, cg = [], []
        for k, (t, y) in pos.items():
            a, z = fit_bands(t, y)
            sg.append(a)
            cg.append(z)
        sg, cg = np.array(sg), np.array(cg)
        res[lab] = dict(slow=float(np.median(sg)), card_abs=float(np.median(np.abs(cg))),
                        card_deg=float(np.degrees(np.angle(np.median(cg.real) + 1j * np.median(cg.imag)))))
        print("%-7s realised gain (centroid, no tracker), median over %d units:  slow %.3f   "
              "cardiac |%.3f| at %+.0f deg" % (lab, len(sg), res[lab]["slow"], res[lab]["card_abs"], res[lab]["card_deg"]))
    print("expected: before 0 / 0 (input static);  after 1.00 / 1.00 at 0 deg (pipeline exact)")
    (HERE / "injection_realised.json").write_text(json.dumps(res, indent=1))


def trackers():
    real = json.loads((HERE / "injection_realised.json").read_text())
    ra = real["after"]
    cal = np.load(HERE / "calib_msit004_syn_card.npz", allow_pickle=True)
    zA = cal["zA"][np.isfinite(cal["zA"])]
    zB = cal["zB"][np.isfinite(cal["zB"])]
    # calib's phase reference is the Hilbert phase of the synthetic cardiac trace, in which
    # h_card has phasor CARD_AMP at 0 deg; content moves by -h, so A = -g * CARD_AMP
    g = -np.median(zA.real) / CARD_AMP + 1j * (-np.median(zA.imag) / CARD_AMP)
    print("cardiac tracker (cardiac_calib): |B| median %.1f um (truth 0); combined gain %.2f at %+.0f deg"
          % (np.median(np.abs(zB)), abs(g), np.degrees(np.angle(g))))
    print("   -> tracker-only gain = combined / realised = %.2f" % (abs(g) / max(ra["card_abs"], 1e-9)))
    tt = np.load(HERE.parent / ("msit004_%s_tracking.npz" % SEG))
    sa = tt["slope_after"][np.isfinite(tt["slope_after"])]
    sb = tt["slope_before"][np.isfinite(tt["slope_before"])]
    print("slow tracker (template_tracking): slope before %.2f (truth 0), after %.2f -> combined gain %.2f"
          % (np.median(sb), np.median(sa), -np.median(sa)))
    print("   -> tracker-only gain = combined / realised = %.2f" % (-np.median(sa) / max(ra["slow"], 1e-9)))


if __name__ == "__main__":
    {"gen": gen, "analyze": analyze, "trackers": trackers}[sys.argv[1]]()
