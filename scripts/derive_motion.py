"""Derive a band-selected, gain-scaled motion estimate from existing DREDge outputs.

    python derive_motion.py <src_tag> <new_tag> "<f_lo,f_hi,gain[,band_src_tag]>[;...]"

e.g.  (set NPX_REC=msit004 first)
    # slow drift + cardiac band at 35 % amplitude, both from one estimate
    python derive_motion.py msit004_nonrigid_800_t805-1745 msit004_n800_slowcard35 "0,0.1,1;0.8,5,0.35"

    # HYBRID: slow drift from the baseline, cardiac band from the LFP-high-passed estimate
    python derive_motion.py msit004_nonrigid_800_t805-1745 msit004_hyb_hp10card \
        "0,0.1,1;0.8,5,1,msit004_n800_hp10_t805-1745"

Each band is "f_lo,f_hi,gain" with an optional 4th field naming the estimate that band is
taken from (default: <src_tag>). f_lo = 0 means low-pass. The output is the sum of the
selected bands, each multiplied by its gain; every other frequency is dropped. All sources
must share window centres, start time and length.

Filtering is zero-phase over the FULL estimate with long even padding, so segments cut from
it later carry no filter edge transients (a 0.1 Hz filter applied per segment would ring for
tens of seconds at each segment boundary).

Rationale: apply only components that can be checked against the AP band - slow drift, and
the cardiac band with its amplitude verified per recording - rather than passing through
content the spikes cannot confirm. On msit004 the 0.1-0.5 Hz band sits on the MSIT trial
rate and cannot be separated from task-locked LFP change, so it is deliberately left out.
"""
import json
import sys

import numpy as np
from scipy.signal import butter, sosfiltfilt

import recinfo

MFS = 250.0


def _load(M, tag):
    D = np.loadtxt(M / (tag + "_displacement.csv"), delimiter=",")
    D = D[None, :] if D.ndim == 1 else D
    c = np.atleast_1d(np.loadtxt(M / (tag + "_window_centers.csv"), delimiter=","))
    t0f = M / (tag + "_t0.txt")
    t0 = float(t0f.read_text().strip()) if t0f.exists() else 0.0
    return D, c, t0


def derive(src, new, bands):
    M = recinfo.get()["motion"]
    base, c0, t00 = _load(M, src)
    cache = {src: base}
    out = np.zeros_like(base)
    pad = min(base.shape[1] - 1, int(60 * MFS))
    for band in bands:
        flo, fhi, gain = float(band[0]), float(band[1]), float(band[2])
        tag = band[3] if len(band) > 3 else src
        if tag not in cache:
            D, c, t0 = _load(M, tag)
            if D.shape != base.shape or not np.allclose(c, c0) or abs(t0 - t00) > 1e-6:
                raise SystemExit("%s does not match %s in windows/start/length" % (tag, src))
            cache[tag] = D
        if flo <= 0:
            sos = butter(3, fhi / (MFS / 2), btype="low", output="sos")
        else:
            sos = butter(3, [flo / (MFS / 2), fhi / (MFS / 2)], btype="band", output="sos")
        out += gain * sosfiltfilt(sos, cache[tag], axis=1, padtype="even", padlen=pad)
    np.savetxt(M / (new + "_displacement.csv"), out, delimiter=",", fmt="%.6f")
    for suffix in ("_window_centers.csv", "_t0.txt"):
        f = M / (src + suffix)
        if f.exists():
            (M / (new + suffix)).write_bytes(f.read_bytes())
    (M / (new + "_config.json")).write_text(json.dumps(dict(source=src, bands=bands)))
    print("%s -> %s   bands %s   median SD %.1f -> %.1f um"
          % (src, new, bands, float(np.median(base.std(1))), float(np.median(out.std(1)))))


def parse(spec):
    bands = []
    for part in spec.split(";"):
        fields = [f.strip() for f in part.split(",") if f.strip()]
        if not fields:
            continue
        band = [float(fields[0]), float(fields[1]), float(fields[2])]
        if len(fields) > 3:
            band.append(fields[3])
        bands.append(band)
    return bands


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    derive(sys.argv[1], sys.argv[2], parse(sys.argv[3]))
