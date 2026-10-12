"""Per-recording gain calibration of the slow and cardiac motion bands, held out by half.

    python calibrate_gains.py

Why gains at all
----------------
The interpolation applies ~100% of whatever displacement it is given (chunk registration on
the real recording: slow x1.00, cardiac x0.96-0.99). Whether that removes motion depends on
whether DREDge's displacement matches how the units actually move. On msit004 the tracked
position of units correlates r ~ +0.8 with DREDge's slow estimate before correction but r ~ -0.5
after - the correction overshoots. This script measures by how much, per band.

Slow band
    template_tracking.py on segment m4_slowonly compares the RAW file with the SLOW-ONLY
    interpolation. Both leave cardiac and fast motion uncorrected, so unit templates are equally
    smeared and the tracker gain g is the same in both files; they differ only by the slow shift:
        slope_after - slope_before = -g
        s = slope_before / g        units move s x DREDge's slow displacement
Cardiac band
    cardiac_calib.py on the drift-free pair (slow-only vs final, which adds 0.6 x cardiac):
        A - B = -g * 0.6 * q0,   T = B / g,   k = <T, q0> / <q0, q0>

Each gain is estimated on one half of the valid range and its consequence predicted on the
other. Writes gains_msit004.json.
"""
import json
import os
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import butter, filtfilt, resample_poly

os.environ.setdefault("NPX_REC", "msit004")
import cardiac_calib as cc  # noqa: E402
import recinfo  # noqa: E402

HERE = recinfo.RESULTS
OUT = recinfo.RESULTS
FS, NT = 30000.0, 128
FSB = FS / NT


def slow_band():
    npz = OUT / "msit004_m4_slowonly_tracking.npz"
    if not npz.exists():
        print("slow: %s not available yet" % npz.name)
        return None
    d = np.load(npz)
    ctr, depth, Yb, Ya = d["ctr"], d["depth"], d["disp_before"], d["disp_after"]
    S = recinfo.get()["ap"] / "m4_slowonly"
    seg_t0 = float((S / "segment_t0.txt").read_text())
    D = np.loadtxt(S / "nonrigid_displacement.csv", delimiter=",")
    cwin = np.loadtxt(S / "dredge_window_centers.csv", delimiter=",")
    h = resample_poly(D.T, 120, 128, axis=0)
    h = filtfilt(*butter(3, recinfo.get()["motion_lowpass_hz"] / (FSB / 2), btype="low"), h, axis=0)
    bt = seg_t0 + np.arange(h.shape[0]) / FSB
    width = float(np.median(np.diff(ctr)))
    edges = np.r_[ctr - width / 2, ctr[-1] + width / 2]
    bi = np.digitize(bt, edges) - 1
    inb = (bi >= 0) & (bi < len(ctr))
    cnt = np.bincount(bi[inb], minlength=len(ctr))[:len(ctr)]
    X = PchipInterpolator(cwin, h, axis=1)(np.clip(depth, cwin[0], cwin[-1]))   # (ntime_batches, nunits)
    mid = 0.5 * (ctr[0] + ctr[-1])
    res = {}
    for hname, sel in (("half1", ctr < mid), ("half2", ctr >= mid)):
        sb, sa, rb, ra = [], [], [], []
        for u in range(len(depth)):
            xb = np.bincount(bi[inb], weights=X[inb, u], minlength=len(ctr))[:len(ctr)] / np.maximum(cnt, 1)
            fit = []
            for Y in (Yb, Ya):
                m = sel & np.isfinite(Y[u]) & (cnt > 0)
                if m.sum() < 8:
                    fit = None
                    break
                xx = xb[m] - xb[m].mean()
                yy = Y[u][m] - Y[u][m].mean()
                fit.append((np.sum(xx * yy) / np.sum(xx * xx), np.corrcoef(xx, yy)[0, 1]))
            if fit is None:
                continue
            sb.append(fit[0][0]); rb.append(fit[0][1])
            sa.append(fit[1][0]); ra.append(fit[1][1])
        sb, sa = np.array(sb), np.array(sa)
        g = float(np.median(-(sa - sb)))
        s = float(np.median(sb) / g) if g > 0 else np.nan
        res[hname] = dict(n=int(len(sb)), tracker_gain=g, s=s, slope_before=float(np.median(sb)),
                          slope_after=float(np.median(sa)), r_before=float(np.median(rb)), r_after=float(np.median(ra)))
        print("slow %s: n=%d  tracker gain g=%.2f  slope %.2f -> %.2f  r %.2f -> %.2f  => units move s=%.2f x DREDge"
              % (hname, len(sb), g, np.median(sb), np.median(sa), np.median(rb), np.median(ra), s))
    s1, s2 = res["half1"]["s"], res["half2"]["s"]
    for cal, ev, sc, se in (("half1", "half2", s1, s2), ("half2", "half1", s2, s1)):
        print("   calibrate on %s (gain %.2f), evaluate on %s: residual fraction %.2f  (gain 1.0 would leave %.2f)"
              % (cal, sc, ev, abs(se - sc) / se, abs(se - 1.0) / se))
    res["recommended_gain"] = float(np.mean([s1, s2]))
    return res


def cardiac_band():
    npz = HERE / "calib_msit004_driftfree_slow_vs_final.npz"
    if not npz.exists():
        print("cardiac: %s not available yet" % npz.name)
        return None
    cal = np.load(npz, allow_pickle=True)
    depth, zB, zA, ok = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
    T0, Tmid, T1 = cal["T"]
    halves = [(T0, Tmid), (Tmid, T1)]
    gain_applied = 0.6
    q0 = cc.cardiac_phasors(str(cal["baseline"]), depth, halves, cal["mt"], cal["phase_u"])
    use = ok & np.isfinite(zA) & np.isfinite(zB) & np.isfinite(q0)
    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = use[:, h]
        qa = gain_applied * q0[u, h]
        g[h] = np.sum(-(zA[u, h] - zB[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
    res, ks = {}, []
    for h in (0, 1):
        u = use[:, h]
        T = zB[u, h] / g[1 - h]
        Q = q0[u, h]
        k = np.sum(T * np.conj(Q)) / np.sum(np.abs(Q) ** 2)
        ks.append(k)
        res["half%d" % (h + 1)] = dict(n=int(u.sum()), tracker_gain=float(abs(g[h])), k_abs=float(abs(k)),
                                       k_deg=float(np.degrees(np.angle(k))), k_real=float(np.real(k)))
        print("cardiac half%d: n=%d  tracker gain |g|=%.2f  => units move |k|=%.2f x DREDge at %+.0f deg (Re k %.2f)"
              % (h + 1, u.sum(), abs(g[h]), abs(k), np.degrees(np.angle(k)), np.real(k)))
    for cal_h, ev_h in ((0, 1), (1, 0)):
        u = use[:, ev_h]
        T = zB[u, ev_h] / g[cal_h]
        Q = q0[u, ev_h]
        c = float(np.real(ks[cal_h]))
        r_c = np.sum(np.abs(T - c * Q) ** 2) / np.sum(np.abs(T) ** 2)
        r_1 = np.sum(np.abs(T - Q) ** 2) / np.sum(np.abs(T) ** 2)
        r_6 = np.sum(np.abs(T - gain_applied * Q) ** 2) / np.sum(np.abs(T) ** 2)
        print("   calibrate on half%d (gain %.2f), evaluate on half%d: residual power %.2f  (gain 0.6: %.2f, gain 1.0: %.2f)"
              % (cal_h + 1, c, ev_h + 1, r_c, r_6, r_1))
    res["recommended_gain"] = float(np.mean([np.real(k) for k in ks]))
    return res


if __name__ == "__main__":
    out = dict(slow=slow_band(), cardiac=cardiac_band())
    (HERE / "gains_msit004.json").write_text(json.dumps(out, indent=1))
    print("-> gains_msit004.json")
