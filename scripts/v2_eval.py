"""Evaluate the calibrated v2 pipeline (slow x0.49, cardiac x0.28) against the earlier runs.

    python v2_eval.py

Slow drift
    SD of each unit's tracked position, raw vs every corrected run, on units tracked in all of
    them (matched by depth). Caveat (README, tracker-gain pitfall): the raw file's templates are
    the most smeared (lowest tracker gain -> SD under-read) and v2's the sharpest, so raw -> v2
    SD reductions are conservative; comparisons between corrected runs are close to like-for-like.
Cardiac
    calib_msit004_slowonly_vs_v2.npz (cardiac_calib.py measure, before = slow-only, after = v2;
    the units are detected on the slow-only file, as in the drift-free calibration). The units'
    cardiac motion T = zB / g uses the drift-free tracker gain g from the other half. The v2
    residual is |zA / gA|^2 / |T|^2 for gA = g (same sharpness) and gA = 1 (upper bound on gain,
    i.e. the smallest residual it could imply).
Writes v2_eval.json.
"""
import json
import os
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

os.environ.setdefault("NPX_REC", "msit004")
import cardiac_calib as cc  # noqa: E402
import recinfo  # noqa: E402

HERE = recinfo.RESULTS
OUT = recinfo.RESULTS
RUNS = [("slow x1.0 only", "m4_slowonly"), ("slow x1 + card x1.0", "m4_full"),
        ("slow x1 + card x0.6", "m4_final"), ("v2: slow x0.49 + card x0.28", "m4_v2")]


def slow_table():
    sd = {}
    for label, seg in RUNS:
        p = OUT / ("msit004_%s_tracking.npz" % seg)
        if not p.exists():
            print("missing %s" % p.name)
            continue
        d = np.load(p)
        for dep, sb, sa in zip(d["depth"], d["sd_before"], d["sd_after"]):
            if np.isfinite(sb) and np.isfinite(sa):
                sd.setdefault(float(dep), {})["raw_" + seg] = float(sb)
                sd[float(dep)][label] = float(sa)
    labels = [l for l, _ in RUNS]
    common = sorted(dep for dep, v in sd.items() if all(l in v for l in labels))
    print("slow drift: units tracked in all runs: %s um" % common)
    print("%-30s %s" % ("run", "  ".join("%6.0f" % c for c in common)) + "   median")
    raw = [np.nanmedian([v for k, v in sd[c].items() if k.startswith("raw_")]) for c in common]
    print("%-30s %s   %6.1f" % ("raw", "  ".join("%6.1f" % r for r in raw), np.median(raw)))
    out = dict(depths=common, raw=raw)
    for l in labels:
        vals = [sd[c][l] for c in common]
        out[l] = vals
        print("%-30s %s   %6.1f" % (l, "  ".join("%6.1f" % v for v in vals), np.median(vals)))
    return out


def cardiac():
    ref = np.load(HERE / "calib_msit004_driftfree_slow_vs_final.npz", allow_pickle=True)
    new = np.load(HERE / "calib_msit004_slowonly_vs_v2.npz", allow_pickle=True)
    if not np.array_equal(ref["depth"], new["depth"]):
        raise SystemExit("unit sets differ between calibrations")
    depth, zB0, zA0, ok0 = ref["depth"], ref["zB"], ref["zA"], ref["ok"]
    T0, Tmid, T1 = ref["T"]
    halves = [(T0, Tmid), (Tmid, T1)]
    q0 = cc.cardiac_phasors(str(ref["baseline"]), depth, halves, ref["mt"], ref["phase_u"])
    use0 = ok0 & np.isfinite(zA0) & np.isfinite(zB0) & np.isfinite(q0)
    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = use0[:, h]
        qa = 0.6 * q0[u, h]
        g[h] = np.sum(-(zA0[u, h] - zB0[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
    zB, zA, ok = new["zB"], new["zA"], new["ok"]
    use = use0 & ok & np.isfinite(zA) & np.isfinite(zB)
    print("\ncardiac: %d unit-halves; slow-only zB identical to drift-free calibration: max diff %.2g"
          % (use.sum(), np.nanmax(np.abs(zB[use] - zB0[use]))))
    Tp, A, F = [], [], []
    for h in (0, 1):
        u = use[:, h]
        Tp.append(zB[u, h] / g[1 - h])
        A.append(zA[u, h] / g[1 - h])
        F.append(zA0[u, h] / g[1 - h])
    Tp, A, F = (np.concatenate(v) for v in (Tp, A, F))
    res = lambda X, s=1.0: float(np.sum(np.abs(X / s) ** 2) / np.sum(np.abs(Tp) ** 2))
    gm = float(np.mean(np.abs(g)))
    p = wilcoxon(np.abs(Tp), np.abs(A), alternative="greater").pvalue
    print("residual cardiac power vs uncorrected (slow-only):")
    print("   final run (card x0.6): %.2f" % res(F))
    print("   v2 (card x0.28):       %.2f  (same tracker gain)   %.2f  (if v2 tracker gain were 1.0)"
          % (res(A), res(A, 1.0 / gm)))
    print("   v2 per unit-half: %d/%d reduced, median |A|/|T| %.2f, Wilcoxon p = %.3g"
          % (int((np.abs(A) < np.abs(Tp)).sum()), len(A), np.median(np.abs(A) / np.abs(Tp)), p))
    return dict(n=int(use.sum()), res_final=res(F), res_v2=res(A), res_v2_gain1=res(A, 1.0 / gm),
                n_reduced=int((np.abs(A) < np.abs(Tp)).sum()), ratio_median=float(np.median(np.abs(A) / np.abs(Tp))),
                wilcoxon_p=float(p), T_um=np.abs(Tp).tolist(), A_um=np.abs(A).tolist(), F_um=np.abs(F).tolist())


if __name__ == "__main__":
    out = dict(slow=slow_table(), cardiac=cardiac())
    (HERE / "v2_eval.json").write_text(json.dumps(out, indent=1))
    print("-> v2_eval.json")
