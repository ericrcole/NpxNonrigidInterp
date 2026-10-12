"""Score candidate DREDge cardiac estimates against units' DRIFT-FREE cardiac motion.

    python score_driftfree.py <tag> [<tag> ...]      (NPX_REC=msit004)

Supersedes the amplitude columns of `cardiac_calib.py score`, whose calibration compared the RAW
file (unit templates smeared by +/-100 um of slow drift, tracker gain ~0.35) with a corrected file
(sharp templates) and so mis-scaled every unit's cardiac motion by ~3x. The phase and slow-drift
columns of that table are gain-invariant and remain valid.

Here the calibration is the drift-free pair (slow-only vs final, which adds 0.6 x the baseline
cardiac estimate q0): tracker gain g is fitted on one half from A - B = -0.6 g q0 (|g| = 0.93-0.96,
matching synthetic ground truth), and the units' cardiac motion T = B / g is used on the other half.
For each candidate cardiac estimate q_c:
    |k|      gain of T on q_c (ideal 1)          phase    timing of T relative to q_c
    RES      residual cardiac power if q_c is applied as-is (<1 helps)
    RES_opt  residual with the best real gain applied to q_c
Writes score_msit004_driftfree.json.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("NPX_REC", "msit004")
import cardiac_calib as cc  # noqa: E402
import recinfo  # noqa: E402

HERE = recinfo.RESULTS
NPZ = HERE / "calib_msit004_driftfree_slow_vs_final.npz"
APPLIED_GAIN = 0.6


def main(tags):
    cal = np.load(NPZ, allow_pickle=True)
    depth, zB, zA, ok = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
    T0, Tmid, T1 = cal["T"]
    halves = [(T0, Tmid), (Tmid, T1)]
    mt, phase_u = cal["mt"], cal["phase_u"]
    base = str(cal["baseline"])
    q0 = cc.cardiac_phasors(base, depth, halves, mt, phase_u)
    use = ok & np.isfinite(zA) & np.isfinite(zB) & np.isfinite(q0)
    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = use[:, h]
        qa = APPLIED_GAIN * q0[u, h]
        g[h] = np.sum(-(zA[u, h] - zB[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
    print("drift-free tracker gain |g| = %.2f / %.2f (halves); %d usable unit-halves" % (abs(g[0]), abs(g[1]), use.sum()))
    print("%-30s %7s %7s | %6s %7s %6s %8s" % ("candidate", "|q|um", "|T|um", "|k|", "phase", "RES", "RES_opt"))
    rows = []
    for tag in [base] + [t for t in tags if t != base]:
        try:
            q = cc.cardiac_phasors(tag, depth, halves, mt, phase_u)
        except OSError as exc:
            print("%-30s missing (%s)" % (tag, exc))
            continue
        per = []
        for h in (0, 1):
            u = use[:, h] & np.isfinite(q[:, h])
            T = zB[u, h] / g[1 - h]
            qc = q[u, h]
            k = np.sum(T * np.conj(qc)) / np.sum(np.abs(qc) ** 2)
            c = float(np.real(k))
            per.append((np.median(np.abs(qc)), np.median(np.abs(T)), abs(k), np.degrees(np.angle(k)),
                        np.sum(np.abs(T - qc) ** 2) / np.sum(np.abs(T) ** 2),
                        np.sum(np.abs(T - c * qc) ** 2) / np.sum(np.abs(T) ** 2)))
        r = np.nanmean(np.array(per), axis=0)
        name = tag.replace("msit004_", "").replace("_t805-1745", "")
        print("%-30s %7.1f %7.1f | %6.2f %+6.0f  %6.2f %8.2f" % (name, *r))
        rows.append(dict(tag=tag, q_um=float(r[0]), T_um=float(r[1]), k=float(r[2]), k_deg=float(r[3]),
                         RES=float(r[4]), RES_opt=float(r[5])))
    (HERE / "score_msit004_driftfree.json").write_text(json.dumps(dict(g=[float(abs(x)) for x in g], candidates=rows), indent=1))
    print("-> score_msit004_driftfree.json")


if __name__ == "__main__":
    main(sys.argv[1:])
