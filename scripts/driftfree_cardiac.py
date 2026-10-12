"""Drift-free cardiac evaluation: slow-only interpolation (before) vs final pipeline (after).

    python driftfree_cardiac.py

Reads calib_msit004_driftfree_slow_vs_final.npz, written by

    NPX_CALIB_BEFORE=<m4_slowonly output>  NPX_CALIB_AFTER=<m4_final output>
    NPX_CALIB_TAG=driftfree_slow_vs_final  python cardiac_calib.py measure m4_final 810 1740

Why this pair
-------------
With the raw recording as "before", each unit's grand-average template is smeared by
+/-100 um of slow drift while the corrected "after" template is sharp, so the two files carry
different tracker gains and a single calibration gain is wrong. Here BOTH files have slow drift
removed identically, and they differ ONLY by the cardiac correction the final pipeline added:
CARD_GAIN x the baseline cardiac estimate q0. Per unit-half:

    A - B = -g * CARD_GAIN * q0         g = tracker gain, now the same in both files
    T     = B / g                       unit's own cardiac motion
    k     = <T, q0> / <q0, q0>          units move k x DREDge's cardiac estimate
    observed residual power = sum|A|^2 / sum|B|^2      (g cancels - no calibration needed)

g is fitted on one half of the recording and applied to the other. Writes driftfree_cardiac.json.
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
NPZ = HERE / "calib_msit004_driftfree_slow_vs_final.npz"
FINAL_TAG = "msit004_n800_slowcard60"


def card_gain():
    cfg = json.loads((recinfo.get()["motion"] / (FINAL_TAG + "_config.json")).read_text())
    for band in cfg["bands"]:
        if float(band[0]) > 0:
            return float(band[2])
    raise SystemExit("no cardiac band in %s" % FINAL_TAG)


def main():
    cal = np.load(NPZ, allow_pickle=True)
    depth, zB, zA, ok = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
    T0, Tmid, T1 = cal["T"]
    halves = [(T0, Tmid), (Tmid, T1)]
    gain = card_gain()
    q0 = cc.cardiac_phasors(str(cal["baseline"]), depth, halves, cal["mt"], cal["phase_u"])
    use = ok & np.isfinite(zA) & np.isfinite(zB) & np.isfinite(q0)
    print("drift-free pair: %d usable unit-halves; final cardiac gain %.2f" % (use.sum(), gain))

    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = use[:, h]
        qa = gain * q0[u, h]
        g[h] = np.sum(-(zA[u, h] - zB[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
        print("  half %d: tracker gain |g| = %.2f at %+.0f deg (n=%d)" % (h + 1, abs(g[h]), np.degrees(np.angle(g[h])), u.sum()))

    Tl, Ql, Bl, Al = [], [], [], []
    for h in (0, 1):
        u = use[:, h]
        gc = g[1 - h]
        Tl.append(zB[u, h] / gc)
        Ql.append(q0[u, h])
        Bl.append(zB[u, h])
        Al.append(zA[u, h])
    T, Q, B, A = (np.concatenate(v) for v in (Tl, Ql, Bl, Al))
    k = np.sum(T * np.conj(Q)) / np.sum(np.abs(Q) ** 2)
    c_opt = float(np.real(k))
    res_pred = lambda c: float(np.sum(np.abs(T - c * Q) ** 2) / np.sum(np.abs(T) ** 2))
    res_obs = float(np.sum(np.abs(A) ** 2) / np.sum(np.abs(B) ** 2))
    try:
        p = wilcoxon(np.abs(B), np.abs(A), alternative="greater").pvalue
    except ValueError:
        p = np.nan
    print("")
    print("DREDge accuracy: units move |k| = %.2f x DREDge's cardiac estimate, phase %+.0f deg" % (abs(k), np.degrees(np.angle(k))))
    print("   median unit cardiac motion |T| %.1f um vs DREDge |q0| %.1f um" % (np.median(np.abs(T)), np.median(np.abs(Q))))
    print("predicted residual cardiac power: gain 1.0 -> %.2f, gain %.2f (applied) -> %.2f, gain %.2f (optimal real) -> %.2f"
          % (res_pred(1.0), gain, res_pred(gain), c_opt, res_pred(c_opt)))
    print("OBSERVED residual cardiac power, final vs slow-only: %.2f  (no calibration involved)" % res_obs)
    print("   per unit-half |A|/|B|: median %.2f; %d/%d reduced; Wilcoxon |B| > |A| p = %.3g"
          % (np.median(np.abs(A) / np.abs(B)), int((np.abs(A) < np.abs(B)).sum()), len(A), p))
    out = dict(n=int(use.sum()), card_gain=gain,
               g=[[float(abs(x)), float(np.degrees(np.angle(x)))] for x in g],
               k_abs=float(abs(k)), k_deg=float(np.degrees(np.angle(k))),
               T_median_um=float(np.median(np.abs(T))), q0_median_um=float(np.median(np.abs(Q))),
               res_pred_gain1=res_pred(1.0), res_pred_applied=res_pred(gain),
               c_opt=c_opt, res_pred_opt=res_pred(c_opt), res_observed=res_obs,
               ratio_median=float(np.median(np.abs(A) / np.abs(B))),
               n_reduced=int((np.abs(A) < np.abs(B)).sum()), wilcoxon_p=float(p))
    (HERE / "driftfree_cardiac.json").write_text(json.dumps(out, indent=1))
    print("-> driftfree_cardiac.json")


if __name__ == "__main__":
    main()
