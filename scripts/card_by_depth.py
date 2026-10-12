"""Per-unit cardiac motion vs depth: is msit004's cardiac motion uniform along the probe?

    python card_by_depth.py

Panel E of the summary is bimodal (units at ~1x DREDge and at ~0 um). Pooled gains (0.27x) mix
the two. For every unit-half of the drift-free calibration this prints the unit's own cardiac motion
T = B/g (tracker gain from the other half) next to DREDge's cardiac estimate at that depth from the
800 um and 200 um window runs, and the per-unit ratio and phase. Writes card_by_depth.png/.json.
"""
import json
import os
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.environ.setdefault("NPX_REC", "msit004")
import cardiac_calib as cc  # noqa: E402
import recinfo  # noqa: E402

HERE = recinfo.RESULTS
cal = np.load(HERE / "calib_msit004_driftfree_slow_vs_final.npz", allow_pickle=True)
depth, zB, zA, ok = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
T0, Tmid, T1 = cal["T"]
halves = [(T0, Tmid), (Tmid, T1)]
q8 = cc.cardiac_phasors(str(cal["baseline"]), depth, halves, cal["mt"], cal["phase_u"])
q2 = cc.cardiac_phasors("msit004_nonrigid_200_t805-1745", depth, halves, cal["mt"], cal["phase_u"])
use = ok & np.isfinite(zA) & np.isfinite(zB) & np.isfinite(q8)
g = np.full(2, np.nan + 0j)
for h in (0, 1):
    u = use[:, h]
    qa = 0.6 * q8[u, h]
    g[h] = np.sum(-(zA[u, h] - zB[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
T = zB / g[::-1][None, :]
print("%7s %5s | %7s %7s | %6s %6s | %7s %6s" % ("depth", "nspk", "|T|h1", "|T|h2", "|q800|", "ratio", "|q200|", "ratio"))
rows = []
for i in np.argsort(depth):
    if not use[i].any():
        continue
    th = [abs(T[i, h]) if use[i, h] else np.nan for h in (0, 1)]
    r8 = np.nanmean([abs(T[i, h]) / abs(q8[i, h]) for h in (0, 1) if use[i, h]])
    r2 = np.nanmean([abs(T[i, h]) / abs(q2[i, h]) for h in (0, 1) if use[i, h]])
    ph = np.nanmean([np.degrees(np.angle(T[i, h] / q8[i, h])) for h in (0, 1) if use[i, h]])
    print("%7.0f %5d | %7.1f %7.1f | %6.1f %6.2f | %7.1f %6.2f   phase vs q800 %+4.0f"
          % (depth[i], cal["nspk"][i], th[0], th[1], np.nanmean(np.abs(q8[i])), r8, np.nanmean(np.abs(q2[i])), r2, ph))
    rows.append(dict(depth=float(depth[i]), T=[float(x) for x in th], q800=float(np.nanmean(np.abs(q8[i]))),
                     q200=float(np.nanmean(np.abs(q2[i]))), ratio800=float(r8), ratio200=float(r2), phase800=float(ph)))

fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
dd = np.array([r["depth"] for r in rows])
for h, mk in ((0, "o"), (1, "s")):
    ax[0].scatter([r["T"][h] for r in rows], dd, marker=mk, color="#762a83", label="unit |T| half %d" % (h + 1))
ax[0].plot([r["q800"] for r in rows], dd, "-", color="0.4", label="DREDge 800 um")
ax[0].plot([r["q200"] for r in rows], dd, "--", color="#d6604d", label="DREDge 200 um")
ax[0].set_xlabel("cardiac-locked motion amplitude (um)")
ax[0].set_ylabel("depth (um)")
ax[0].legend(fontsize=8)
ax[0].set_title("msit004 cardiac motion along the probe")
ax[1].scatter([r["ratio800"] for r in rows], dd, color="0.3", label="unit / DREDge 800")
ax[1].scatter([r["ratio200"] for r in rows], dd, color="#d6604d", marker="x", label="unit / DREDge 200")
ax[1].axvline(1, color="k", lw=0.8)
ax[1].set_xlabel("ratio")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(HERE.parent / "card_by_depth.png", dpi=120)
(HERE / "card_by_depth.json").write_text(json.dumps(rows, indent=1))
print("-> card_by_depth.png / .json")
