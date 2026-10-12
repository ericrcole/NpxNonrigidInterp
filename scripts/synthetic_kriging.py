"""Ground-truth test of the kriging step on its own, independent of DREDge.

    python synthetic_kriging.py

A spike-like footprint exp(-distance / lambda) is sampled at the true NP1016 site
positions, the tissue is moved by s um in depth, and the pipeline's kriging operator
(same kernel, sig = 30 um, nugget 0.05) corrects it with shifts = -s. The corrected
footprint is compared with the unmoved truth. Repeated for the correct NP1016 channel map
(aligned columns x = {27, 59}, from settings.xml) and for chanMap_default.mat (staggered,
wrong for this probe on 192/384 sites), averaged over source positions and footprint widths.

Writes interp_eval/synthetic_kriging.png.
"""
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SIG = 30.0
NUG = 0.05
PROJ = Path(r"C:/Users/Data1/Documents/Projects/NPIX_MSIT")
import recinfo

OUT = recinfo.RESULTS


def kern(A, B):
    dx = np.abs(A[:, 0][:, None] - B[:, 0][None, :])
    dy = np.abs(A[:, 1][:, None] - B[:, 1][None, :])
    return np.exp(-dx / SIG - dy / (1.5 * SIG))


y = np.repeat(np.arange(192) * 20.0, 2)
P_true = np.c_[np.tile([27.0, 59.0], 192), y]
with h5py.File(PROJ / "chanMap_default.mat", "r") as f:
    P_wrong = np.c_[np.array(f["xcoords"]).ravel(), np.array(f["ycoords"]).ravel()]
MAPS = {"correct NP1016 map": P_true, "chanMap_default (wrong)": P_wrong}
IK = {k: np.linalg.inv(kern(P, P) + NUG * np.eye(len(P))) for k, P in MAPS.items()}

shifts = np.arange(0.0, 121.0, 2.5)
src_x = [27.0, 35.0, 43.0, 51.0, 59.0]
lams = [25.0, 40.0]
src_y = [1503.0, 1911.0, 2317.0]


def footprint(P, sx, sy, lam):
    return np.exp(-np.hypot(P[:, 0] - sx, P[:, 1] - sy) / lam)


def centroid(v, sy):
    m = np.abs(P_true[:, 1] - sy) < 160
    w = np.clip(v[m], 0, None)
    return float(np.sum(w * P_true[m, 1]) / max(w.sum(), 1e-12))


res = {k: dict(cen=[], peak=[], err=[]) for k in list(MAPS) + ["uncorrected"]}
for s in shifts:
    acc = {k: dict(cen=[], peak=[], err=[]) for k in res}
    ops = {}
    for k, P in MAPS.items():
        Q = P.copy()
        Q[:, 1] = Q[:, 1] + s          # shifts = -s  ->  yp = y - shifts = y + s
        ops[k] = kern(Q, P) @ IK[k]
    for sx in src_x:
        for sy in src_y:
            for lam in lams:
                truth = footprint(P_true, sx, sy, lam)
                moved = footprint(P_true, sx, sy + s, lam)
                pk = int(np.argmax(truth))
                ct = centroid(truth, sy)
                outs = {k: ops[k] @ moved for k in MAPS}
                outs["uncorrected"] = moved
                for k, o in outs.items():
                    acc[k]["cen"].append(centroid(o, sy) - ct)
                    acc[k]["peak"].append(o[pk] / truth[pk])
                    acc[k]["err"].append(np.linalg.norm(o - truth) / np.linalg.norm(truth))
    for k in res:
        for m in ("cen", "peak", "err"):
            res[k][m].append(np.mean(acc[k][m]))

colors = {"correct NP1016 map": "#1b7837", "chanMap_default (wrong)": "#d6604d", "uncorrected": "0.5"}
fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
for k, c in colors.items():
    ls = "--" if k == "uncorrected" else "-"
    axes[0].plot(shifts, res[k]["cen"], ls, color=c, lw=2, label=k)
    if k != "uncorrected":
        axes[1].plot(shifts, res[k]["peak"], ls, color=c, lw=2, label=k)
    axes[2].plot(shifts, res[k]["err"], ls, color=c, lw=2, label=k)
axes[0].set_ylabel("Residual depth error of footprint centroid (um)")
axes[0].set_title("Position is restored")
axes[1].set_ylabel("Amplitude on the true peak site (1 = preserved)")
axes[1].set_title("Amplitude cost (dips at half-pitch shifts)")
axes[1].set_ylim(0.5, 1.05)
axes[2].set_ylabel("Relative waveform error vs unmoved truth")
axes[2].set_title("Overall footprint error")
for ax in axes:
    ax.set_xlabel("Tissue displacement corrected (um)")
    ax.grid(alpha=0.3)
    for p in np.arange(10, 121, 20):
        ax.axvline(p, color="0.85", lw=0.8, zorder=0)
axes[0].legend(fontsize=9)
fig.suptitle("Kriging step alone, synthetic spike footprints on NP1016 geometry "
             "(grey lines = half-pitch shifts)", fontsize=11)
plt.tight_layout()
fig.savefig(OUT / "synthetic_kriging.png", dpi=140, bbox_inches="tight")

print("%-26s %14s %14s %14s" % ("", "max |centroid|", "min peak", "mean error"))
for k in res:
    print("%-26s %11.2f um %14.2f %14.3f" % (
        k, np.max(np.abs(res[k]["cen"])), np.min(res[k]["peak"]), np.mean(res[k]["err"])))
on = np.isin(shifts % 20, [0.0])
half = np.isin(shifts % 20, [10.0])
for k in MAPS:
    pk = np.array(res[k]["peak"])
    print("%-26s peak on-grid %.2f   at half-pitch %.2f" % (k, pk[on].mean(), pk[half].mean()))
print("figure -> %s" % (OUT / "synthetic_kriging.png"))
