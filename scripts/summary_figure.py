"""One-page summary: what the interpolation does, and where DREDge is or is not right.

    python summary_figure.py

Reads only caches written by the other scripts (no raw-data I/O):
  ../msit004_m4_final_tracking.npz          template_tracking.py on the final msit004 output
  ../_np150_tracking.npz                    template_tracking.py on msit001 (correct geometry)
  calib_msit004_driftfree_slow_vs_final.npz cardiac_calib.py measure, slow-only vs final (preferred)
  calib_msit004_m4_full.npz                 cardiac_calib.py measure, raw vs baseline (fallback)
  score_msit004_m4_full.json                cardiac_calib.py score (candidate table, raw-based)
  injection_realised.json                   injection_test.py analyze (tracker-free ground truth)
  calib_msit004_syn_card.npz                cardiac_calib.py measure on the synthetic segment
  ../msit004_syn_card_tracking.npz          template_tracking.py on the synthetic segment
Writes ../summary_interpolation.png

Row 1  slow drift: tracked positions; paired SD box plot; slope decomposition separating
       interpolation fidelity (after-slope minus before-slope = -1) from DREDge accuracy.
Row 2  cardiac: measured phasor change vs applied shift (fidelity); units' calibrated motion vs
       DREDge's amplitude (accuracy); predicted residual per candidate.
Row 3  observed cardiac amplitude before vs after the final correction; ground-truth
       injection test; settled facts.
"""
import json
import os
from pathlib import Path

import matplotlib
import numpy as np
from scipy.stats import wilcoxon

matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.environ.setdefault("NPX_REC", "msit004")
import cardiac_calib as cc  # noqa: E402
import recinfo  # noqa: E402

HERE = recinfo.RESULTS
OUT = recinfo.RESULTS
CARD_AMP_SYN = 50.0
FINAL_CARD_GAIN = 0.6


def paired_p(a, b):
    try:
        return wilcoxon(a, b, alternative="greater").pvalue
    except ValueError:
        return np.nan


def missing(ax, what):
    ax.axis("off")
    ax.text(0.5, 0.5, "%s\nnot available yet" % what, ha="center", va="center", color="0.4")


fig = plt.figure(figsize=(20, 17))
gs = fig.add_gridspec(3, 3, hspace=0.5, wspace=0.32)

# ============================ row 1: slow drift ============================
# Prefer the like-for-like pair (raw vs slow-only): both files keep cardiac/fast motion, so unit
# templates are equally smeared and SD / slope comparisons are not biased by differing tracker gain.
# A/B show raw vs the calibrated v2 pipeline (slow x0.49 + cardiac x0.28); C keeps the slow-only pair,
# whose slopes are regressed on the same displacement before and after.
for cand in ("msit004_m4_v3_tracking.npz", "msit004_m4_v2_tracking.npz", "msit004_m4_final_tracking.npz"):
    t4 = OUT / cand
    if t4.exists():
        break
t4c = OUT / "msit004_m4_slowonly_tracking.npz"
t1 = OUT / "_np150_tracking.npz"
if t4.exists():
    d4 = np.load(t4)
    sub = gs[0, 0].subgridspec(2, 1, hspace=0.08)
    lim = np.nanmax(np.abs(np.concatenate([d4["disp_before"], d4["disp_after"]])))
    for k, (key, lab) in enumerate([("disp_before", "before"), ("disp_after", "after")]):
        ax = fig.add_subplot(sub[k])
        for row in d4[key]:
            ax.plot(d4["ctr"], row, lw=1.0, alpha=0.85)
        ax.axhline(0, color="k", ls=":", lw=0.8)
        ax.set_ylim(-lim * 1.05, lim * 1.05)
        ax.text(0.01, 0.88, lab.upper(), transform=ax.transAxes, fontsize=10, fontweight="bold")
        ax.set_ylabel("position (um)")
        if k == 0:
            ax.set_title("A  msit004 raw vs %s: best-correlated position of each unit" % "_".join(t4.name.split("_")[1:3]))
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Time (s)")
    ax = fig.add_subplot(gs[0, 1])
    sb, sa = d4["sd_before"], d4["sd_after"]
    ok = np.isfinite(sb) & np.isfinite(sa)
    bp = ax.boxplot([sb[ok], sa[ok]], labels=["before", "after"], widths=0.55,
                    patch_artist=True, showfliers=False)
    for patch, col in zip(bp["boxes"], ["#9ecae1", "#fdae6b"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    for i in np.where(ok)[0]:
        ax.plot([1, 2], [sb[i], sa[i]], "-o", color="0.35", ms=4, lw=0.9, alpha=0.8)
    ax.set_ylabel("SD of tracked position (um)")
    ax.set_title("B  msit004 slow drift: per-unit variability, paired\n"
                 "median %.1f -> %.1f um, %d/%d improve, Wilcoxon p = %.3g"
                 % (np.median(sb[ok]), np.median(sa[ok]), int((sb[ok] > sa[ok]).sum()), ok.sum(),
                    paired_p(sb[ok], sa[ok])), fontsize=10)
else:
    missing(fig.add_subplot(gs[0, 0]), t4.name)
    missing(fig.add_subplot(gs[0, 1]), t4.name)

ax = fig.add_subplot(gs[0, 2])
xs = np.linspace(-2.2, 2.2, 10)
ax.plot(xs, xs - 1, "--", color="0.5", lw=1.2, label="after = before - 1  (interpolation exact)")
ax.axhline(0, color="#1b7837", lw=1.2, label="after = 0  (motion fully removed)")
ax.axvline(1, color="0.8", lw=0.8)
notes = []
for path, lab, mk, col in [(t4c, "msit004", "o", "#2166ac"), (t1, "msit001", "s", "#b2182b")]:
    if not path.exists():
        continue
    d = np.load(path)
    m = np.isfinite(d["slope_before"]) & np.isfinite(d["slope_after"])
    if not m.any():
        continue
    ax.scatter(d["slope_before"][m], d["slope_after"][m], marker=mk, s=45,
               facecolor=col if lab == "msit004" else "none", edgecolor=col, label=lab, zorder=3)
    dl = d["slope_after"][m] - d["slope_before"][m]
    notes.append("%s: after-before %.2f +/- %.2f (n=%d)" % (lab, dl.mean(), dl.std(), m.sum()))
ax.set_xlabel("slope BEFORE (tracked position vs applied displacement)\n1 = unit moved as DREDge says")
ax.set_ylabel("slope AFTER")
ax.set_title("C  Slow: interpolation fidelity vs DREDge accuracy\n" + "\n".join(notes), fontsize=9)
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.3)

# ============================ row 2: cardiac ============================
drift_path = HERE / "calib_msit004_driftfree_slow_vs_final.npz"
base_path = HERE / "calib_msit004_m4_full.npz"
if drift_path.exists():
    cal, scale = np.load(drift_path, allow_pickle=True), FINAL_CARD_GAIN
    design = "drift-free pair: slow-only vs final (x%.1f cardiac)" % FINAL_CARD_GAIN
elif base_path.exists():
    cal, scale = np.load(base_path, allow_pickle=True), 1.0
    design = "raw vs baseline run - template-smear caveat, see README"
else:
    cal = None
g_used = None
if cal is not None:
    depth, zB, zA, okc = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
    T0, Tmid, T1 = cal["T"]
    halves = [(T0, Tmid), (Tmid, T1)]
    q0 = cc.cardiac_phasors(str(cal["baseline"]), depth, halves, cal["mt"], cal["phase_u"])
    use = okc & np.isfinite(q0) & np.isfinite(zA) & np.isfinite(zB)
    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = use[:, h]
        qa = scale * q0[u, h]
        g[h] = np.sum(-(zA[u, h] - zB[u, h]) * np.conj(qa)) / np.sum(np.abs(qa) ** 2)
    meas, pred, T, Q, Bv, Av = [], [], [], [], [], []
    for h in (0, 1):
        gc = g[1 - h]
        u = use[:, h]
        meas.append((zA[u, h] - zB[u, h]) / abs(gc))
        pred.append(-scale * q0[u, h] * gc / abs(gc))
        T.append(zB[u, h] / gc)
        Q.append(q0[u, h])
        Bv.append(np.abs(zB[u, h]) / abs(gc))
        Av.append(np.abs(zA[u, h]) / abs(gc))
    meas, pred, T, Q, Bv, Av = (np.concatenate(v) for v in (meas, pred, T, Q, Bv, Av))
    g_used = g

    ax = fig.add_subplot(gs[1, 0])
    ang = np.degrees(np.angle(meas / pred))
    sc = ax.scatter(np.abs(pred), np.abs(meas), c=ang, cmap="coolwarm", vmin=-90, vmax=90, s=45, zorder=3)
    mx = 1.1 * max(np.abs(pred).max(), np.abs(meas).max())
    ax.plot([0, mx], [0, mx], "--", color="0.5")
    fig.colorbar(sc, ax=ax, label="phase error (deg)")
    ax.set_xlabel("cardiac shift APPLIED at unit depth (um)")
    ax.set_ylabel("measured change |A-B|/|g| (um)")
    ax.set_title("D  Cardiac: interpolation fidelity on real units\n%s\n"
                 "|measured|/|applied| median %.2f, |phase error| median %.0f deg, n=%d"
                 % (design, np.median(np.abs(meas) / np.abs(pred)), np.median(np.abs(ang)), len(meas)), fontsize=8.5)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 1])
    ax.scatter(np.abs(Q), np.abs(T), s=45, color="#762a83", zorder=3)
    mx = 1.1 * max(np.abs(Q).max(), np.abs(T).max())
    ax.plot([0, mx], [0, mx], "--", color="0.5", label="DREDge amplitude correct")
    k = np.sum(T * np.conj(Q)) / np.sum(np.abs(Q) ** 2)
    ax.plot([0, mx], [0, abs(k) * mx], "-", color="#762a83", label="fitted |k| = %.2f" % abs(k))
    ax.set_xlabel("DREDge cardiac amplitude at unit depth, |q0| (um)")
    ax.set_ylabel("unit's own cardiac motion, |B/g| (um)")
    ax.set_title("E  Cardiac: DREDge accuracy\n"
                 "tracker gain |g| = %.2f / %.2f (halves); units move %.2fx DREDge at %+.0f deg"
                 % (abs(g[0]), abs(g[1]), abs(k), np.degrees(np.angle(k))), fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[2, 0])
    bpl = ax.boxplot([Bv, Av], labels=["before", "after"], widths=0.55, patch_artist=True, showfliers=False)
    for patch, col in zip(bpl["boxes"], ["#9ecae1", "#fdae6b"]):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    for i in range(len(Bv)):
        ax.plot([1, 2], [Bv[i], Av[i]], "-o", color="0.35", ms=4, lw=0.9, alpha=0.8)
    ax.set_ylabel("cardiac-locked motion of unit (um, tracker-calibrated)")
    ax.set_title("G  Cardiac band, observed before vs after correction\n%s\n"
                 "power ratio after/before %.2f; %d/%d reduced; Wilcoxon p = %.3g"
                 % (design, np.sum(Av ** 2) / np.sum(Bv ** 2), int((Av < Bv).sum()), len(Bv), paired_p(Bv, Av)),
                 fontsize=8.5)
else:
    for pos in (gs[1, 0], gs[1, 1], gs[2, 0]):
        missing(fig.add_subplot(pos), "cardiac calibration cache")

ax = fig.add_subplot(gs[1, 2])
score_path = HERE / "score_msit004_driftfree.json"
if score_path.exists():
    rows = sorted(json.loads(score_path.read_text())["candidates"], key=lambda r: r["RES"])
    names = [r["tag"].replace("msit004_", "").replace("_t805-1745", "") for r in rows]
    res = [min(r["RES"], 3.0) for r in rows]
    # green = meaningful help, grey = essentially none, red = would add motion
    ax.barh(range(len(res)), res, color=["#1b7837" if v < 0.9 else ("0.65" if v < 1.0 else "#d6604d") for v in res])
    ax.axvline(1, color="k", lw=1)
    ax.set_yticks(range(len(res)))
    ax.set_yticklabels(names, fontsize=8)
    ax.invert_yaxis()
    for i, r in enumerate(rows):
        ax.text(res[i] + 0.03, i, "units %.2fx at %+.0f deg; best gain -> %.2f" % (r["k"], r["k_deg"], r["RES_opt"]),
                va="center", fontsize=7)
    ax.set_xlim(0, 4.6)
    ax.set_xlabel("residual cardiac power if applied as-is (<1 helps; bars clipped at 3)")
    ax.set_title("F  Candidate cardiac estimates vs units' drift-free cardiac motion\n"
                 "(tracker gain fitted on the other half)", fontsize=9)
else:
    missing(ax, score_path.name)

# ============================ row 3: ground truth + facts ============================
ax = fig.add_subplot(gs[2, 1])
inj = HERE / "injection_realised.json"
syn_cal = HERE / "calib_msit004_syn_card.npz"
syn_tt = OUT / "msit004_syn_card_tracking.npz"
if inj.exists():
    real = json.loads(inj.read_text())["after"]
    labels = ["pipeline\nslow", "pipeline\ncardiac"]
    vals = [real["slow"], real["card_abs"]]
    cols = ["#1b7837", "#1b7837"]
    if syn_cal.exists():
        zs = np.load(syn_cal, allow_pickle=True)["zA"]
        zs = zs[np.isfinite(zs)]
        gsyn = abs(complex(-np.median(zs.real), -np.median(zs.imag)) / CARD_AMP_SYN)
        labels.append("cardiac\ntracker")
        vals.append(gsyn / max(real["card_abs"], 1e-9))
        cols.append("#762a83")
    if syn_tt.exists():
        sa_ = np.load(syn_tt)["slope_after"]
        sa_ = sa_[np.isfinite(sa_)]
        labels.append("slow\ntracker")
        vals.append(-np.median(sa_) / max(real["slow"], 1e-9))
        cols.append("#762a83")
    ax.bar(range(len(vals)), vals, color=cols)
    ax.axhline(1, color="k", lw=1)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.02, "%.2f" % v, ha="center", fontsize=9)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0, max(1.25, max(vals) + 0.15))
    ax.set_ylabel("gain (1 = exact)")
    ax.set_title("H  Ground truth: static synthetic units + known motion through\n"
                 "the real MATLAB pipeline (green) and each tracker (purple)", fontsize=9)
else:
    missing(ax, inj.name)

ax = fig.add_subplot(gs[2, 2])
ax.axis("off")
facts = [
    "I  Settled facts (details in README.md)",
    "",
    "- Real-data realisation (chunk registration): the pipeline applied",
    "  slow x1.01 and cardiac x0.96 of the specified shift (robust fit).",
    "- DREDge-LFP amplitude vs units (msit004, held out by half):",
    "  slow drift 0.49 / 0.48x, cardiac 0.30 / 0.26x (timing within 12 deg).",
    "  At face value slow correction overshoots (units end up moving",
    "  -0.5x DREDge) and cardiac correction triples cardiac motion.",
    "- DREDge-LFP recovers injected LFP shifts at 1.00x: real LFP-vs-unit gap.",
    "- v2 (slow x0.49, cardiac x0.28): slow SD 23.5 -> 14.8 um, 5/5 units;",
    "  slope bracket puts the best slow gain at ~0.6. Cardiac residual",
    "  1.19-1.35 (p 0.71): no cardiac gain helped on msit004.",
    "- Probe NP1016 has ALIGNED columns (settings.xml); chanMap_default.mat",
    "  is staggered and wrong on 192/384 sites. Use chanMap_np1016.mat.",
    "- Kriging, synthetic footprints: centroid error <= 0.06 um for any shift;",
    "  amplitude 0.96 on-grid, 0.84 at half-pitch.",
    "- DC-offset mixing removed: <150 Hz std 108 -> 29 (raw 30);",
    "  spike band unchanged (0.9% RMS).",
    "- msit004 0.1-0.5 Hz motion sits on the MSIT trial rate (0.221 Hz):",
    "  not applied. 2-5 Hz is the cardiac 2nd harmonic (PLV 0.96).",
    "- LFP high-pass before DREDge does NOT fix cardiac: phase error",
    "  +72 deg and slow-drift agreement 0.72 -> 0.18.",
    "- DREDge temporal prior (lambda_t 10, 100) leaves cardiac unchanged.",
    "- estimate_motion() overrides dredge_lfp's 800 um windows with 200 um.",
]
ax.text(0.0, 1.0, "\n".join(facts), va="top", ha="left", fontsize=8.5, family="monospace")

fig.suptitle("DREDge + kriging on NP1016 recordings: the interpolation applies what it is given; "
             "whether that helps depends on DREDge", fontsize=14)
fig.savefig(OUT / "summary_interpolation.png", dpi=120, bbox_inches="tight")
print("figure -> %s" % (OUT / "summary_interpolation.png"))
