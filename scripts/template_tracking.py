"""Waveform-template spatial cross-correlation tracking, before vs after interpolation.

    python template_tracking.py <segdir> <t0_s> <t1_s> [bin_s]

Produces, for a prepared segment:
  <seg>_examples.png   raw traces with detected spikes, mean waveform, and the
                       space x time cross-correlation map for example units
  <seg>_tracking.png   best-correlated channel position vs time (all units), and a
                       box plot of its variability before vs after
  <seg>_tracking.txt   per-unit numbers

Design
------
Spike TIMES are identical in the before and after files - the interpolation resamples
across channels and never touches timing. So the same event train is used for both
files, giving a PAIRED comparison of the same unit rather than independently selected
units that may not correspond. Each file builds its own template from the same reference
window, so displacement is measured relative to each file's own baseline; we are
comparing stability, not absolute position.
"""
import io
import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
from scipy.signal import butter, filtfilt

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import recinfo
_R = recinfo.get()
AP = _R["ap"]
OUT = recinfo.RESULTS

NCH = _R["nch"]
NCHN = _R["nchan_neural"]
FS = 30000.0
THR = 6.0
MERGE_MS = 0.5
MERGE_CH = 10        # rows (=200 um)

TMPL_CH = 6          # template half-height, ROWS (=120 um)
SEARCH_CH = 25       # +/- search range, ROWS (=500 um)
WF_PRE, WF_POST = 20, 40
MAX_SAMP_LAG = 3     # tolerate small timing jitter
MIN_SPK_BIN = 4
MIN_PEAK_R = 0.7     # reject ambiguous per-bin matches (good bins run 0.95+)
RESTRICT_UM = 250.0  # peak must lie within this of the unit's own median position
EDGE_ROWS = 3        # reject a bin whose NCC peak sits this close to the search edge
MAX_EDGE_FRAC = 0.15 # reject a unit hitting the edge more often than this
MIN_P2P_UM = float(os.environ.get("NPX_MIN_P2P_UM", "15.0"))  # before-side gate fixed a priori (15 um); override only for synthetic truth tests
N_UNITS = 20
N_CAND = 40          # screen this many, keep N_UNITS by bin coverage
MIN_RATE = 0.3       # Hz
MIN_SNR = 4.0
MIN_SEP_CH = 5       # rows between distinct units (=100 um)

seg = sys.argv[1]
T0, T1 = float(sys.argv[2]), float(sys.argv[3])
BIN = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0
S = AP / seg
SEG_T0 = float((S / "segment_t0.txt").read_text().strip())
buf = io.StringIO()
PRE = recinfo.key() + "_"   # output files are named <recording>_<segment>_*


def say(m):
    print(m)
    buf.write(m + "\n")


yc_ch, xc_ch = recinfo.chanmap_yx()
yc, RMAT = recinfo.row_matrix()   # track in row space, not channel index
NROW = len(yc)
ROW_OF = np.searchsorted(yc, yc_ch[:_R["nchan_neural"]])   # channel index -> row index
UM_PER_CH = float(np.median(np.diff(yc)))   # um per ROW (uniform pitch)
B = S / "continuous.dat"
A = S / ("interp" + seg) / ("interp" + seg + ".dat")
nfr = B.stat().st_size // (NCH * 2)
mm = {"before": np.memmap(B, dtype="int16", mode="r", shape=(nfr, NCH)),
      "after": np.memmap(A, dtype="int16", mode="r", shape=(nfr, NCH))}
bp = butter(3, [300 / (FS / 2), 6000 / (FS / 2)], btype="band")
say("segment %s  t=%.0f-%.0f s  bin=%.2f s  %.2f um/channel" % (seg, T0, T1, BIN, UM_PER_CH))


def load(which, s0, n):
    x = np.asarray(mm[which][s0:s0 + n, :NCHN]).astype(np.float32)
    x = filtfilt(*bp, x, axis=0)
    return x - np.median(x, axis=1, keepdims=True)


F0, F1 = int((T0 - SEG_T0) * FS), int((T1 - SEG_T0) * FS)
BLK = int(10 * FS)


def detect_before():
    """Merged threshold crossings on the BEFORE file; times reused for both files."""
    sd = sdc = None
    TT, CC, AA = [], [], []
    for s in range(F0, F1, BLK):
        e = min(F1, s + BLK)
        x = load("before", s, e - s)
        if sdc is None:
            sdc = np.median(np.abs(x), axis=0) / 0.6745
            sdc[sdc < 1e-6] = np.inf
            sd = sdc @ RMAT           # per-row noise
        mid = x[1:-1]
        mask = (mid < -THR * sdc) & (mid < x[:-2]) & (mid < x[2:])
        ii, ch = np.nonzero(mask)
        ii = ii + 1 + s
        amp = x[ii - s, ch]
        ch = ROW_OF[ch]                       # collapse to rows
        o = np.lexsort((ch, ii))
        ii, ch, amp = ii[o], ch[o], amp[o]
        wnd = int(MERGE_MS * 1e-3 * FS)
        keep = np.ones(len(ii), bool)
        for j in range(len(ii)):
            if not keep[j]:
                continue
            g = np.nonzero((ii > ii[j] - wnd) & (ii < ii[j] + wnd) & (np.abs(ch - ch[j]) <= MERGE_CH))[0]
            if len(g) > 1:
                bst = g[np.argmin(amp[g])]
                keep[g] = False
                keep[bst] = True
        TT.append(ii[keep])
        CC.append(ch[keep])
        AA.append(amp[keep])
        del x
    return np.concatenate(TT), np.concatenate(CC), np.concatenate(AA), sd


t_all, c_all, a_all, chan_sd = detect_before()
dur = T1 - T0
say("detected %d merged events on the before file" % len(t_all))

# ---- pick example units: high rate, high SNR, well separated in depth ----
lim = TMPL_CH + SEARCH_CH
cand = []
for c in range(lim, NROW - lim):
    k = np.abs(c_all - c) <= 3
    n = int(k.sum())
    if n / dur < MIN_RATE:
        continue
    snr = float(np.median(-a_all[k]) / chan_sd[c])
    if snr < MIN_SNR:
        continue
    cand.append(dict(chan=c, n=n, rate=n / dur, snr=snr, t=t_all[k]))
cand.sort(key=lambda d: -d["snr"] * np.sqrt(d["n"]))
units = []
for d in cand:
    if all(abs(d["chan"] - u["chan"]) >= MIN_SEP_CH for u in units):
        units.append(d)
    if len(units) >= N_CAND:
        break
units.sort(key=lambda d: d["chan"])
say("\n%d example units:" % len(units))
say("  %5s %9s %8s %7s %8s" % ("chan", "depth_um", "rate_Hz", "SNR", "n_spikes"))
for u in units:
    say("  %5d %9.0f %8.1f %7.1f %8d" % (u["chan"], yc[u["chan"]], u["rate"], u["snr"], u["n"]))
if not units:
    raise SystemExit("no units passed selection")

# ---- accumulate per-bin mean waveforms for each unit, in each file ----
edges = np.arange(T0, T1 + BIN, BIN)
ctr = edges[:-1] + BIN / 2
nb = len(ctr)
nband = 2 * (TMPL_CH + SEARCH_CH) + 1
nsamp = WF_PRE + WF_POST + 1
ACC = {w: np.zeros((len(units), nb, nband, nsamp), np.float32) for w in mm}
CNT = np.zeros((len(units), nb), np.int64)

for w in mm:
    for s in range(F0, F1, BLK):
        e = min(F1, s + BLK)
        a0 = max(0, s - WF_PRE - 5)
        a1 = min(nfr, e + WF_POST + 5)
        x = load(w, a0, a1 - a0) @ RMAT
        for ui, u in enumerate(units):
            band = np.arange(u["chan"] - lim, u["chan"] + lim + 1)
            st = u["t"][(u["t"] >= s) & (u["t"] < e)]
            st = st[(st - WF_PRE >= a0) & (st + WF_POST < a1)]
            if not len(st):
                continue
            idx = (st - a0)[:, None] + np.arange(-WF_PRE, WF_POST + 1)[None, :]
            sn = x[:, band][idx].transpose(0, 2, 1)          # (nspk, nband, nsamp)
            bi = np.clip(((st / FS + SEG_T0 - T0) / BIN).astype(int), 0, nb - 1)
            for b in np.unique(bi):
                m = bi == b
                ACC[w][ui, b] += sn[m].sum(0)
                if w == "before":
                    CNT[ui, b] += int(m.sum())
        del x
cov = (CNT >= MIN_SPK_BIN).mean(axis=1)
keep = np.argsort(-cov)[:N_UNITS]
keep = np.array(sorted(k for k in keep if cov[k] >= 0.25))
if len(keep) == 0:
    keep = np.argsort(-cov)[:min(3, len(units))]
    keep = np.array(sorted(keep))
say("bin coverage per candidate: " + " ".join("ch%d=%.0f%%" % (units[i]["chan"], 100 * cov[i])
                                              for i in np.argsort(-cov)[:10]))
say("keeping %d units with >=25%% of bins filled" % len(keep))
units = [units[i] for i in keep]
CNT = CNT[keep]
ACC = {w: ACC[w][keep] for w in ACC}
MEAN = {w: ACC[w] / np.maximum(CNT, 1)[:, :, None, None] for w in mm}
say("\nmedian spikes per %.1f s bin: %d" % (BIN, int(np.median(CNT[CNT > 0]))))


def ncc_profile(T, W):
    """Normalised cross-correlation of template T against window W over depth lags."""
    nt, ns = T.shape
    core = slice(MAX_SAMP_LAG, ns - MAX_SAMP_LAG)
    Tc = T[:, core]
    Tz = Tc - Tc.mean()
    Tn = np.linalg.norm(Tz)
    nl = W.shape[0] - nt + 1
    out = np.full(nl, np.nan)
    if Tn == 0:
        return out
    for i in range(nl):
        best = -np.inf
        for sft in range(-MAX_SAMP_LAG, MAX_SAMP_LAG + 1):
            P = W[i:i + nt, MAX_SAMP_LAG + sft:ns - MAX_SAMP_LAG + sft]
            Pz = P - P.mean()
            pn = np.linalg.norm(Pz)
            if pn > 0:
                best = max(best, float((Pz * Tz).sum() / (pn * Tn)))
        out[i] = best
    return out


def parabolic(y):
    i = int(np.nanargmax(y))
    if 0 < i < len(y) - 1 and np.all(np.isfinite(y[i - 1:i + 2])):
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        if abs(den) > 1e-12:
            return i + 0.5 * (a - c) / den, float(b)
    return float(i), float(y[i])


tsl = slice(SEARCH_CH, SEARCH_CH + 2 * TMPL_CH + 1)
nlag = nband - (2 * TMPL_CH + 1) + 1
lag_um = (np.arange(nlag) - SEARCH_CH) * UM_PER_CH
K_REF = max(1, int(round(10.0 / BIN)))       # 10 s reference window

PROF = {w: np.full((len(units), nb, nlag), np.nan) for w in mm}
DISP = {w: np.full((len(units), nb), np.nan) for w in mm}
PEAKR = {w: np.full((len(units), nb), np.nan) for w in mm}
REFB = np.zeros(len(units), int)
EDGE = {w: np.zeros((len(units), nb)) for w in mm}

for ui in range(len(units)):
    cs = np.convolve(CNT[ui], np.ones(K_REF), "valid")
    b0 = int(np.argmax(cs))
    REFB[ui] = b0
    sl = slice(b0, b0 + K_REF)
    wt = CNT[ui][sl].astype(float)
    if wt.sum() == 0:
        continue
    for w in mm:
        T = np.tensordot(wt, MEAN[w][ui, sl], axes=(0, 0)) / wt.sum()
        T = T[tsl]
        # pass 1: unconstrained profiles, to locate this unit's typical position
        profs = {}
        for b in range(nb):
            if CNT[ui, b] < MIN_SPK_BIN:
                continue
            pr = ncc_profile(T, MEAN[w][ui, b])
            PROF[w][ui, b] = pr
            profs[b] = pr
        if not profs:
            continue
        c0 = int(np.median([int(np.nanargmax(pr)) for pr in profs.values()]))
        # pass 2: the NCC ridge is broad (~200 um), so an unconstrained argmax can jump
        # between local maxima and fabricate huge excursions. Restrict it to a window
        # around the unit's own median position - applied identically to both files.
        half = int(round(RESTRICT_UM / UM_PER_CH))
        lo, hi = max(0, c0 - half), min(nlag, c0 + half + 1)
        for b, pr in profs.items():
            i, v = parabolic(pr[lo:hi])
            i += lo
            if v < MIN_PEAK_R:
                continue
            if i < EDGE_ROWS or i > nlag - 1 - EDGE_ROWS:
                EDGE[w][ui, b] = 1.0     # peak pinned at the search boundary
                continue
            DISP[w][ui, b] = (i - SEARCH_CH) * UM_PER_CH
            PEAKR[w][ui, b] = v

# ---- reject units the tracker cannot follow honestly ----
efrac = {w: EDGE[w].mean(axis=1) for w in mm}
p2p_b = np.nanmax(DISP["before"], axis=1) - np.nanmin(DISP["before"], axis=1)
good = np.ones(len(units), bool)
for ui in range(len(units)):
    why = []
    if max(efrac["before"][ui], efrac["after"][ui]) > MAX_EDGE_FRAC:
        why.append("edge-pinned %.0f%%/%.0f%%" % (100*efrac["before"][ui], 100*efrac["after"][ui]))
    if not np.isfinite(p2p_b[ui]) or p2p_b[ui] < MIN_P2P_UM:
        why.append("before-trace static (p2p %.1f um) - degenerate NCC peak" % p2p_b[ui])
    if why:
        good[ui] = False
        say("  DROP ch%d (%.0f um): %s" % (units[ui]["chan"], yc[units[ui]["chan"]], "; ".join(why)))
if good.sum() == 0:
    raise SystemExit("no units survived tracking quality gates")
units = [u for u, g in zip(units, good) if g]
for w in mm:
    DISP[w] = DISP[w][good]
    PROF[w] = PROF[w][good]
    PEAKR[w] = PEAKR[w][good]
CNT = CNT[good]
MEAN = {w: MEAN[w][good] for w in MEAN}
REFB = REFB[good]
say("%d units passed tracking quality gates" % len(units))
DC = {w: DISP[w] - np.nanmedian(DISP[w], axis=1, keepdims=True) for w in mm}
say("\ntracked bins: " + "  ".join("%s %d/%d (%.0f%%)" % (
    w, np.isfinite(DISP[w]).sum(), DISP[w].size, 100 * np.isfinite(DISP[w]).mean()) for w in mm))

# ================== FIGURE 1: example units ==================
show = list(range(min(3, len(units))))
fig, axes = plt.subplots(len(show), 4, figsize=(21, 4.6 * len(show)), squeeze=False)
for r, ui in enumerate(show):
    u = units[ui]
    b0 = REFB[ui]
    nwin = int(0.25 * FS)
    cw = np.arange(F0, max(F0 + 1, F1 - nwin), nwin)
    cn = np.histogram(u["t"], bins=np.append(cw, cw[-1] + nwin))[0]
    s0 = int(cw[int(np.argmax(cn))])
    vis = np.arange(max(0, u["chan"] - 12), min(NROW, u["chan"] + 13))
    for ci, w in enumerate(["before", "after"]):
        ax = axes[r][ci]
        x = (load(w, s0, nwin) @ RMAT)[:, vis]
        sp = 9 * np.std(x)
        tt = (s0 + np.arange(nwin)) / FS + SEG_T0
        for i in range(len(vis)):
            ax.plot(tt, x[:, i] + i * sp, color="k", lw=0.5)
        ev = u["t"][(u["t"] >= s0) & (u["t"] < s0 + nwin)]
        if len(ev):
            rel = ev - s0
            pk = np.argmin(np.stack([x[max(0, q - 2):q + 3].min(0) for q in rel]), axis=1)
            ax.scatter(tt[rel], x[rel, pk] + pk * sp, color="tab:red", s=26, zorder=5)
        ax.set_yticks(np.arange(len(vis))[::4] * sp)
        ax.set_yticklabels(vis[::4], fontsize=7)
        ax.set_title("%s  unit ch%d  (%d spikes shown)" % (w, u["chan"], len(ev)), fontsize=10)
        ax.set_xlabel("Time (s)")
        if ci == 0:
            ax.set_ylabel("Channel")
    for ci, w in enumerate(["before", "after"]):
        ax = axes[r][2 + ci]
        pr = PROF[w][ui]
        im = ax.imshow(pr.T, aspect="auto", origin="lower", cmap="magma", vmin=0, vmax=1,
                       extent=[T0, T1, lag_um[0], lag_um[-1]])
        ax.plot(ctr, DISP[w][ui], ".", color="#66ff66", ms=4)
        ax.axhline(0, color="w", ls=":", lw=0.8, alpha=0.7)
        ax.axvspan(ctr[b0], ctr[min(b0 + K_REF - 1, nb - 1)], color="cyan", alpha=0.12)
        sd = np.nanstd(DC[w][ui])
        ax.set_title("%s  NCC(depth, time)  SD=%.1f um" % (w, sd), fontsize=10)
        ax.set_xlabel("Time (s)")
        if ci == 0:
            ax.set_ylabel("Depth lag from template (um)")
        fig.colorbar(im, ax=ax, shrink=0.85, label="norm. xcorr")
plt.tight_layout()
fig.savefig(OUT / (PRE + seg + "_examples.png"), dpi=135, bbox_inches="tight")

# ================== FIGURE 2: tracking + box plot ==================
fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
for ax, w in zip(axes[:2], ["before", "after"]):
    for ui in range(len(units)):
        ax.plot(ctr, DC[w][ui], lw=1.3, label="ch%d" % units[ui]["chan"])
    ax.axhline(0, color="k", ls=":", lw=0.8)
    ax.set_title("%s - best-correlated channel position" % w.upper())
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Displacement from template (um)")
    ax.legend(fontsize=7, ncol=2)
lim_y = np.nanmax(np.abs(np.concatenate([DC["before"], DC["after"]])))
for ax in axes[:2]:
    ax.set_ylim(-lim_y * 1.05, lim_y * 1.05)

sd = {w: np.nanstd(DC[w], axis=1) for w in mm}
ok = np.isfinite(sd["before"]) & np.isfinite(sd["after"])
ax = axes[2]
bpl = ax.boxplot([sd["before"][ok], sd["after"][ok]], labels=["before", "after"],
                 widths=0.55, patch_artist=True, showfliers=False)
for p, col in zip(bpl["boxes"], ["#9ecae1", "#fdae6b"]):
    p.set_facecolor(col)
    p.set_alpha(0.75)
for i in range(ok.sum()):
    ax.plot([1, 2], [sd["before"][ok][i], sd["after"][ok][i]], "-o", color="0.35", ms=4, lw=0.9, alpha=0.8)
ax.set_ylabel("SD of tracked position (um)")
ax.set_title("Per-unit variability (paired, n=%d)" % ok.sum())
plt.tight_layout()
fig.savefig(OUT / (PRE + seg + "_tracking.png"), dpi=140, bbox_inches="tight")

say("\n%5s %9s %10s %10s %8s %8s %8s" % ("chan", "depth_um", "SD_before", "SD_after", "ratio", "p2p_bef", "p2p_aft"))
for ui in np.where(ok)[0]:
    p2b = np.nanmax(DC["before"][ui]) - np.nanmin(DC["before"][ui])
    p2a = np.nanmax(DC["after"][ui]) - np.nanmin(DC["after"][ui])
    say("%5d %9.0f %10.1f %10.1f %8.2f %8.1f %8.1f" % (
        units[ui]["chan"], yc[units[ui]["chan"]],
        sd["before"][ui], sd["after"][ui], sd["before"][ui] / max(sd["after"][ui], 1e-9), p2b, p2a))
if ok.sum():
    from scipy.stats import wilcoxon
    mb, ma = np.median(sd["before"][ok]), np.median(sd["after"][ok])
    say("\nmedian SD  %.1f -> %.1f um   (%.2fx reduction)" % (mb, ma, mb / max(ma, 1e-9)))
    if ok.sum() >= 5:
        try:
            st, p = wilcoxon(sd["before"][ok], sd["after"][ok], alternative="greater")
            say("Wilcoxon signed-rank (paired, before > after): W=%.1f, p=%.4g" % (st, p))
        except ValueError as exc:
            say("Wilcoxon not computed: %s" % exc)
    say("median peak correlation: before %.2f, after %.2f"
        % (np.nanmedian(PEAKR["before"]), np.nanmedian(PEAKR["after"])))

# ---- does each unit's tracked position follow the displacement that was APPLIED? ----
# Regress tracked position on the applied displacement at the unit's depth, binned like the
# tracker. Before correction the slope should be ~1 (unit moves with the tissue as DREDge
# says) and after ~0 (that motion removed). Unlike SD this is immune to tracker noise, which
# adds variance but cannot correlate with the applied shift. The applied displacement is
# reproduced as applyDredge does it: resample to the batch rate, then the same low-pass.
from scipy.interpolate import PchipInterpolator as _Pchip
from scipy.signal import resample_poly as _rp
_D = np.loadtxt(S / "nonrigid_displacement.csv", delimiter=",")
_D = _D[None, :] if _D.ndim == 1 else _D
_cc = np.atleast_1d(np.loadtxt(S / "dredge_window_centers.csv", delimiter=","))
_fsb = FS / 128.0
_h = _rp(_D.T, 120, 128, axis=0)
_lp = _R.get("motion_lowpass_hz")
if _lp:
    _h = filtfilt(*butter(3, _lp / (_fsb / 2), btype="low"), _h, axis=0)
_bt = SEG_T0 + np.arange(_h.shape[0]) / _fsb
_bi = np.digitize(_bt, edges) - 1
_inb = (_bi >= 0) & (_bi < nb)
REG = {}
for ui, u in enumerate(units):
    dep = yc[u["chan"]]
    hu = _Pchip(_cc, _h, axis=1)(np.clip(dep, _cc[0], _cc[-1])) if _h.shape[1] > 1 else _h[:, 0]
    cnt_ = np.bincount(_bi[_inb], minlength=nb)[:nb]
    hb = np.bincount(_bi[_inb], weights=hu[_inb], minlength=nb)[:nb] / np.maximum(cnt_, 1)
    hb[cnt_ == 0] = np.nan
    for w in mm:
        y_ = DC[w][ui]
        m = np.isfinite(y_) & np.isfinite(hb)
        if m.sum() < 8:
            REG[(w, ui)] = (np.nan, np.nan)
            continue
        xx = hb[m] - hb[m].mean()
        yy = y_[m] - y_[m].mean()
        REG[(w, ui)] = (float(np.sum(xx * yy) / np.sum(xx * xx)), float(np.corrcoef(xx, yy)[0, 1]))
say("")
say("=== does the tracked position follow the APPLIED displacement? (slope 1 = moves with it, 0 = removed) ===")
say("%5s %9s %9s %9s %7s %7s" % ("chan", "depth_um", "slope_bef", "slope_aft", "r_bef", "r_aft"))
_sb, _sa = [], []
for ui in range(len(units)):
    b_ = REG[("before", ui)]
    a_ = REG[("after", ui)]
    say("%5d %9.0f %9.2f %9.2f %7.2f %7.2f" % (units[ui]["chan"], yc[units[ui]["chan"]], b_[0], a_[0], b_[1], a_[1]))
    if np.isfinite(b_[0]) and np.isfinite(a_[0]):
        _sb.append(b_[0])
        _sa.append(a_[0])
if len(_sb) >= 3:
    _sb, _sa = np.array(_sb), np.array(_sa)
    say("median slope  before %.2f -> after %.2f   (n=%d)" % (np.median(_sb), np.median(_sa), len(_sb)))
    if len(_sb) >= 5:
        from scipy.stats import wilcoxon as _wx
        try:
            say("Wilcoxon |slope| before > after: p=%.4g" % _wx(np.abs(_sb), np.abs(_sa), alternative="greater").pvalue)
        except ValueError as exc:
            say("Wilcoxon not computed: %s" % exc)
np.savez(OUT / (PRE + seg + "_tracking.npz"), ctr=ctr, chan=np.array([u["chan"] for u in units]),
         depth=np.array([yc[u["chan"]] for u in units]), disp_before=DC["before"], disp_after=DC["after"],
         slope_before=np.array([REG[("before", i)][0] for i in range(len(units))]),
         slope_after=np.array([REG[("after", i)][0] for i in range(len(units))]),
         r_before=np.array([REG[("before", i)][1] for i in range(len(units))]),
         r_after=np.array([REG[("after", i)][1] for i in range(len(units))]),
         sd_before=sd["before"], sd_after=sd["after"])

(OUT / (PRE + seg + "_tracking.txt")).write_text(buf.getvalue(), encoding="utf-8")
print("\nfigures -> %s_examples.png, %s_tracking.png" % (seg, seg))
