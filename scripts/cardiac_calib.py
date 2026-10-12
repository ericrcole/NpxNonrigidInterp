"""Tracker-calibrated scoring of DREDge motion estimates against AP-band spikes.

    python -u cardiac_calib.py measure <segdir> <t0_s> <t1_s>
    python -u cardiac_calib.py score   <segdir> <candidate_tag> [<candidate_tag> ...]

(select the recording with NPX_REC, as for the other scripts)

measure  (heavy, run once) Detect units on the BEFORE file. For each unit, and for each
         half of [t0, t1), measure the cardiac-locked depth modulation as a complex
         phasor in the before file and in the baseline-interpolated after file, and a
         10 s time-binned slow depth trace in both files. Cached to calib_<seg>.npz.
score    (fast) Project each candidate DREDge trace onto the same phase reference and
         report how well its cardiac and slow components match the AP-measured motion.

Why a calibration is needed
---------------------------
Spike-measured cardiac amplitude (~15 um) looked 3.5x smaller than DREDge's (~53 um),
but the overshooting units ended at ~25 um rather than the ~37 um pure overcorrection
predicts - so the NCC tracker probably under-reads amplitude. The baseline after-file
fixes that internally: interpolation moved content by a KNOWN displacement q0 (the
applied estimate, verified independently by profile cross-correlation), so for every
unit the measured phasor change must satisfy

        A - B = -g * q0

where g is the tracker gain (complex least squares across units; angle should be ~0).
The true AP-side cardiac displacement of a unit is then T = B / g, and a candidate
estimate q_c is judged by
        k      gain of T on q_c (ideal 1; the baseline is expected well below 1)
        RES    |T - q_c|^2 / |T|^2  residual motion if q_c were applied (<1 helps)
        RESk   same after rescaling q_c by its best gain (shape/timing match only)
g is fitted on one half of the recording and the metrics reported on the other half,
in both directions, so the calibration is never evaluated on the data it was fitted to.

Phasor convention: a curve y(phi) is fitted as Re(z * exp(i*phi)). The phase reference
phi(t) is the Hilbert phase of the depth-averaged 0.8-2 Hz component of the BASELINE
estimate, whose timing was validated (before-phase vs applied shift offset = pi).
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import butter, filtfilt, hilbert

import recinfo

_R = recinfo.get()
AP = _R["ap"]
NCH = _R["nch"]
NCHN = _R["nchan_neural"]
FS = 30000.0
MFS = 250.0
BASELINE_TAG = _R.get("baseline_motion_tag")   # set in recordings.json after the DREDge run

THR = 6.0
MERGE_MS = 0.5
MERGE_CH = 10          # rows
TMPL_CH = 6            # rows
SEARCH_CH = 25         # rows
WF_PRE, WF_POST = 20, 40
MAX_SAMP_LAG = 3
NPHASE = 12
MIN_SPK_PHASE = 8      # per phase bin, per half
MIN_PHASE_BINS = 10    # of 12, tracked in BOTH files, for a unit-half to be used
MIN_PEAK_R = 0.7
RESTRICT_UM = 250.0
EDGE_ROWS = 3
MIN_LOCK_UM = 2.0      # before-side only: p2p below this = locked NCC, not a unit
N_UNITS = 24
MIN_RATE = 0.3
MIN_SNR = 4.0
MIN_SEP_CH = 5
SLOW_BIN = 10.0
CARD_BAND = (0.8, 2.0)

yc_ch, _ = recinfo.chanmap_yx()
yc, RMAT = recinfo.row_matrix()
NROW = len(yc)
ROW_OF = np.searchsorted(yc, yc_ch[:NCHN])
UM_PER_ROW = float(np.median(np.diff(yc)))
bp = butter(3, [300 / (FS / 2), 6000 / (FS / 2)], btype="band")
bcard = butter(3, [CARD_BAND[0] / (MFS / 2), CARD_BAND[1] / (MFS / 2)], btype="band")


def phasor(y, ph):
    m = np.isfinite(y) & np.isfinite(ph)
    if m.sum() < 6:
        return np.nan + 0j
    X = np.c_[np.cos(ph[m]), np.sin(ph[m]), np.ones(m.sum())]
    c, *_ = np.linalg.lstsq(X, y[m], rcond=None)
    return complex(c[0], -c[1])


def ncc(T, W):
    nt, ns = T.shape
    core = slice(MAX_SAMP_LAG, ns - MAX_SAMP_LAG)
    Tz = T[:, core] - T[:, core].mean()
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


def parab(y):
    i = int(np.nanargmax(y))
    if 0 < i < len(y) - 1 and np.all(np.isfinite(y[i - 1:i + 2])):
        a, b, c = y[i - 1], y[i], y[i + 1]
        den = a - 2 * b + c
        if abs(den) > 1e-12:
            return i + 0.5 * (a - c) / den, float(b)
    return float(i), float(y[i])


def track(template, means, counts, min_count):
    """Restricted-argmax NCC depth track. means: (nbins, nband, nsamp)."""
    T = template[SEARCH_CH:SEARCH_CH + 2 * TMPL_CH + 1]
    profs = {}
    for b in range(means.shape[0]):
        if counts[b] >= min_count:
            profs[b] = ncc(T, means[b])
    out = np.full(means.shape[0], np.nan)
    if not profs:
        return out
    nlag = len(next(iter(profs.values())))
    c0 = int(np.median([int(np.nanargmax(p)) for p in profs.values()]))
    half = int(round(RESTRICT_UM / UM_PER_ROW))
    lo, hi = max(0, c0 - half), min(nlag, c0 + half + 1)
    for b, p in profs.items():
        i, v = parab(p[lo:hi])
        i += lo
        if v < MIN_PEAK_R or i < EDGE_ROWS or i > nlag - 1 - EDGE_ROWS:
            continue
        out[b] = (i - SEARCH_CH) * UM_PER_ROW
    return out


# ============================== measure ==============================
def measure(seg, T0, T1):
    S = AP / seg
    SEG_T0 = float((S / "segment_t0.txt").read_text().strip())
    # NPX_CALIB_BEFORE / NPX_CALIB_AFTER override the file pair. Used to compare two DRIFT-
    # CORRECTED files (slow-only vs slow+cardiac): with raw "before" data the 930 s grand-
    # average template is smeared by +/-100 um of slow drift while the corrected "after"
    # template is sharp, so the two files carry different tracker gains and a single g is wrong.
    B = Path(os.environ["NPX_CALIB_BEFORE"]) if os.environ.get("NPX_CALIB_BEFORE") else S / "continuous.dat"
    A = (Path(os.environ["NPX_CALIB_AFTER"]) if os.environ.get("NPX_CALIB_AFTER")
         else S / ("interp" + seg) / ("interp" + seg + ".dat"))
    print("before file: %s" % B, flush=True)
    print("after file : %s" % A, flush=True)
    nfr = B.stat().st_size // (NCH * 2)
    mm = {"before": np.memmap(B, dtype="int16", mode="r", shape=(nfr, NCH)),
          "after": np.memmap(A, dtype="int16", mode="r", shape=(nfr, NCH))}

    D = np.loadtxt(S / "nonrigid_displacement.csv", delimiter=",")
    D = D[None, :] if D.ndim == 1 else D
    sig = filtfilt(*bcard, D, axis=1).mean(0)
    phase_u = np.unwrap(np.angle(hilbert(sig)))
    mt = SEG_T0 + np.arange(len(sig)) / MFS
    f0 = float(np.median(np.diff(phase_u)) * MFS / (2 * np.pi))
    print("reference cardiac rhythm %.2f Hz (%.0f bpm)" % (f0, 60 * f0), flush=True)

    def load(w, s0, n):
        x = np.asarray(mm[w][s0:s0 + n, :NCHN]).astype(np.float32)
        x = filtfilt(*bp, x, axis=0)
        return x - np.median(x, axis=1, keepdims=True)

    F0, F1 = int((T0 - SEG_T0) * FS), int((T1 - SEG_T0) * FS)
    BLK = int(10 * FS)

    # ---- detect on the before file ----
    sdc = sdr = None
    TT, CC, AA = [], [], []
    for s in range(F0, F1, BLK):
        e = min(F1, s + BLK)
        x = load("before", s, e - s)
        if sdc is None:
            sdc = np.median(np.abs(x), axis=0) / 0.6745
            sdc[sdc < 1e-6] = np.inf
            sdr = sdc @ RMAT
        mid = x[1:-1]
        ii, ch = np.nonzero((mid < -THR * sdc) & (mid < x[:-2]) & (mid < x[2:]))
        ii = ii + 1 + s
        amp = x[ii - s, ch]
        ch = ROW_OF[ch]
        o = np.lexsort((ch, ii))
        ii, ch, amp = ii[o], ch[o], amp[o]
        wnd = int(MERGE_MS * 1e-3 * FS)
        keep = np.ones(len(ii), bool)
        for j in range(len(ii)):
            if keep[j]:
                g = np.nonzero((ii > ii[j] - wnd) & (ii < ii[j] + wnd) & (np.abs(ch - ch[j]) <= MERGE_CH))[0]
                if len(g) > 1:
                    bst = g[np.argmin(amp[g])]
                    keep[g] = False
                    keep[bst] = True
        TT.append(ii[keep])
        CC.append(ch[keep])
        AA.append(amp[keep])
        del x
    t_all, c_all, a_all = np.concatenate(TT), np.concatenate(CC), np.concatenate(AA)
    dur = T1 - T0
    print("detected %d merged events" % len(t_all), flush=True)

    lim = TMPL_CH + SEARCH_CH
    cands = []
    for c in range(lim, NROW - lim):
        k = np.abs(c_all - c) <= 3
        n = int(k.sum())
        if n / dur < MIN_RATE:
            continue
        snr = float(np.median(-a_all[k]) / sdr[c])
        if snr >= MIN_SNR:
            cands.append(dict(row=c, n=n, t=t_all[k]))
    cands.sort(key=lambda d: -d["n"])
    units = []
    for d in cands:
        if all(abs(d["row"] - u["row"]) >= MIN_SEP_CH for u in units):
            units.append(d)
        if len(units) >= N_UNITS:
            break
    units.sort(key=lambda d: d["row"])
    nu = len(units)
    print("%d units" % nu, flush=True)

    # ---- accumulate phase-binned (per half) and slow time-binned waveforms ----
    Tmid = 0.5 * (T0 + T1)
    nband, nsamp = 2 * lim + 1, WF_PRE + WF_POST + 1
    nslow = int(np.ceil(dur / SLOW_BIN))
    pedges = np.linspace(-np.pi, np.pi, NPHASE + 1)
    ACC = {w: np.zeros((nu, 2, NPHASE, nband, nsamp), np.float32) for w in mm}
    CNT = np.zeros((nu, 2, NPHASE), np.int64)
    ACS = {w: np.zeros((nu, nslow, nband, nsamp), np.float32) for w in mm}
    CNS = np.zeros((nu, nslow), np.int64)
    for w in mm:
        print("accumulating %s ..." % w, flush=True)
        for s in range(F0, F1, BLK):
            e = min(F1, s + BLK)
            a0, a1 = max(0, s - WF_PRE - 5), min(nfr, e + WF_POST + 5)
            x = load(w, a0, a1 - a0) @ RMAT
            for ui, u in enumerate(units):
                band = np.arange(u["row"] - lim, u["row"] + lim + 1)
                st = u["t"][(u["t"] >= s) & (u["t"] < e)]
                st = st[(st - WF_PRE >= a0) & (st + WF_POST < a1)]
                if not len(st):
                    continue
                idx = (st - a0)[:, None] + np.arange(-WF_PRE, WF_POST + 1)[None, :]
                sn = x[:, band][idx].transpose(0, 2, 1)
                tt = st / FS + SEG_T0
                ph = np.angle(np.exp(1j * np.interp(tt, mt, phase_u)))
                pb = np.clip(np.digitize(ph, pedges) - 1, 0, NPHASE - 1)
                hb = (tt >= Tmid).astype(int)
                sb = np.clip(((tt - T0) / SLOW_BIN).astype(int), 0, nslow - 1)
                for h in (0, 1):
                    for b in range(NPHASE):
                        m = (hb == h) & (pb == b)
                        if m.any():
                            ACC[w][ui, h, b] += sn[m].sum(0)
                            if w == "before":
                                CNT[ui, h, b] += int(m.sum())
                for b in np.unique(sb):
                    m = sb == b
                    ACS[w][ui, b] += sn[m].sum(0)
                    if w == "before":
                        CNS[ui, b] += int(m.sum())
            del x

    # ---- track ----
    pc = pedges[:-1] + np.diff(pedges) / 2
    zB = np.full((nu, 2), np.nan + 0j)
    zA = np.full((nu, 2), np.nan + 0j)
    ntrk = np.zeros((nu, 2), int)
    lock = np.zeros((nu, 2))
    slowB = np.full((nu, nslow), np.nan)
    slowA = np.full((nu, nslow), np.nan)
    for ui in range(nu):
        for w, zz in (("before", zB), ("after", zA)):
            wt = CNT[ui].sum(0).astype(float)
            tot = ACC[w][ui].sum(0)                                  # (NPHASE, nband, nsamp)
            tmpl = tot.sum(0) / max(wt.sum(), 1)
            for h in (0, 1):
                means = ACC[w][ui, h] / np.maximum(CNT[ui, h], 1)[:, None, None]
                y = track(tmpl, means, CNT[ui, h], MIN_SPK_PHASE)
                if w == "before":
                    lock[ui, h] = np.nanmax(y) - np.nanmin(y) if np.isfinite(y).any() else 0.0
                    ntrk[ui, h] = int(np.isfinite(y).sum())
                else:
                    ntrk[ui, h] = min(ntrk[ui, h], int(np.isfinite(y).sum()))
                zz[ui, h] = phasor(y - np.nanmean(y), pc) if np.isfinite(y).sum() >= 6 else np.nan
            means_s = ACS[w][ui] / np.maximum(CNS[ui], 1)[:, None, None]
            tmpl_s = ACS[w][ui].sum(0) / max(CNS[ui].sum(), 1)
            ys = track(tmpl_s, means_s, CNS[ui], 4)
            (slowB if w == "before" else slowA)[ui] = ys - np.nanmean(ys)
        print("  unit row %3d (%4.0f um, %5d spk): |B| %5.1f %5.1f   |A| %5.1f %5.1f   bins %2d/%2d"
              % (units[ui]["row"], yc[units[ui]["row"]], units[ui]["n"],
                 abs(zB[ui, 0]), abs(zB[ui, 1]), abs(zA[ui, 0]), abs(zA[ui, 1]),
                 ntrk[ui, 0], ntrk[ui, 1]), flush=True)
    ok = (ntrk >= MIN_PHASE_BINS) & (lock >= MIN_LOCK_UM)
    out = recinfo.RESULTS / ("calib_%s_%s.npz" % (recinfo.key(), os.environ.get("NPX_CALIB_TAG", seg)))
    np.savez(out, rows=np.array([u["row"] for u in units]),
             depth=np.array([yc[u["row"]] for u in units]),
             nspk=np.array([u["n"] for u in units]), zB=zB, zA=zA, ok=ok, ntrk=ntrk,
             lock=lock, slowB=slowB, slowA=slowA,
             slow_ctr=T0 + (np.arange(nslow) + 0.5) * SLOW_BIN,
             T=np.array([T0, Tmid, T1]), mt=mt[::5], phase_u=phase_u[::5],
             seg_t0=SEG_T0, baseline=BASELINE_TAG)
    print("usable unit-halves: %d / %d  -> %s" % (ok.sum(), ok.size, out.name), flush=True)


# ============================== score ==============================
def candidate_at(tag, depth):
    M = _R["motion"]
    D = np.loadtxt(M / (tag + "_displacement.csv"), delimiter=",")
    D = D[None, :] if D.ndim == 1 else D
    c = np.atleast_1d(np.loadtxt(M / (tag + "_window_centers.csv"), delimiter=","))
    t0 = float((M / (tag + "_t0.txt")).read_text().strip())
    if D.shape[0] == 1:
        X = np.repeat(D, len(depth), axis=0)
    else:
        X = PchipInterpolator(c, D, axis=0)(np.clip(depth, c[0], c[-1]))
    return X, t0 + np.arange(D.shape[1]) / MFS


def cardiac_phasors(tag, depth, halves, mt, phase_u):
    """Cardiac phasor of a candidate estimate at each unit depth, per half (shared with plots)."""
    X, tx = candidate_at(tag, depth)
    Xc = filtfilt(*bcard, X, axis=1)
    q = np.full((len(depth), 2), np.nan + 0j)
    for h, (a, b) in enumerate(halves):
        idx = np.nonzero((tx >= a + 2) & (tx < b - 2))[0][::5]
        ph = np.angle(np.exp(1j * np.interp(tx[idx], mt, phase_u)))
        for u in range(len(depth)):
            q[u, h] = phasor(Xc[u, idx], ph)
    return q


def score(seg, tags):
    cal = np.load(recinfo.RESULTS / ("calib_%s_%s.npz" % (recinfo.key(), seg)),
                  allow_pickle=True)
    depth, zB, zA, ok = cal["depth"], cal["zB"], cal["zA"], cal["ok"]
    T0, Tmid, T1 = cal["T"]
    mt, phase_u = cal["mt"], cal["phase_u"]
    halves = [(T0, Tmid), (Tmid, T1)]
    base = str(cal["baseline"])

    def card_phasors(tag):
        X, tx = candidate_at(tag, depth)
        Xc = filtfilt(*bcard, X, axis=1)
        q = np.full((len(depth), 2), np.nan + 0j)
        for h, (a, b) in enumerate(halves):
            m = (tx >= a + 2) & (tx < b - 2)
            idx = np.nonzero(m)[0][::5]
            ph = np.angle(np.exp(1j * np.interp(tx[idx], mt, phase_u)))
            for u in range(len(depth)):
                q[u, h] = phasor(Xc[u, idx], ph)
        return q, X, tx

    q0, X0, tx0 = card_phasors(base)

    # tracker gain from the baseline after-file, per half
    g = np.full(2, np.nan + 0j)
    for h in (0, 1):
        u = ok[:, h] & np.isfinite(q0[:, h]) & np.isfinite(zA[:, h]) & np.isfinite(zB[:, h])
        d = -(zA[u, h] - zB[u, h])
        g[h] = np.sum(d * np.conj(q0[u, h])) / np.sum(np.abs(q0[u, h]) ** 2)
    print("\ntracker gain g from the baseline after-file (A - B = -g q0):")
    for h in (0, 1):
        print("  half %d: |g| = %.2f, angle %+.0f deg  (n=%d unit-halves)"
              % (h + 1, abs(g[h]), np.degrees(np.angle(g[h])), ok[:, h].sum()))

    # model check: does the calibrated model predict what the interpolation actually left?
    print("\nmodel check on held-out half (calibrate g on the other half):")
    for h in (0, 1):
        gc = g[1 - h]
        u = ok[:, h] & np.isfinite(q0[:, h])
        pred = gc * (zB[u, h] / gc - q0[u, h])
        obs = zA[u, h]
        rel = np.sqrt(np.sum(np.abs(obs - pred) ** 2) / np.sum(np.abs(obs) ** 2))
        print("  half %d: predicted vs observed after-phasor, relative error %.2f; "
              "median |obs| %.1f um, median |pred| %.1f um"
              % (h + 1, rel, np.median(np.abs(obs)), np.median(np.abs(pred))))

    print("\n%-34s %7s %7s | %6s %8s %6s %6s | %6s %6s"
          % ("candidate", "|q|um", "|T|um", "|k|", "k_ang", "RES", "RESk", "slow_r", "slow_k"))
    RESULTS = []
    for tag in [base] + [t for t in tags if t != base]:
        try:
            q, X, tx = card_phasors(tag)
        except OSError as exc:
            print("%-34s missing (%s)" % (tag, exc))
            continue
        rows = []
        for h in (0, 1):
            gc = g[1 - h]                                   # held-out calibration
            u = ok[:, h] & np.isfinite(q[:, h])
            T = zB[u, h] / gc
            qc = q[u, h]
            k = np.sum(T * np.conj(qc)) / np.sum(np.abs(qc) ** 2)
            res = np.sum(np.abs(T - qc) ** 2) / np.sum(np.abs(T) ** 2)
            resk = np.sum(np.abs(T - k * qc) ** 2) / np.sum(np.abs(T) ** 2)
            rows.append((np.median(np.abs(qc)), np.median(np.abs(T)), abs(k),
                         np.degrees(np.angle(k)), res, resk))
        r = np.nanmean(np.array(rows), axis=0)
        # slow: candidate <0.05 Hz, bin-averaged, vs before slow trace (tracker-uncalibrated)
        Xs = filtfilt(*butter(3, 0.05 / (MFS / 2), btype="low"), X, axis=1)
        ctr = cal["slow_ctr"]
        rr, kk = [], []
        for u in range(len(depth)):
            y = cal["slowB"][u]
            v = np.array([Xs[u, (tx >= c - 5) & (tx < c + 5)].mean() for c in ctr])
            m = np.isfinite(y) & np.isfinite(v)
            if m.sum() >= 10 and np.std(v[m]) > 0:
                vv = v[m] - v[m].mean()
                yy = y[m] - y[m].mean()
                rr.append(np.corrcoef(yy, vv)[0, 1])
                kk.append(np.sum(yy * vv) / np.sum(vv * vv))
        print("%-34s %7.1f %7.1f | %6.2f %+7.0f %6.2f %6.2f | %6.2f %6.2f"
              % (tag.replace(recinfo.key() + "_", ""), r[0], r[1], r[2], r[3], r[4], r[5],
                 np.median(rr) if rr else np.nan, np.median(kk) if kk else np.nan))
        RESULTS.append(dict(tag=tag, q_um=float(r[0]), T_um=float(r[1]), k=float(r[2]),
                            k_deg=float(r[3]), RES=float(r[4]), RESk=float(r[5]),
                            slow_r=float(np.median(rr)) if rr else None,
                            slow_k=float(np.median(kk)) if kk else None))
        (recinfo.RESULTS / ("score_%s_%s.json" % (recinfo.key(), seg))).write_text(
            json.dumps(dict(g=[[float(abs(x)), float(np.degrees(np.angle(x)))] for x in g],
                            candidates=RESULTS), indent=1))
    print("\n|k|: gain of true AP cardiac motion on the candidate (ideal 1).  RES: residual")
    print("fraction if applied as-is (<1 helps, 0 perfect).  RESk: residual after optimal")
    print("rescaling (timing/shape match only).  slow_r/slow_k: agreement with the slow")
    print("before-trace (slow_k is NOT tracker-calibrated).  Metrics use held-out halves.")


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "measure":
        measure(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
    elif mode == "score":
        score(sys.argv[2], sys.argv[3:])
    else:
        raise SystemExit(__doc__)
