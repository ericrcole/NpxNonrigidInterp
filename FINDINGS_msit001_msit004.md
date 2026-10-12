# Findings: validating this pipeline on two human cortical recordings (msit001, msit004)

This is the working record of the validation, kept in the order it happened. It includes
approaches that turned out to be wrong, each marked where it was corrected or withdrawn - that
history is the point, because several of the wrong turns are easy to repeat. For the conclusions
and the procedure to follow, read README.md first; use this document for evidence and detail.

Script paths below refer to interp_eval/scripts/ as the work was done; the same files are now in
scripts/ in this package, and per-recording constants moved to recordings.json.

---

# Interpolation-after-DREDge: fix and evaluation

Session of 2026-08-28 / 08-31. Recording: `Neuropix_2026-03-05_11-20-16`, ProbeA,
385-channel Open Ephys AP stream (384 neural + 1 sync), 36,114,264 frames @ 30 kHz
(1203.81 s). Motion: `nonrigid_displacement.csv`, 11 windows @ 250 Hz spanning
160-2160 um (in-brain only; above ~2160 um is saline).

## 1. Root-cause bug (fixed)

`shift_batch_on_disk2_dredge.m` used `ops.Nchan` (384) for byte offsets and `fread`
shapes on a **raw** file that is `ops.NchanTOT` (385) wide. Kilosort's original
operates on `ops.fproc`, which really is `Nchan` wide; repointing it at the raw file
without changing the count sheared every batch diagonally across the probe, mixed the
sync line into the neural channels, and slipped the write offset 2 bytes per batch.

Evidence: sync channel 384 was `{0,1}` (std 0.4) before and 1089 distinct values
(std 210) after. `manualDREDgeKS_NT128/NT300/NT320` are all invalid and should be
regenerated.

Fixes, all in `InterpolationAfterDREDge/Scripts/AP/`:
- `Kilosort25Configuration_Dredge.m` - `ops.NchanTOT` from `dredgeParams.NchanRaw` (385)
- `Dredge/shift_batch_on_disk2_dredge.m` - `NchanTOT` for file geometry, `Nchan` for
  kriging; sync channel written back untouched; short reads clamped (no EOF overrun);
  `dredgeParams.maxInterpDepth` leaves saline channels raw
- `Dredge/applyDredge.m` - `size(hread,1)` not `length()`; batch count clamped to file
  length; `refDisplacement` (explicit reference so a segment run matches a full run);
  `batchRange`, `skipCopy`, `skipDriftMap` options
- `Dredge/RealignRecordingUsingDredgeAndKilosort2point5.m` - bypasses
  `preprocessDataSub` when the drift map is skipped (it only feeds the drift map)

`pre_fix_uncommitted.patch` holds the pre-fix working-tree state.

## 2. Verification

Structural, on t=790-870 s: output size identical to input, sync channel bit-identical,
234 Hz batch-edge excess +0.3% (was large), neighbour correlation 0.205 -> 0.218
(was 0.24 -> 0.80).

Functional: a unit on channel 81 relocates to channel 85; applied shift at that depth
was +52.9 um, predicting +5.3 channels. Cost of kriging: peak amplitude 0.74x,
footprint power 0.86x, AP-band channel MAD 39.0 -> 29.9.

## 3. Choosing evaluation windows (`window_scan.py`)

Drift is concentrated in the first ~200 s (settling). Whole-recording median 60 s
slow excursion is only 11.8 um.

| window | slow p2p | pulse SD (0.7-1.6 Hz) |
|---|---|---|
| 0-60 s | 174 um | 23 um |
| **100-200 s** | **~80 um** | 5 um |
| **1080-1200 s** | ~30 um | **13-15.5 um** |
| 790-870 s | 22 um | 10.8 um |

t=790-870 s (the window in `compare_files_beforeafterinterp_xcorr3.ipynb`) is a quiet
stretch and a poor choice for the demonstration.

## 4. Result, high-drift window t=100-210 s

Drift map: the band at ~1400 um descends to ~1230 um over 110 s in the raw data and is
flat at ~1220 um after correction (`driftmap_merged_100-210s.png`).

| depth band | corr(before, DREDge) | corr(after, DREDge) | SD before | SD after |
|---|---|---|---|---|
| 600-800 um | -0.41 | -0.04 | 40.7 | 30.5 |
| 800-1000 um | +0.51 | +0.08 | 37.7 | **16.8** |
| 1200-1400 um | **+0.87** | **-0.11** | 45.8 | 29.3 |

Median depth SD 40.7 -> 29.3 um (1.39x). Split-half reliability of the depth trace at
1200-1400 um: 0.89 -> 0.17, i.e. the *reproducible* motion component is removed and
only tracker noise remains. That reliability statistic is the cleanest single measure,
since it needs no reference to DREDge.

## 5. Vertical streaks in the drift map = task-locked population bursts

Not a pipeline artifact. Burst modulation measured from spikes is 0.21 Hz; the MSIT
trial rate from the TTL log is 0.222 Hz (median ITI 4.51 s). 41-42% of all events fall
in the busiest 10% of 50 ms bins (Poisson: 10-15%), and event counts correlate 0.67
between before and after.

Secondary contributors: duplicate multi-channel detections (fixed by the 0.5 ms /
20-channel merge in `eval150b.py`: same-ms duplicates 30.8% -> 3.0% before, 61% ->
4.9% after), and noisier centroids for marginal events after kriging (spread-amplitude
correlation -0.20 before vs -0.40 after).

Consequence: spike availability is modulated at 0.22 Hz, which aliases into short-bin
depth estimates. The pulsatile analysis at t=1080-1200 s must weight bins by spike
count and be checked against high-rate bins only.

## 6. Open / next

1. Taper the shift to zero between ~1960-2160 um instead of the hard `maxInterpDepth`
   cut - the hard boundary makes a spurious band at ~2150 um (visible in
   `driftmap_highdrift_100-210s.png`).
2. Full-length run (~64 min at 13.6 ms/batch, plus a 27.8 GB copy).
3. Pulsatile analysis at t=1080-1200 s with bins short enough for 1 Hz.
4. Optimisation, after correctness is settled: hoist the constant `Kxx` out of the
   batch loop, cache `M` while `shifts` changes by < ~0.5 um, keep the file handle open
   instead of `fclose all`/`fopen` per batch.

### Dead ends - do not repeat
- Running `dredge_lfp` on the **AP** band bandpassed 0.5-300 Hz. The AP band is
  hardware high-passed at ~300 Hz, so this estimates motion from noise.
- `frame_slice`-ing a recording before `dredge_lfp`. It is an online estimator; a 200 s
  slice returned implausible 125-437 um SDs.
- Per-channel AP-band RMS profiles for tracking depth. Dominated by static per-channel
  gain; the argmax pins at lag 0 or at the search-band edge.

---

# Final pipeline (2026-09-04)

## Configuration

DREDge (`scripts/rerun_dredge.py`), on the LFP stream, in-brain channels only (<=2320 um):

    estimate_motion(rec, method='dredge_lfp', rigid=False,
                    win_step_um=800.0, win_scale_um=850.0)

NOT the defaults. `estimate_motion` defaults to `win_step_um=200, win_scale_um=300,
rigid=False`, which OVERRIDE `dredge_online_lfp`'s own defaults of 800/850/rigid=True.
The original estimate used the former and gave 11 windows over 2000 um.

MATLAB (`scripts/run_interp_final.m`):

    dredgeParams.NchanRaw        = 385;    % raw file width incl. sync channel
    dredgeParams.BatchSamplesNT  = 128;
    dredgeParams.maxInterpDepth  = 2160;   % brain/saline boundary
    dredgeParams.taperUm         = 200;    % ramp shift to zero over 1960-2160 um
    dredgeParams.motionLowpassHz = 0.05;   % drop estimator noise above the motion band
    dredgeParams.mTolUm          = 0.5;    % reuse kriging matrix below this change
    dredgeParams.refDisplacement = <median of the full-recording resampled trace>

Runs by streaming source -> dest with no pre-copy. 23-46 s per 130-150 s segment;
the full 1204 s recording extrapolates to roughly 6-8 min.

## Why these window settings (scripts/score_motion_estimates.py)

Estimates scored against depth traces measured from raw AP spikes, weighted by the
split-half reliability of that measurement. No interpolation run needed.

| config        | win | SD um | strain      | RMS>0.1Hz | r raw | r after <0.05Hz |
|---------------|-----|-------|-------------|-----------|-------|-----------------|
| current 200um |  11 |  89.1 | 0.26 (21%)  |     20.1  |  0.31 |            0.41 |
| rigid         |   1 | 127.3 | 0.00        |    114.3  |  0.19 |            0.30 |
| nonrigid_800  |   2 | 150.6 | 0.06 (0.2%) |    105.8  |  0.18 |            0.40 |
| nonrigid_400  |   5 | 192.8 | 0.15 (17%)  |    107.6  |  0.20 |            0.41 |

Once low-passed all nonrigid configs agree with the AP data equally well, so the extra
spatial freedom of 11 windows buys nothing - but only the 800 um config has physically
plausible strain. Rigid is genuinely worse, so some nonrigidity does help.

## Results

Slow-drift window t=100-210 s, the config change plus low-pass improves on the original:

| motion estimate                   | median depth SD | reliability |
|-----------------------------------|-----------------|-------------|
| original 200 um, no low-pass      | 40.7 -> 29.3 um (1.39x) | 0.22 -> 0.30 |
| nonrigid_800 + 0.05 Hz low-pass   | 39.2 -> 22.1 um (**1.78x**) | 0.22 -> 0.03 |

Quiet window t=1070-1200 s still degrades (median SD 23.9 -> 38.8 um). The drift map
shows the corrected bands acquiring a slow ~100 um upward curve over 130 s where the raw
bands are flat: the estimate contains slow drift there that the AP band does not
corroborate, so a temporal low-pass passes it straight through. More spatial
regularisation does not fix this - it is the accuracy of the estimate in quiet stretches,
not its smoothness.

Caveat on the quiet-window table: the correction shifts populations wholesale between
depth bands (1800-2000 um goes from 16,602 to 920 events, 2000-2200 um from 8,484 to
22,408), so fixed-depth-band before/after comparison there compares different neurons.
Unit-template tracking would be the right instrument.

## Verdict on this recording

Usable for the slow-drift demonstration (t=100-210 s) and nothing else. Whole-recording
median 60 s slow excursion is 11.8 um; drift is concentrated in the first ~200 s. The
cardiac/breathing-band content of the estimate is not credible - frequency wanders
0.9-1.7 Hz (86 +/- 16 bpm), amplitude grows toward the saline boundary, and there is no
ECG or pulse channel on the NI-DAQ to check against (AI0/AI2 are TTL-like, AI6/AI7 have
no periodic content). Use a different recording for the pulsatile demonstration.

---

# NPIX_MSIT_004 (OneBox) — the convincing demonstration (2026-09-11)

## Why this recording

| | msit001 | msit004 |
|---|---|---|
| median 60 s slow excursion | 11.8 um | **60.8 um** |
| nonrigid strain (800 um windows) | 0.06 | **0.008** |
| cardiac peak freq stability | 86 +/- 16 bpm (noise) | **92 +/- 0 bpm** |
| cardiac coherence across depth | +0.11 | **+1.00** |
| cardiac amplitude | ~20 um, edge-weighted | **38-53 um** |

msit004 has genuine cardiac pulsation; msit001 did not.

## Setup differences from msit001 (all verified, all matter)

- Valid range is **805-1745 s** file-relative, NOT 1870. From ~1745 s the LFP amplitude
  jumps 6-10x (withdrawal). Both streams start at absolute t=125.029 s.
- AP stream is **384 channels with no sync channel** (sync is on OneBox-ADC), so
  NchanRaw = 384. On msit001 it is 385 including sync.
- Probe geometry differs: x = {0, 32} aligned columns vs msit001's staggered
  {11, 27, 43, 59}. chanMap_msit004.mat was generated for it; the maps are NOT
  interchangeable because the kriging kernel uses x.
- motionLowpassHz = **5**, not 0.05. Here the 1-2 Hz band is real signal.
- No saline boundary is evident on this probe, so maxInterpDepth is left unset.

## Tracker corrections (each found by a result that could not be true)

1. **Column-swap.** Shifting the NCC template by one channel INDEX swaps the probe's
   two x-columns, which carry different amplitudes for the same unit. Negligible for
   msit001's ~3-index shifts, fatal for msit004's ~34-index ones. Now tracking in ROW
   space (uniform 20 um pitch), matching what average_across_direction does for DREDge.
2. **Degenerate peaks.** Some units gave a before-trace with <1 um p2p (a locked NCC)
   that exploded to 800+ um after. Gated on a BEFORE-side criterion only, so the
   exclusion cannot bias the before/after comparison.
3. **Ridge jumping.** The NCC ridge is ~200 um wide; an unconstrained argmax hops
   between local maxima and fabricates excursions. Neighbouring units were moving
   oppositely, which is impossible. Now constrained to +/-250 um of each unit's own
   median, applied identically to both files. This alone took ch118 from 0.43x to 3.29x.

## Result 1 — slow drift, full valid range, PASSES

Paired per-unit SD of the best-correlated channel position (same spike times in both
files, since interpolation never changes timing):

    n = 11 units    median 22.9 -> 7.9 um    2.90x reduction
    10 / 11 units improve    Wilcoxon signed-rank p = 0.027
    per-unit ratio: median 2.68x, IQR 1.75-3.74

(12 units tracked; ch150 excluded for a before-SD of 1.1 um - a locked tracker, not a
stable unit, in a recording whose estimate says 60+ um of drift. Including it:
2.48x, p = 0.117.)

## Result 2 — cardiac, FAILS, and the reason is measurable

Phase-binned tracking (12 bins, ~1400 cardiac cycles pooled) gives no net improvement:
median cardiac-locked amplitude 14.9 -> ~16 um. But the diagnostic is unambiguous:

    applied cardiac shift (DREDge) : amplitude 52.6 um, phase +3.14 rad
    measured in the AP band        : amplitude 14.9 um (IQR 6.8-18.6), phase ~0
    mean phase offset              : -3.07 rad = -176 deg

A phase offset of pi is the CORRECT relationship (the applied shift is ref - h, i.e.
minus the displacement), so the timing is right. The amplitude is not:

    implied cardiac gain = 14.9 / 52.6 = 0.28   ->   DREDge overestimates ~3.5x

The overcorrection signature is explicit in the units that worsen - their after-phase
flips to +/-pi, i.e. pushed through zero and out the other side:

    ch97   0.8 um @ +0.24  ->  23.6 um @ -3.12
    ch102 10.3 um @ -0.53  ->  25.1 um @ -3.13
    ch107  9.9 um @ -0.20  ->  26.0 um @ -3.10

while units with genuinely larger motion improve (ch46 29.1 -> 6.5, ch62 21.0 -> 6.6).

Likely cause: LFP-based motion estimation conflates true tissue displacement with
cardiac-locked LFP changes that are not motion (blood volume / pulsation artifact),
inflating the apparent displacement in that band. The slow band does not show this -
correction there works at ~2.9x.

## Next step

Add a band-split gain to applyDredge: keep the slow component at gain 1 and scale the
>0.5 Hz component by ~0.3, then re-run the cardiac analysis. Prediction: the units that
currently overshoot to phase +/-pi should collapse toward zero amplitude instead.

---

# Validation of the kriging step itself, and three pipeline defects (2026-09-12)

## Probe geometry: chanMap_default.mat is wrong for these recordings

Both recordings use probe **NP1016 (Neuropixels 1.0 Sapiens)**. Both `settings.xml` files store
ALIGNED columns, x = 27, 59, 27, 59 ... (y = 0, 0, 20, 20, ...), and probeinterface's official
NP1016 geometry agrees. `chanMap_default.mat` has the staggered x = {11, 27, 43, 59} layout of a
different NP1.0 variant and is wrong on **192 of 384 sites** (16 um in x). Every msit001
interpolation so far, and any Kilosort run that used that file, used the wrong geometry.
`chanMap_np1016.mat` (positions straight from settings.xml) is now used for msit001;
`chanMap_msit004.mat` (x = {0, 32}) is a constant offset of the true layout, which is equivalent
for kriging.

## Kriging validated against ground truth (scripts/synthetic_kriging.py)

Spike-like footprints exp(-distance/lambda) sampled at the true site positions, tissue moved by
0-120 um, corrected with the pipeline's own operator, compared with the unmoved truth. Averaged
over 5 lateral source positions, 3 depths and 2 footprint widths:

| | max centroid error | peak amplitude, on-grid / half-pitch | mean waveform error |
|---|---|---|---|
| correct NP1016 map | **0.06 um** | 0.96 / 0.84 | 7.5% |
| chanMap_default (wrong) | 0.14 um | 0.90 / 0.80 | 11.9% |
| no correction | 105 um | - | 88% |

Position is restored essentially exactly for any shift. The only cost is amplitude, worst (-16%)
when the shift lands half-way between rows. The wrong map does not bias depth - which is why
tracking results on msit001 are unchanged when it is re-run with the correct map (median SD ratio
1.06x vs 1.07x, same units) - but it adds waveform error and amplitude loss. The correct map
recovers fully every 20 um (one row); the wrong staggered map only every 40 um.

## DC-offset mixing artifact (fixed)

On msit004 the interpolated file's broadband std rose 1.7x and neighbour correlation fell. Split
by band, the change is **entirely below 150 Hz** (std 29.7 -> 108.2, correlation 0.83 -> 0.24);
the spike band behaves exactly as kriging should (std 0.88x, correlation 0.55 -> 0.78).

Cause: raw AP channels carry different DC offsets (SD 194 counts across channels). Kriging with
time-varying weights mixes those offsets, producing a slow, channel-specific signal that follows
the shift trajectory. Predicted from the measured offsets and the applied shift: 104.0; observed
extra <300 Hz power: 103.2; correlation over time r = 0.98.

Fix: `shift_batch_on_disk2_dredge.m` subtracts each channel's offset before kriging and restores
the output channel's own offset afterwards; `applyDredge.m` estimates offsets once from 20 probes
across the file (`dredgeParams.removeOffset`, default true). Validated on m4_card: <150 Hz std
back to 28.9 (raw 29.7), correlation 0.90; spike band (300-6000 Hz) new vs old differs by 0.9% RMS,
so all tracking results above are unaffected. It sat below Kilosort's high-pass (150-300 Hz band
was already unchanged, 0.97x), so it was never visible to sorting - but it was an artifact.

## Two latent bugs fixed (neither affected any run so far)

- `applyDredge.m` took the frame count from the OUTPUT file when one existed. With streaming output,
  a stale or partial output from an earlier run would silently truncate processing. Now always
  taken from the source.
- A rigid estimate (one window) passed a scalar shift into the depth taper, which would error.
  Now expanded to every channel first.

## Tracking caveat on msit001

With the row-space template tracker, msit001's high-drift window (t=100-210 s) gives a mixed
result: 7 units, 2 improve strongly (3.98x, 2.52x), 4 worsen, median SD 1.07x, median peak
correlation only 0.83 (msit004: 0.95-0.97). Earlier centroid-based methods on the same window gave
1.4-1.9x with split-half reliability collapsing 0.96 -> -0.06. The window has only 3744 events in
110 s, and its SD-based result depends on the method. msit004 carries the slow-drift
demonstration; msit001 is supporting evidence at best.

## What the msit004 motion estimate contains, band by band

Baseline estimate (nonrigid_800, LFP 0.5-250 Hz), RMS per band, median over the 4 depth windows:

| band | RMS | coherent across depth | verdict |
|---|---|---|---|
| < 0.05 Hz | 61.2 um | yes | slow drift - validated against spikes (2.9x SD reduction) |
| 0.05-0.1 Hz | 11.9 um | yes | slow drift |
| 0.1-0.5 Hz | 18.3 um | r = 1.00 | **confounded**: peak 0.20-0.25 Hz (13.4 +/- 1.0 /min) sits on the MSIT trial rate, 0.221 Hz (Stimulus_Shown, n=138, median ITI 4.53 s). Breathing and task-locked LFP change cannot be separated, and spikes cannot validate it. Do not apply. |
| 0.5-0.8 Hz | 3.0 um | - | negligible |
| 0.8-2 Hz | 36.3 um | r = 1.00 | cardiac fundamental, 1.53 Hz (92 bpm, stable) - real, but amplitude inflated relative to units |
| 2-5 Hz | 18.5 um | r = 0.99 | **cardiac 2nd harmonic**: peak 3.00 Hz (2 f0 = 3.06), phase-locked to twice the cardiac phase at 0.96 (0.31 to three times). Treat with the fundamental: cardiac band = 0.8-5 Hz, one gain. |
| 5-20 Hz | 8.8 um | - | not tissue motion |

The LFP high-pass 10 Hz variant (n800_hp10) is not a targeted fix: it shrinks EVERY band 2-3x
(slow 61 -> 37 um, cardiac 36 -> 15 um, 0.1-0.5 Hz 18 -> 6.5 um) and adds large noise above 5 Hz
(5-20 Hz 20.6 um, >20 Hz 53.4 um). It looks like a generally attenuated, noisier estimator, which
would under-correct slow drift. Scored against spikes below before any conclusion.

`scripts/derive_motion.py` builds band-selected, gain-scaled candidates from any estimate
(zero-phase over the full trace, so segments cut from it have no filter edge transients).

## The decisive per-unit test: does the tracked position follow the APPLIED displacement?

SD conflates two separate questions. Regressing each unit's tracked position on the displacement
actually applied at its depth (binned like the tracker) separates them:

- **after-slope minus before-slope** tests the interpolation alone: it should be -1 (the
  interpolation removes exactly the displacement it was given), whatever DREDge's accuracy.
- **before-slope** tests DREDge: +1 means the unit really moved as the estimate says.

msit001, t=100-210 s, correct NP1016 geometry (scripts/template_tracking.py, 5 s bins):

| unit | depth | slope before | slope after | r before | r after |
|---|---|---|---|---|---|
| ch55 | 1100 um | **1.10** | **-0.03** | 0.99 | -0.12 |
| ch65 | 1300 um | **1.30** | **0.23** | 0.97 | 0.48 |
| ch41 | 820 um | 0.31 | -0.59 | 0.45 | -0.71 |
| ch49 | 980 um | -0.16 | -0.87 | -0.14 | -0.58 |
| ch35 | 700 um | -0.86 | -1.86 | -0.87 | -0.97 |
| ch71 | 1420 um | -0.63 | -1.72 | -0.50 | -0.87 |

**Interpolation: after-before slope = -0.99 +/- 0.14 (range -1.14 to -0.71), n=6** - every unit shifted by
the applied displacement, at every depth. **DREDge: right for 2 units (1100-1300 um, which end at
~0), wrong for 3** (they moved against the estimate, and correction pushed them further). The
same depth pattern was seen independently with the centroid method earlier (600-800 um band
anti-correlated with DREDge, r = -0.41; 1200-1400 um, r = +0.87). msit001's mixed SD result is
therefore a DREDge accuracy problem at some depths in this window, not an interpolation problem.

## DREDge solver prior does not set the cardiac amplitude

n800_lt10 (thomas_kw lambda_t = 10, ten times the default temporal smoothness prior) is
indistinguishable from the baseline: cardiac-band RMS 36.2 vs 36.3 um, RMS above 0.1 Hz 45.9 vs
46.0 um, strain 0.008. The inflated cardiac amplitude is not a regularisation artifact of the
solver - it is in the pairwise LFP displacement measurements themselves.

## Scoring DREDge variants against spikes (scripts/cardiac_calib.py)

The cardiac tracker is calibrated from the baseline interpolation: content was moved by a known
cardiac shift q0, so each unit's measured phasor change must equal -g q0.

    tracker gain g: half 1 |g| = 0.41 at -13 deg (n=11), half 2 |g| = 0.30 at -6 deg (n=8)
    held-out model check: median |observed| vs |predicted| after-phasor 11.0 / 10.4 um and
    14.5 / 14.7 um - right on average, but per-unit relative error ~1 (noisy, n = 19 unit-halves)

Candidates, metrics on held-out halves (k = gain of the units' calibrated cardiac motion on the
candidate, RES = predicted residual cardiac power if applied, slow r = slow-trace agreement):

| candidate | cardiac k | phase | RES | RES (optimal complex rescale) | slow r |
|---|---|---|---|---|---|
| baseline nonrigid_800 (LFP 0.5 Hz) | 0.82 | +14 deg | 0.83 | 0.44 | 0.72 |
| LFP high-pass 10 Hz | 5.42 | **+72 deg** | 0.99 | 0.44 | **0.18** |
| LFP high-pass 20 Hz | 22.2 | +12 deg | 0.98 | 0.67 | **0.10** |
| lambda_t 10 / 100 | 0.82 / 0.86 | +14 / +15 deg | 0.83 / 0.81 | 0.44 / 0.43 | 0.72 |
| slow only | - | - | 1.00 | - | 0.72 |
| slow + cardiac x0.25 / 0.35 / 0.50 / 0.60 / 0.75 / 1.00 | | +14 deg | 0.70 / 0.63 / **0.57 / 0.57** / 0.61 / 0.83 | 0.44 | 0.72 |
| hybrid: baseline slow + hp10 / hp20 cardiac | | | 0.98 / 0.98 | | 0.72 |

Robust conclusions, independent of how the tracker gain is interpreted:
- **High-passing the LFP before DREDge is not a fix - it breaks the estimate.** Cardiac timing goes
  wrong (+72 deg at 10 Hz) and slow-drift agreement collapses (0.72 -> 0.18 / 0.10). Its smaller
  cardiac amplitude was lost signal, not removed contamination.
- **The DREDge solver prior does not matter** (lambda_t x10, x100 are identical to baseline).
- **The best achievable cardiac correction on msit004 is modest.** Even optimally rescaled, the
  residual is 0.44 of the power (~34% RMS reduction), limited by a consistent ~14 deg (~25 ms)
  LFP-vs-spike timing offset and by unit-to-unit differences in cardiac motion. A gain of
  0.5-0.6 on the cardiac band sits on the optimum plateau (RES 0.57); the unscaled estimate
  helps only a little (0.83).

Pending (injection test): whether |g| ~ 0.35 is the TRACKER under-reading - in which case units
move ~0.8x DREDge's cardiac amplitude, not 3.5x less as first measured - or the interpolation
under-realising cardiac-frequency shifts.

## Ground-truth injection test (scripts/injection_test.py, run_interp_syn.m)

120 s synthetic recording in the msit004 layout: 12 STATIC units (800-3000 um, biphasic spikes,
exp(-d/30 um) footprint, trough 400 counts, 20 Hz), Gaussian noise (sd 40), per-channel DC offsets
(sd 190). The real MATLAB pipeline is given a known motion identical at all depths,
h(t) = 40 sin(2 pi 0.02 t) + 50 sin(2 pi 1.53 t) um, so every unit in the output must sit at y0 - h(t).

Position measured at the known spike times with a high-SNR centroid - no template, no NCC.
The centroid estimator itself returns gain 1.000 on noiseless shifted footprints.

| | slow gain | cardiac gain | cardiac phase |
|---|---|---|---|
| input (before) | 0.001 | 0.007 | - |
| **pipeline output (after)** | **0.959** | **0.958** | **-2 deg** |
| exact | 1 | 1 | 0 deg |

**The interpolation realises cardiac-frequency shifts exactly as well as slow ones**: ~96% of the
applied displacement in both bands, with correct timing. The ~4% shortfall is identical across
bands, so it is frequency-independent (kriging shrinkage, not temporal smoothing).

Consequence: the gain |g| ~ 0.35 measured on msit004 cannot be the interpolation. It must be the
phase-binned NCC tracker under-reading amplitude - confirmed or refuted by running both trackers
on this same synthetic segment (below).

### Injection test, tracker stage

On the same synthetic segment (static units, known motion), with the pipeline's realised gains
(0.959 slow, 0.958 cardiac) divided out:

| tracker | measured before (truth 0) | combined gain | tracker-only gain |
|---|---|---|---|
| cardiac_calib.py (phase-binned NCC) | 0.0 um | 0.99 at -2 deg | **1.04** |
| template_tracking.py (time-binned NCC, slope) | slope 0.00 | slope -1.00 | **1.04** |

Both trackers read true displacement at gain ~1 on high-SNR synthetic units (trough ~10x noise).
This CONTRADICTS the interpretation that the real-data cardiac calibration |g| ~ 0.35 is generic
tracker attenuation. Since the interpolation code and trackers are both verified here, the
real-data discrepancy is either (a) tracking attenuating on real, lower-SNR units - under which
the calibrated estimate "units move 0.82x DREDge" still holds, because the attenuation cancels - or
(b) the real interpolation run applying less cardiac shift than intended (a path the synthetic test,
which wrote its own motion file, did not exercise). Resolved next with a tracker-free per-event test
on the real data.

### Real-data realisation, bias-free (chunk registration)

Each 20 ms chunk of the after file is cross-correlated over depth lags with the SAME chunk of the
before file. The output should be the input translated by the applied shift - including its noise,
which removes the noise-weight bias that attenuated a per-spike centroid version of this test
(-0.53 slow / -0.68 cardiac, larger displacement -> more attenuation, the signature of estimator
bias). Realised shift regressed on the applied slow (<0.1 Hz) and cardiac (0.8-5 Hz) components:

| | slow | cardiac | exact |
|---|---|---|---|
| synthetic positive control | -0.996 | -0.996 | -1 |
| **real msit004 baseline run** | **-1.017** | **-1.066** | -1 |

(The real fit's residual SD is inflated by a few parabolic-fit outliers on flat profiles; those pull
the coefficient toward 0, not toward -1. Clean re-run with clipped fits: see below.)

**The interpolation applied the full slow AND cardiac displacement on the real recording.** The
motion slice for the segment matches the source estimate exactly (max diff 0.000 um, best lag 0).
So the real-data cardiac calibration |g| ~ 0.35 is a property of the cardiac TRACKER on real data.

Identified cause: `cardiac_calib.py` builds each unit's template as the grand average over the whole
930 s, which on the raw "before" file is smeared by +/-100 um of slow drift. The synthetic units
had no drift, so the tracker read gain 1 there. On real data the raw "before" (smeared template)
and corrected "after" (sharp template) carry DIFFERENT tracker gains, so one calibration gain g does
not describe both, and the cardiac accuracy numbers above (units move 0.82x DREDge; best gain
0.5-0.6) inherit an unknown bias. Redesign: measure cardiac motion between two DRIFT-CORRECTED files
- a slow-only interpolation as "before", the final pipeline (slow + 0.6 x cardiac) as "after" - so
both templates are equally sharp and the files differ only in the cardiac correction.

Clean re-run (fits clipped, peaks on the search edge excluded, 12,608 of 13,500 chunk-blocks):

    real msit004 baseline run   realised = -1.007 * slow -0.961 * cardiac  (robust 3-MAD fit)
                                ordinary fit -0.996 / -0.967, residual SD 22 um

These match the synthetic pipeline gains (0.959 / 0.958) almost exactly. The interpolation applies
~100% of slow and ~96% of cardiac displacement on real data; the earlier -1.07 cardiac figure was
outlier-inflated.

    real msit004 FINAL run      realised = -0.999 * slow -0.988 * cardiac   residual SD 4.5 um
    (applied: slow 0-0.1 Hz band, SD 101.6 um; cardiac 0.8-5 Hz band x0.6, SD 24.7 um)

The final band-selected pipeline applied exactly the displacement it was specified, verified on the
real recording with no tracker involved.

## CARDIAC CONCLUSION, CORRECTED (drift-free calibration) - supersedes the cardiac numbers above

Measured between two drift-corrected files - slow-only interpolation vs the final pipeline, which
adds 0.6 x the baseline cardiac estimate - so unit templates are equally sharp in both
(`scripts/driftfree_cardiac.py`, `scripts/calibrate_gains.py`, 29 usable unit-halves):

    tracker gain |g| = 0.93 (half 1, -4 deg) and 0.96 (half 2, -5 deg)  ~ synthetic truth 1.04

**The template-smear diagnosis is confirmed**: once templates are equally sharp the tracker reads
correctly. The raw-vs-baseline calibration's g ~ 0.35 was the smear artifact, and everything
derived from it is withdrawn - in particular "units move 0.82x DREDge's cardiac estimate" and the
absolute RES values in the candidate table (its phase and slow-drift columns are gain-invariant
and remain valid).

**The original finding was right: DREDge's cardiac amplitude is ~3.7x too large on msit004.**

| | half 1 | half 2 |
|---|---|---|
| units move k x DREDge's cardiac estimate | **0.30** at +4 deg | **0.26** at +12 deg |
| median unit cardiac motion / DREDge | 14.1 um / 54.1 um (both halves) | |

Timing is right to within ~10 deg (~20 ms); the amplitude is not.

**Observed, no calibration involved** - final pipeline (cardiac gain 0.6) vs slow-only:
cardiac-locked power after/before = **1.20** (model prediction 1.16), 11/29 unit-halves reduced,
Wilcoxon p = 0.94. Gain 0.6 still overcorrects.

**Held-out gain calibration** (estimate on one half, evaluate on the other):

| gain applied to DREDge's cardiac band | residual cardiac power (half 2 / half 1) |
|---|---|
| 1.0 (DREDge as-is) | 3.52 / 2.36 - roughly triples cardiac motion |
| 0.6 (final run) | 1.30 / 1.01 |
| **calibrated 0.30 / 0.25** | **0.70 / 0.71** |

So: the interpolation applies cardiac shifts exactly (registration x0.96-0.99), but on this recording
DREDge's LFP-derived cardiac displacement is ~3.7x larger than the units' actual cardiac motion.
Applied at face value it adds cardiac motion; calibrated per recording against spikes it removes
~30% of cardiac power (~16% RMS). Cardiac correction is therefore not usable without a per-recording
check of DREDge's cardiac amplitude against the AP band.

## SLOW-DRIFT GAIN, like-for-like (raw vs slow-only interpolation) - `scripts/calibrate_gains.py`

Both files keep cardiac/fast motion; they differ only by DREDge's slow displacement (gain 1.0).
Tracked unit position regressed on that displacement, 10 units per half:

| | half 1 | half 2 |
|---|---|---|
| slope before -> after | 0.46 -> -0.46 | 0.46 -> -0.44 |
| r before -> after (gain-invariant) | 0.88 -> -0.82 | 0.70 -> -0.74 |
| tracker gain g = -(after - before) | 0.95 | 0.96 |
| **units move s x DREDge's slow drift** | **0.49** | **0.48** |

A gain calibrated on one half leaves ~1% of the slow drift on the other half; gain 1.0 leaves
104-107% (the same magnitude as before correction, with the sign flipped). The negative r after
correction proves the overshoot without any tracker calibration. The tracker gain matches
the drift-free cardiac pair (0.93-0.96) and synthetic truth.

Checked and excluded as a units bug: SpikeInterface's LFP geometry and the chanMap both have a
20 um row pitch spanning 0-3820 um, and the injection test verifies MATLAB + tracker scale end to end.
`scripts/lfp_scale_test.py` checks whether DREDge-LFP itself recovers a known injected shift at
full amplitude.

msit001 is not a single scale factor: slopes are ~1.1-1.3 before correction at 1100-1300 um
(DREDge right) but negative at 700/980/1420 um (DREDge wrong in that window), so a per-recording
scalar gain is only meaningful where depths agree, as on msit004.

Drift-free cardiac scores of the DREDge variants (`scripts/score_driftfree.py`):
nonrigid_800 k 0.28 at +8 deg, RES 2.94 (0.69 at best gain); lambda_t 10/100 identical; 200 um
windows k 0.14 at +38 deg, RES 4.90; LFP high-pass 10/20 Hz collapses the cardiac estimate
(1.5-10 um) with a +66 deg phase error at 10 Hz, and even at best gain removes <=10%. No DREDge
configuration tried measures cardiac amplitude correctly; the timing of nonrigid_800 is right.

## V2 PIPELINE (slow x0.49 + cardiac x0.28, gains calibrated held-out) - `scripts/v2_eval.py`

`derive_motion.py msit004_nonrigid_800_t805-1745 msit004_n800_v2_s49c28 "0,0.1,0.49;0.8,5,0.28"`,
segment m4_v2 (806-1744 s, hard-linked source), `run_interp_m4v2.m` (933 s).

**DREDge-LFP measures shifts at the right scale.** `scripts/lfp_scale_test.py` adds a known
40 um / 60 s + 10 um / 1.5 Hz shift to the real preprocessed LFP (900-1200 s) and reruns the same
nonrigid_800 estimate: recovered gain 1.00, r 1.00, in every window and both bands. The ~0.5-0.6x
slow and ~0.27x cardiac ratios are real disagreements between LFP-pattern motion and unit motion,
not an estimator or units bug.

**Slow drift, raw vs v2** (template_tracking m4_v2 810-1740): per-unit SD 23.5 -> 14.8 um (1.59x),
5/5 units improve, Wilcoxon p = 0.031; slope on the applied displacement 1.10 -> 0.30, r ~+0.85 -> ~+0.4.
The slope stays positive after correction, so gain 0.49 slightly UNDER-corrects.

**Bracketing the slow gain** (same 5 units, slope per unit of DREDge's slow displacement, after):
gain 1.0 -> -0.44 (r -0.8, over), gain 0.49 -> +0.15 (r +0.4, under). Assuming a similar tracker
gain in both corrected files, the zero crossing is at gain ~0.6. **Units move ~0.6x DREDge-LFP's slow
estimate on msit004**; 0.49 from the slow-only fit alone assumed tracker gain 0.95, while the bracket
implies ~1.1.

**Cardiac at gain 0.28 did NOT help**: residual cardiac power vs slow-only 1.35 (same tracker
gain) / 1.19 (if v2's tracker gain were 1.0); 15/29 unit-halves reduced; median |A|/|T| 0.91;
p = 0.71. The held-out prediction was 0.70. Final run (x0.6): 1.21. So on msit004 no gain applied
to DREDge-LFP's cardiac band gives a measurable reduction of the units' cardiac motion (~12 um).
The per-unit scatter around proportionality, plus the residual-slow interaction, removes the
modest gain the linear model predicted.

**Caution - SD across corrected runs:** on the 5 units tracked in every run, median SD is: raw 21.8,
slow-only 21.6, full-baseline (cardiac x1.0 and the 0.1-0.8 Hz band) 7.8, final 16.2, v2 14.8. The
7.8 for the baseline run is NOT a better correction: it is the run with the most added cardiac
motion (x2.9 power), whose smeared templates lower the tracker gain and under-read SD. Compare runs
only through gain-invariant quantities (sign of r, slope bracketing, same-gain pairs).

## CARDIAC MOTION IS HETEROGENEOUS ACROSS UNITS - `scripts/card_by_depth.py`, `card_by_depth.png`

The pooled "units move 0.27x DREDge" is an average over qualitatively different units
(drift-free calibration, T = B/g, per half):

| group | units (depth um) | own cardiac motion | vs DREDge (800 or 200 um windows) |
|---|---|---|---|
| move WITH DREDge | 920, 1540, 2400 (and 1320 at 34 um) | 47-52 um, stable across halves | ratio 0.93-1.10, phase +7..+20 deg |
| essentially still | 1720, 1920, 2760, 2860, 2960, 3060 | 0.3-1.1 um (<=~5 um allowing tracker quantisation) | ratio ~0.01, random phase |
| partial | 620, 720, 1140, 2020, 2300, 2500, 3180 | 10-34 um | 0.3-0.4 |

Neighbours 100 um apart fall in different groups (1920 still vs 2020 partial; 2300 partial vs
2400 full). The 200 um window estimate reproduces the same ratios, so finer DREDge windows do not
resolve it. This is not a tracker failure on the still units: e.g. the 1920 um unit (42k spikes)
tracks slow drift normally (slope 0.39, r 0.80, the same as the co-moving 920/2400 um units), and
the applied 17-30 um cardiac shifts are read correctly on it.

Consequence: no smooth displacement field - rigid, or nonrigid at 800 or 200 um - can remove
msit004's cardiac motion. Correction at gain 1 fixes the co-moving units and adds ~50 um to the still
ones; lower gains trade one against the other, so pooled cardiac power never falls (x0.28: 1.19-1.35,
x0.6: 1.21, x1.0: ~2.9). Slow drift, by contrast, is homogeneous (per-unit slopes 0.38-0.71 at all
depths, ~0.5-0.6x DREDge) and is corrected. A plausible physical reading - not proven here - is that
pulsatile motion is local, e.g. some neurons near the shaft move with the probe while others move
with the surrounding tissue, whereas slow drift is a bulk displacement.

## Spike-based DREDge cross-check - inconclusive (`scripts/dredge_ap_check.py`, `dredge_ap_vs_lfp.png`)

SpikeInterface preset "dredge" (locally_exclusive detection threshold 8, monopolar triangulation,
dredge_ap) on the AP stream 805-1745 s with the same 800 um windows. The estimate is not usable
as an independent amplitude reference on msit004:
- all four windows are nearly identical (effectively rigid);
- long straight interpolated stretches after ~1300 s, where peaks were too sparse to fit;
- -200 um transients and a +100 um offset in the first ~200 s.

Agreement with DREDge-LFP's slow component is r -0.15..0.24 at zero lag. Note the motion
object's temporal bins follow the OneBox stream clock (first bin 930.5 s for a slice starting at
805 s); `file_relative()` corrects this. The ~0.6x slow finding therefore rests on the unit-tracking
bracket (gain 1.0 over, 0.49 under) plus the v3 prediction test below, not on this check.

## V3 = RECOMMENDED msit004 OUTPUT: slow drift only, gain 0.6 - prediction test PASSED

`derive_motion.py msit004_nonrigid_800_t805-1745 msit004_n800_v3_s60 "0,0.1,0.6"`, segment m4_v3
(806-1744 s), `run_interp_m4v3.m` (452 s); output
`.../OneBox-122.ProbeA-AP/m4_v3/interpm4_v3/interpm4_v3.dat` (21,611,520,000 bytes = source size).

Gain 0.6 was chosen BEFORE this run from the bracket (1.0 over, 0.49 under). Raw vs v3
(template_tracking m4_v3 810-1740, same 5 units as v2):

| | gain 1.0 (slow-only) | gain 0.49 (v2) | **gain 0.6 (v3)** |
|---|---|---|---|
| median slope after, per applied shift | -0.45 | +0.30 | **-0.05** |
| median r after | ~-0.8 | ~+0.4 | **-0.09** |
| per-unit SD raw -> after (same 5 units) | 21.8 -> 21.6 | 23.5 -> 14.8 | **23.5 -> 13.8 um (1.70x, 5/5, p = 0.031)** |

Per-unit r after correction scatters around zero (+0.15, -0.41, +0.40, -0.09, -0.47): the remaining
slow residual is unit-specific rather than a shared under- or over-correction. v3 applies no cardiac
correction, because DREDge-LFP's cardiac band did not reduce unit cardiac motion at any gain.
