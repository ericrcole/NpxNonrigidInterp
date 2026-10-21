# Neuropixels motion correction: DREDge + kriging interpolation, verified against spikes

Estimate probe/tissue motion from the LFP with DREDge, apply it to the AP band by kriging
interpolation in MATLAB, and - the part that matters - **check against the spikes that the
correction actually removed motion before you sort the recording**.

The checking is not optional. On the two human cortical recordings this was built for, the
interpolation applied the requested shift essentially exactly (x0.96-1.01, verified by registration
and by injecting known motion), yet applying DREDge's estimate at face value made units move
*backwards* by almost as much as they had moved forwards: units there moved only ~0.6x the slow
displacement DREDge reported. One number per recording (`slow_gain`) fixed that. See
[FINDINGS_msit001_msit004.md](FINDINGS_msit001_msit004.md) for the full evidence.

Built on [InterpolationAfterDREDge](https://github.com/williamunoz/InterpolationAfterDREDge) by
Richard Hardstone and William Munoz (CC0 1.0), which adapts Kilosort 2.5, with motion estimated by
DREDge through SpikeInterface. This is a standalone internal copy, not a fork: see
[CREDITS.md](CREDITS.md) for the full provenance and what is new here, and check upstream for fixes.

## What you get

| | |
|---|---|
| Input | an Open Ephys recording with an AP stream and an LFP stream |
| Output | `interp<segment>.dat` - the AP data with motion removed, same size and format as the input, ready for Kilosort |
| Verification | per-unit figures and numbers showing whether motion actually went down |

## Requirements

- **MATLAB** (developed on R2024b) and a clone of
  [InterpolationAfterDREDge](https://github.com/williamunoz/InterpolationAfterDREDge), which
  supplies the supporting functions. Four of its files are modified here and shadow the
  originals automatically (see [matlab/README.md](matlab/README.md)).
- **Python 3.11+** with `spikeinterface`, `torch` (DREDge needs it), `numpy`, `scipy`, `matplotlib`,
  `h5py`. A conda env is easiest:
  `conda create -n dredge python=3.11 && pip install spikeinterface[full] torch h5py`
- **Disk**: each interpolated output is the same size as the segment you interpolate (~21 GB per
  15 min of 384-channel data). You need room for the segment plus its output.
- **Time** on a 32-core desktop: DREDge ~10 min per recording, interpolation ~15-30 min per
  15 min segment (the ADC phase shift roughly doubles it), unit tracking ~15 min per segment.

## Setup

1. Clone/copy this folder and the InterpolationAfterDREDge repo.
2. Nothing to copy: `run_interp.m` puts this package's `matlab/` folder ahead of the repo on the
   MATLAB path, so its modified files shadow the repo's originals (see
   [matlab/README.md](matlab/README.md)).
3. Set `INTERP_REPO` to the repo path, or pass it as the 2nd argument to `run_interp`.
4. Add your recording to [recordings.json](recordings.json) (copy an existing entry):

| field | meaning |
|---|---|
| `ap` | folder containing the AP `continuous.dat` |
| `lfp_base`, `lfp_stream` | what SpikeInterface opens for the LFP, e.g. `Record Node 123#OneBox-122.ProbeA-LFP` |
| `motion` | where motion estimates are written |
| `chanmap` | `.mat` with `xcoords`/`ycoords`; relative paths resolve inside this package |
| `nch` | channels **in the file**, sync channel included (385 on PXI, 384 on OneBox) |
| `nchan_neural` | channels to interpolate (384) |
| `valid` | `[t0, t1]` in file-relative seconds - the part of the recording where the probe was in the brain and stationary enough to be worth correcting |
| `max_depth_um` | ignore sites above this (e.g. in saline); `null` to use all |
| `baseline_motion_tag` | filled in after step 1 below |
| `slow_gain`, `cardiac_gain` | filled in after calibration; this is the output of the whole procedure |

5. Check it before running anything heavy:

```bash
NPX_REC=mynewrec python scripts/check_setup.py
```

It verifies the file size matches `nch`, prints the probe geometry, and opens the LFP stream.
**Confirm the printed column x positions against your probe's `settings.xml`** - a wrong channel map
silently corrupts everything downstream. Every script picks the recording from `NPX_REC`, and
everything they write goes to `results/`.

## Running it

All commands run from the package root. Replace `mynewrec`, times and names.

### 1. Estimate motion from the LFP

```bash
NPX_REC=mynewrec python -u scripts/dredge_recording.py mynewrec 805 1745 nonrigid_800
```

DREDge on the LFP stream (bandpass 0.5-250 Hz -> phase shift -> resample to 250 Hz -> average the
probe's two columns), restricted to your valid range. Writes
`<motion>/mynewrec_nonrigid_800_t805-1745_{displacement,window_centers,t0,config}.csv/json`.
Put that tag in `recordings.json` as `baseline_motion_tag`.

Estimate motion from the **LFP**, never the AP band: the AP band is hardware high-passed at 300 Hz,
and DREDge on it returns noise. Use the `nonrigid_800` config: SpikeInterface's `estimate_motion`
silently overrides DREDge's 800 µm windows with 200 µm, which gave implausible tissue strain.

### 2. Select the bands and the gain

```bash
NPX_REC=mynewrec python scripts/derive_motion.py mynewrec_nonrigid_800_t805-1745 mynewrec_g100 "0,0.1,1.0"
```

`"f_lo,f_hi,gain;..."` selects frequency bands and scales each one; `f_lo = 0` means low-pass.
`"0,0.1,1.0"` is slow drift only, at full amplitude - the right starting point. Filtering happens
on the whole estimate, so segments cut out of it later carry no filter edge transients.

Apply only what you can verify. Slow drift is checkable and correctable. Cardiac motion was not
correctable on our recordings (see "Cardiac" below). Anything else - e.g. a 0.1-0.5 Hz band sitting
on your task's trial rate - is task-locked LFP change you cannot distinguish from motion; leave it out.

### 3. Cut a segment to work on

```bash
NPX_REC=mynewrec python scripts/make_segment.py 806 1744 seg_g100 mynewrec_g100
```

Writes `<ap>/seg_g100/` with `continuous.dat`, the matching motion slice, `chanMap.mat`,
`ref_displacement.csv` and `interp_params.json`. The displacement reference comes from the whole
recording, so a segment run equals the corresponding batches of a full-length run.

For calibration, a few minutes of data is enough and much faster; use the full valid range for the
final output. To make a second segment over the **same** time range without copying 20 GB again,
hard-link the data: `NPX_LINK_FROM=seg_g100 python scripts/make_segment.py 806 1744 seg_g50 mynewrec_g50`.

Those files are the complete input to MATLAB, and you can write them from any source - see
**Input file formats** below. Check them first:

```bash
python scripts/check_segment.py "<ap folder>/seg_g100"
```

### 4. Interpolate (MATLAB)

The ADC phase shift is applied to the AP data before kriging,
automatically (see **ADC phase shift** below).

```bash
matlab -batch "addpath('C:/path/to/npix_motion_pipeline/matlab'); run_interp('C:/.../seg_g100', 'C:/path/to/InterpolationAfterDREDge')"
```

All parameters come from `interp_params.json` in the segment folder, so no MATLAB file needs editing
per run. Output: `<segment>/interpseg_g100/interpseg_g100.dat`, byte-for-byte the same size as the
input, with any sync channel copied through untouched.

### 5. Check it against the spikes

```bash
NPX_REC=mynewrec python -u scripts/template_tracking.py seg_g100 810 1740
```

Detects units on the before file, builds each unit's mean waveform, cross-correlates it across depth
and time in **both** files, and tracks the best-correlated depth over time. Writes to `results/`:

- `mynewrec_seg_g100_examples.png` - raw traces with detected spikes, and the space x time
  correlation map per unit, before vs after
- `mynewrec_seg_g100_tracking.png` - tracked position of every unit, before vs after, plus a paired
  box plot of per-unit variability
- `mynewrec_seg_g100_tracking.txt` / `.npz` - the numbers

**How to read it.** The key line is the regression of each unit's tracked position on the
displacement that was applied:

| after-correction slope / r | meaning |
|---|---|
| ~0 | motion removed - this is what you want |
| negative | **over-corrected**: DREDge's amplitude is too big, lower the gain |
| positive (and lower than before) | under-corrected: raise the gain |

The sign is trustworthy even though the absolute size is not: tracker sensitivity depends on how
sharp each unit's template is, so **never compare SD between two files whose residual motion
differs** - compare signs, or compare files that differ in only one thing.

### 6. Calibrate the gain, then produce the final file

If the slope after step 5 is not ~0, run steps 2-5 again at a different gain (e.g. `"0,0.1,0.5"`,
segment `seg_g50`), then:

```bash
NPX_REC=mynewrec python scripts/bracket_gain.py seg_g100 1.0 seg_g50 0.5
```

The after-correction slope is linear in the gain, so two runs give the crossing - the fraction of
DREDge's displacement your units actually move. Put it in `recordings.json` as `slow_gain`, derive
the motion again with that gain over the full valid range, interpolate, and re-run tracking to
confirm slope and r are ~0. On msit004: gain 1.0 gave slope -0.45, gain 0.49 gave +0.30, the
bracket predicted 0.6, and the run at 0.6 gave -0.05 with per-unit SD 23.5 -> 13.8 µm (5/5 units).

### 7. Sort

Feed `interp<segment>.dat` to Kilosort with `NchanTOT` = `nch` (sync included) and `Nchan` =
`nchan_neural`, using the same channel map. The file has identical geometry to the input.

## Input file formats

[example/](example/) holds real msit004 files, trimmed to 20 s so they can be opened in a text
editor, with [example/README.md](example/README.md) describing every field:

- `example/motion_estimate/` - what DREDge writes: `<tag>_displacement.csv` (µm, **rows = depth
  windows, columns = time at 250 Hz**), `<tag>_window_centers.csv` (depth of each row),
  `<tag>_t0.txt`, `<tag>_config.json`.
- `example/segment_seg_example/` - what MATLAB reads: the sliced displacement, window centres,
  `ref_displacement.csv`, `segment_t0.txt`, `chanMap.mat`, `inter_sample_shift.csv` and
  `interp_params.json`. The raw
  `continuous.dat` (int16, channel-interleaved, width = `nch` including sync, no header) is the only
  file not included, for size.

`scripts/check_segment.py` validates a folder against those rules - row counts matching between the
three CSVs, data length covering the motion trace, file size divisible by the channel count, channel
map and window depths in the same coordinates - and its docstring is the format reference.

**Nothing requires DREDge.** To apply a displacement from another estimator, write those CSVs in the
shapes above (one row = a rigid estimate), run `check_segment.py`, then `run_interp.m`, then verify
with `template_tracking.py` exactly as below.

## ADC phase shift

Neuropixels digitises channels at staggered times within one sample period (up to 28 us on
NP1.0). Kriging mixes neighbouring channels, so this has to be corrected **before**
interpolation - mixing time-skewed copies distorts the waveform and nothing downstream,
Kilosort included, can undo it. Kilosort does not do this correction, so if the pipeline does
not, it never happens.

It is **on by default** and needs nothing from you: `make_segment.py` writes
`inter_sample_shift.csv` from the recording's metadata, and MATLAB applies it per block before
the kriging (`matlab/apply_frequency_shift.m`, `matlab/phase_shift_reader.m`). The
implementation is a port of `spikeinterface.preprocessing.phase_shift`, matching it to 6e-5 ADC
counts on real data - verify on your own recording with:

```bash
python scripts/check_phase_shift.py "<ap folder>/seg_g100"
```

Motion estimation already used phase-shifted data (`dredge_recording.py` phase-shifts the LFP
before DREDge). To interpolate without the correction, set `"phase_shift": false` in the
segment's `interp_params.json`; a segment missing the CSV errors out rather than skipping it
quietly.

## Cardiac (1-2 Hz) motion

**Default: off** (`"0,0.1,<gain>"` selects slow drift only). On msit004, DREDge-LFP reported ~50 µm
of cardiac motion, but units disagreed with it and with each other: three units moved ~50 µm in
phase with it, six barely moved at all, and neighbours 100 µm apart fell in different groups. No
smooth displacement field can fix that, and in practice every cardiac gain we tried made unit
cardiac motion worse or left it unchanged (x1.0: ~2.9x worse; x0.6: 1.21x; x0.28: 1.19-1.35x,
p = 0.71). Finer DREDge windows (200 µm) did not resolve it either.

If you want to test it on your own recording before deciding, that is what `scripts/cardiac_calib.py`,
`scripts/driftfree_cardiac.py` and `scripts/card_by_depth.py` do; the method and its pitfalls are in
[FINDINGS_msit001_msit004.md](FINDINGS_msit001_msit004.md). The short version: compare two files that
differ **only** in the cardiac correction, and look at per-unit amplitudes rather than a pooled gain.

## Gotchas that cost us real time

- **Channel count.** Reading a 385-channel file as 384 shears every sample and corrupts the sync
  channel. `nch` is the file's width; `nchan_neural` is what gets interpolated. PXI probe streams
  have a sync channel, OneBox probe streams do not.
- **Channel map geometry.** Our probe (NP1016, 1.0 Sapiens) has *aligned* columns; the default
  Kilosort map is staggered and wrong on half the sites. Check `settings.xml`.
- **Valid range.** Exclude insertion, driving to target and withdrawal. Motion there is huge and
  meaningless, and it distorts the estimate for the rest.
- **`estimate_motion` window override** (see step 1).
- **OneBox clocks** do not start at 0; SpikeInterface motion objects carry that offset while the
  displacement CSVs here are indexed from the start of the slice.
- **Disk**: use `NPX_LINK_FROM` for repeat segments over the same range, and delete superseded
  `interp*` outputs - their analysis caches in `results/` are small and are all you need later.

## What is in here

```
recordings.json        per-recording config - the only file you normally edit
example/               real input files, trimmed, with a README describing every format
scripts/               check_setup, dredge_recording, derive_motion, make_segment, check_segment,
                       template_tracking, bracket_gain  (the pipeline, in order)
                       check_phase_shift  (verifies the MATLAB ADC shift against spikeinterface)
                       cardiac_calib, driftfree_cardiac, card_by_depth, calibrate_gains,
                       score_driftfree, v2_eval, summary_figure  (cardiac + reporting)
                       injection_test, synthetic_kriging, lfp_scale_test  (validation of the
                       method itself: known motion in, measured motion out)
matlab/                run_interp.m driver, the ADC phase shift (apply_frequency_shift.m,
                       phase_shift_reader.m) + the 4 modified InterpolationAfterDREDge files
chanmaps/              channel maps for the NP1016 probes used here
results/               everything the scripts write
CREDITS.md             provenance: upstream repos, licenses, and what originated here
FINDINGS_*.md          full record of the validation, including approaches that failed
summary_interpolation.png, card_by_depth.png   headline figures from those recordings
```

Every script's docstring explains what it does and why; `python scripts/<name>.py` with no arguments
prints usage.
