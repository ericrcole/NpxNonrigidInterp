# MATLAB side

`run_interp.m` is the driver you call. The other four files are **modified versions of files in
[InterpolationAfterDREDge](https://github.com/williamunoz/InterpolationAfterDREDge)** by Richard
Hardstone and William Munoz, used and redistributed here under its CC0 1.0 public domain
dedication; each carries a header saying so. They adapt Kilosort 2.5, and the original per-file
credits are intact. Full provenance: [../CREDITS.md](../CREDITS.md).

These four must replace the originals in your clone:

| file here | replaces in the repo |
|---|---|
| `shift_batch_on_disk2_dredge.m` | `Scripts/AP/Dredge/shift_batch_on_disk2_dredge.m` |
| `applyDredge.m` | `Scripts/AP/Dredge/applyDredge.m` |
| `RealignRecordingUsingDredgeAndKilosort2point5.m` | `Scripts/AP/Dredge/RealignRecordingUsingDredgeAndKilosort2point5.m` |
| `Kilosort25Configuration_Dredge.m` | `Scripts/AP/Kilosort25_Configs/Kilosort25Configuration_Dredge.m` |

**You do not have to copy them in.** `run_interp.m` puts this folder at the front of the
MATLAB path after it adds the repo, so these copies shadow the repo's originals and the
phase-shift helpers resolve. Copy them over the originals only if you call the repo's
entry points directly instead of going through `run_interp.m`; keep a copy of the originals
first (or work on a branch).

That shadowing matters: `run_interp.m` calls `restoredefaultpath`, so anything the caller
added is dropped. Before this was handled, a run silently used the repo's unmodified
`applyDredge` and `shift_batch_on_disk2_dredge` - same arguments, no error, no phase shift -
and produced a byte-identical file to an uncorrected run.

To run:

```bash
matlab -batch "addpath('<this folder>'); run_interp('<segment folder>', '<InterpolationAfterDREDge folder>')"
```

The second argument can be omitted if the `INTERP_REPO` environment variable is set. Every other
parameter comes from `interp_params.json`, which `make_segment.py` writes into the segment folder.

`run_interp.m` calls two helpers that are new here, not from upstream: `apply_frequency_shift.m`
and `phase_shift_reader.m` (the ADC phase shift, below). `dump_phase_shift.m` exists only for the
equivalence test.

## What was changed, and why

**File geometry and the sync channel.** The raw file's width (`NchanTOT`, e.g. 385 = 384 neural +
sync) is now kept separate from the number of channels that get interpolated (`Nchan`, 384).
Reading a 385-channel file as 384 shears every sample and destroys the sync channel. The sync
column is now copied through verbatim, and writes are clamped so the output can never grow past the
input's length.

**Streaming instead of copying.** Input and output are opened as streams (`fidIn`/`fidOut`) and the
frame count is taken from the *source* file, so no multi-GB copy is made first and a half-written
output cannot silently define the length.

**DC offset mixing (a real artifact, now fixed).** Kriging mixes channels with time-varying weights.
Because each channel carries its own DC offset, that mixing injected power below 150 Hz - the
standard deviation of the <150 Hz band rose from 30 (raw) to 108. Offsets are now removed before
interpolation and the output channel's own offset is added back (`removeOffset`, default true);
<150 Hz std drops to 29 while the spike band changes by 0.9% RMS.

**Nonrigid shifts.** Per-depth shifts are interpolated onto channels with `pchip` and clamped at the
edges; a rigid (scalar) estimate is expanded to all channels rather than silently mis-tapering.
Batches are crossfaded with a raised cosine, as before.

**Speed.** The kriging matrix `Kxx` depends only on probe geometry, so it is computed once rather
than per batch; the full operator is cached and rebuilt only when the shift moves by more than
`mTolUm` (0.5 µm). A 15-minute 384-channel segment takes ~8-15 minutes.

**Reference displacement.** `refDisplacement` is passed in explicitly (from `ref_displacement.csv`,
computed over the whole recording) instead of taken as the median of whatever segment is being
processed, so interpolating a segment gives the same result as interpolating the whole file.

**Neuropixels ADC phase shift, applied before kriging (ON by default).** The probe multiplexes
channels onto a few ADCs, so within one sample period each channel is digitised at a slightly
different time - `inter_sample_shift`, 0 to ~0.85 samples (28 us) on NP1.0. Kriging builds every
output channel from a weighted sum of its neighbours, so without this it sums waveforms skewed in
time by a fraction of a sample, and because the weights follow the motion the distortion varies
over time. It cannot be fixed afterwards: the mixing has already happened. Kilosort does not do
this correction either, so if it is not done here it never happens.

`apply_frequency_shift.m` is a port of `spikeinterface.preprocessing.phase_shift`'s
`apply_frequency_shift` (which follows IBL and SpikeGLX): rFFT, multiply by
`exp(-1j * 2*pi*rfftfreq(N) * shift)`, inverse. `phase_shift_reader.m` reproduces the same chunk
handling - 40 ms margin, zero padding past the file edges, raised-cosine taper on each margin -
but reads one large block at a time (60000 frames by default) and serves the 128-sample batches
from it, so the margin is not recomputed per batch. Measured cost on a 120 s segment of
msit004: 23 s without, 49 s with, so budget roughly double the interpolation time.

What it changes in the output: samples move by a fraction of a sample, which is a large
change in the sample values themselves (23.6% of signal SD on msit004) but almost none in
spike amplitude (mean-waveform peak-to-peak -0.25% median, 10th-90th percentile -0.89% to
+0.08%). The point is the inter-channel TIMING that kriging, CAR and template matching all
depend on, not amplitude.

`scripts/check_phase_shift.py` verifies the port against SpikeInterface on real data. On msit004
the two agree to **6e-5 ADC counts** (1.7e-7 of the correction's own size), far below the
+/-0.5 count rounding of SpikeInterface's own int16 output.

Per-channel delays come from `inter_sample_shift.csv` in the segment folder, which
`make_segment.py` writes from the recording's metadata. Turn the correction off only
deliberately, with `"phase_shift": false` in `interp_params.json`; a segment without the CSV
fails with a message saying so rather than silently skipping it.

**`skipDriftMap`.** Motion comes from DREDge, so the Kilosort drift-map step is skipped and `rez` is
built directly from the channel map.
