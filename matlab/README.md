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

Keep a copy of the originals first (or work on a branch), then:

```bash
matlab -batch "addpath('<this folder>'); run_interp('<segment folder>', '<InterpolationAfterDREDge folder>')"
```

The second argument can be omitted if the `INTERP_REPO` environment variable is set. Every other
parameter comes from `interp_params.json`, which `make_segment.py` writes into the segment folder.

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

**`skipDriftMap`.** Motion comes from DREDge, so the Kilosort drift-map step is skipped and `rez` is
built directly from the channel map.
