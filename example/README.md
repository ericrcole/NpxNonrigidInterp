# Example inputs

Real files from msit004 (NPIX_MSIT_004, Neuropixels 1.0 Sapiens on a OneBox), trimmed to the first
20 s so they are small enough to read in a text editor. Use them as the reference for what every
input must look like. To check your own files against the same rules:

```bash
python scripts/check_segment.py example/segment_seg_example
python scripts/check_segment.py "D:/.../OneBox-122.ProbeA-AP/my_segment"
```

Everything is plain comma-separated text except `chanMap.mat` (MATLAB) and `continuous.dat` (raw
int16). Nothing has a header row.

## `motion_estimate/` - what DREDge writes

`dredge_recording.py` writes these per estimate; `derive_motion.py` writes the same four files for
each band/gain variant it derives. The name before the suffix is the "tag" other scripts refer to.

| file | contents |
|---|---|
| `<tag>_displacement.csv` | displacement in µm. **Rows = depth windows, columns = time at 250 Hz.** Here 4 × 5000 = 4 windows × 20 s |
| `<tag>_window_centers.csv` | depth in µm of each row: `710, 1510, 2310, 3110` |
| `<tag>_t0.txt` | file-relative time of column 0, in seconds: `805.0` (this recording's valid range starts there) |
| `<tag>_config.json` | provenance: `{"source": "msit004_nonrigid_800_t805-1745", "bands": [[0.0, 0.1, 0.6]]}` = slow drift below 0.1 Hz, scaled to 0.6 |

The example is `msit004_n800_v3_s60`, the final msit004 estimate: slow drift only, at the calibrated
gain of 0.6. A rigid estimate is a single row.

## `segment_seg_example/` - what MATLAB reads

`make_segment.py` writes one of these per interpolation run. The only file not included here is
`continuous.dat` (2.8 GB for this segment); everything else is byte-for-byte what was used.

| file | contents |
|---|---|
| `continuous.dat` | **(not included)** raw AP data: int16, little-endian, channel-interleaved, no header. Width = `nchan_raw`, **sync channel included**. Its size must be divisible by `2 * nchan_raw` |
| `nonrigid_displacement.csv` | the estimate above, sliced to this segment: 4 × 5000, µm, rows = windows, columns = time at 250 Hz |
| `dredge_window_centers.csv` | depth of each row, µm - same coordinates as the channel map's ycoords |
| `ref_displacement.csv` | one value per row, subtracted before applying: `-31.1, -28.5, -21.5, -19.6`. Taken over the WHOLE recording, so a segment run equals the matching batches of a full-length run |
| `segment_t0.txt` | `1000.0` - the segment's start, in seconds from the start of the raw recording |
| `chanMap.mat` | Kilosort-style map with `xcoords`/`ycoords`; here 384 sites, y 0-3820 µm at 20 µm pitch, columns at x = 0 and 32 |
| `interp_params.json` | `motion_fs` 250 (rate of the displacement), `batch_samples` 128 (interpolation batch, in 30 kHz samples), `nchan_raw` 384 (file width), `m_tol_um` 0.5 (rebuild the kriging operator when the shift moves this far), `out_name` (output folder and `.dat` name) |

Sign convention: positive displacement means the tissue sits higher up the probe. Only
`displacement - ref_displacement` is applied, so the absolute zero is arbitrary - but the units are
µm and must match the channel map's.

## Using your own motion estimate

Nothing here requires DREDge. To apply a displacement from another source, write
`nonrigid_displacement.csv`, `dredge_window_centers.csv` and `ref_displacement.csv` in the shapes
above (rows = depth windows at 250 Hz; one row = rigid), drop in `chanMap.mat`, `segment_t0.txt`
and `interp_params.json`, run `check_segment.py`, then `matlab/run_interp.m`.

Verify the result the same way regardless of where the estimate came from: `template_tracking.py`,
then the slope/r rule in the main README. The pipeline applies whatever it is given; only that check
tells you whether it helped.
