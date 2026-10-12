"""Validate a prepared segment folder - every input the MATLAB interpolation reads.

    python check_segment.py <segment folder>

e.g.  python check_segment.py example/segment_seg_example      (no continuous.dat: inputs only)
      NPX_REC=msit004 python check_segment.py "D:/.../OneBox-122.ProbeA-AP/seg_g100"

A segment folder is what make_segment.py writes, and is the complete input to
`matlab/run_interp.m`. Run this after preparing a segment by hand, or after writing a motion
estimate with something other than dredge_recording.py.

The files, and what each one must contain
-----------------------------------------
continuous.dat
    The raw AP data: int16, little-endian, CHANNEL-INTERLEAVED (sample 0 ch 0, sample 0 ch 1, ...),
    exactly as Open Ephys writes it. Width = `nch` in recordings.json, INCLUDING the sync channel
    if the stream has one (385 on a PXI probe stream, 384 on a OneBox probe stream). Reading a
    385-channel file as 384 shears every sample; there is no header to catch the mistake, so the
    file size must be divisible by 2 * nch. The interpolation writes an output of identical size,
    copying any sync channel through unchanged.

nonrigid_displacement.csv
    The motion to apply, in MICROMETRES, as plain comma-separated text.
    ROWS = depth windows (1 row = a rigid estimate), COLUMNS = time at 250 Hz
    (`motion_fs` in interp_params.json). Column 0 is the time in segment_t0.txt. Positive values
    mean the tissue sits higher up the probe; the correction applied is
    (displacement - ref_displacement), so the absolute zero does not matter.
    Length must cover the segment: columns >= segment duration * motion_fs.

dredge_window_centers.csv
    One value per ROW of the displacement, the depth in micrometres of that window's centre, in the
    same coordinates as the channel map's ycoords. Shifts are interpolated across these with pchip
    and clamped outside the range.

ref_displacement.csv
    One value per row: the reference position subtracted before applying. make_segment.py takes it
    over the WHOLE recording, so interpolating a segment matches interpolating the whole file. If
    you prepare a segment by hand, a per-window median over the full estimate is the right value.

segment_t0.txt
    One number: the segment's start time in seconds, relative to the start of the raw recording.
    Used to line the displacement up with spike times during verification.

chanMap.mat
    Kilosort-style channel map with xcoords/ycoords (and chanMap/connected). Must describe the
    probe that recorded this file - the kriging weights come from these coordinates.

interp_params.json
    motion_fs (sample rate of the displacement, 250), batch_samples (interpolation batch length in
    30 kHz samples, 128), nchan_raw (= nch, sync included), m_tol_um (rebuild the kriging operator
    when the shift moves this far, 0.5), out_name (output folder and .dat name).
"""
import json
import sys
from pathlib import Path

import numpy as np

FS = 30000.0


def main(seg):
    seg = Path(seg)
    if not seg.is_dir():
        raise SystemExit("not a folder: %s" % seg)
    bad = 0
    print("segment %s" % seg)
    p = json.loads((seg / "interp_params.json").read_text())
    print("  interp_params.json  %s" % p)
    t0 = float((seg / "segment_t0.txt").read_text().strip())
    D = np.loadtxt(seg / "nonrigid_displacement.csv", delimiter=",", ndmin=2)
    wc = np.atleast_1d(np.loadtxt(seg / "dredge_window_centers.csv", delimiter=","))
    ref = np.atleast_1d(np.loadtxt(seg / "ref_displacement.csv", delimiter=","))
    dur = D.shape[1] / p["motion_fs"]
    print("  displacement        %d window(s) x %d samples = %.1f s at %d Hz, starting at t = %.1f s"
          % (D.shape[0], D.shape[1], dur, p["motion_fs"], t0))
    print("                      range %.1f to %.1f um, per-window SD %s um"
          % (D.min(), D.max(), np.round(D.std(1), 1)))
    print("  window centres      %s um" % np.round(wc, 0))
    print("  reference           %s um" % np.round(ref, 1))
    for name, arr in (("dredge_window_centers.csv", wc), ("ref_displacement.csv", ref)):
        if len(arr) != D.shape[0]:
            print("  !! %s has %d values but the displacement has %d rows" % (name, len(arr), D.shape[0]))
            bad += 1
    if D.shape[0] > 1:
        strain = np.abs(np.diff(D, axis=0) / np.median(np.diff(wc))).max()
        print("  max strain          %.3f between neighbouring windows%s"
              % (strain, "  (>0.5 is implausible - check the window size)" if strain > 0.5 else ""))
    dat = seg / "continuous.dat"
    if dat.exists():
        size = dat.stat().st_size
        nfr = size / (2 * p["nchan_raw"])
        print("  continuous.dat      %d bytes = %.0f frames x %d ch = %.1f s"
              % (size, nfr, p["nchan_raw"], nfr / FS))
        if size % (2 * p["nchan_raw"]):
            print("  !! not a whole number of %d-channel int16 samples - wrong channel count" % p["nchan_raw"])
            bad += 1
        elif nfr / FS > dur + 1e-6:
            print("  !! data is %.1f s but the motion trace covers only %.1f s" % (nfr / FS, dur))
            bad += 1
    else:
        print("  continuous.dat      absent (inputs-only example; the real thing is the raw AP data)")
    try:
        from scipy.io import loadmat
        m = loadmat(seg / "chanMap.mat")
        y, x = np.asarray(m["ycoords"]).ravel(), np.asarray(m["xcoords"]).ravel()
    except Exception:
        import h5py
        with h5py.File(seg / "chanMap.mat", "r") as f:
            y, x = np.array(f["ycoords"]).ravel(), np.array(f["xcoords"]).ravel()
    print("  chanMap.mat         %d sites, y %.0f-%.0f um (pitch %s), x %s"
          % (len(y), y.min(), y.max(), np.unique(np.diff(np.unique(y))), np.unique(x)))
    if wc.min() < y.min() - 500 or wc.max() > y.max() + 500:
        print("  !! window centres and channel-map depths disagree - different coordinate systems?")
        bad += 1
    print("FAILED %d check(s)" % bad if bad else "OK - ready for matlab/run_interp.m")
    return 1 if bad else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    sys.exit(main(sys.argv[1]))
