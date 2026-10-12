"""Check a recordings.json entry before running anything heavy.

    set NPX_REC=<name>
    python check_setup.py

Verifies the AP file exists and its size matches nch (the commonest setup error: reading a
385-channel file as 384 shears every sample and corrupts the sync channel), that the channel map
matches nchan_neural, and that the LFP stream opens. Prints the probe geometry so you can confirm
the channel map is the right one for the probe - compare the column x positions with settings.xml.
"""
import os
import sys

import numpy as np

import recinfo


def main():
    k = recinfo.key()
    r = recinfo.get(k)
    print("recording %s" % k)
    dat = r["ap"] / "continuous.dat"
    print("  AP file      %s" % dat)
    print("  %d channels in file, %d neural, %d frames = %.1f s at 30 kHz"
          % (r["nch"], r["nchan_neural"], r["nframes"], r["nframes"] / 30000.0))
    print("  valid range  %.1f - %.1f s" % tuple(r["valid"]))
    y, x = recinfo.chanmap_yx(k)
    print("  chanMap      %s" % r["chanmap"])
    print("               %d sites, y %.0f-%.0f um, pitch %s, columns at x = %s"
          % (len(y), y.min(), y.max(), np.unique(np.diff(np.unique(y))), np.unique(x)))
    if len(y) < r["nchan_neural"]:
        print("  !! channel map has fewer sites (%d) than nchan_neural (%d)" % (len(y), r["nchan_neural"]))
        return 1
    yr, _ = recinfo.row_matrix(k)
    print("               %d depth rows after collapsing columns" % len(yr))
    if r.get("motion"):
        print("  motion dir   %s%s" % (r["motion"], "" if r["motion"].exists() else "  (will be created)"))
    try:
        import spikeinterface.full as si
        rec = si.read_openephys(r["lfp_base"], stream_name=r["lfp_stream"], load_sync_channel=False)
        ly = rec.get_channel_locations()[:, 1]
        print("  LFP stream   %s: %d ch at %.0f Hz, %.0f s, y %.0f-%.0f um"
              % (r["lfp_stream"], rec.get_num_channels(), rec.get_sampling_frequency(),
                 rec.get_num_frames() / rec.get_sampling_frequency(), ly.min(), ly.max()))
        if abs(np.median(np.diff(np.unique(ly))) - np.median(np.diff(np.unique(y)))) > 1e-6:
            print("  !! LFP row pitch differs from the channel map - displacements would be mis-scaled")
    except Exception as exc:
        print("  !! could not open the LFP stream: %r" % exc)
        return 1
    print("OK - ready to run dredge_recording.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
