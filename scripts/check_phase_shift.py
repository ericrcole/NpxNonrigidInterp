"""Verify the MATLAB ADC phase shift against spikeinterface.phase_shift on real data.

    python check_phase_shift.py <segment folder> [f0_frame] [n_frames]

matlab/apply_frequency_shift.m is a port of spikeinterface's apply_frequency_shift, and
matlab/phase_shift_reader.m reproduces its chunk handling (40 ms margin, raised-cosine taper on
the margin, zero padding at file edges). This runs both on the same window of the same file and
reports the difference, in ADC counts, against the size of the correction itself.

The two are not expected to be bit-identical: both split the recording into chunks and the chunk
boundaries differ, which is exactly the small chunking error the 40 ms margin exists to bound.
A pass means the residual is a tiny fraction of the correction, and far below one ADC count.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

MATLAB = r"C:/Program Files/MATLAB/R2024b/bin/matlab.exe"
HERE = Path(__file__).resolve().parent
FS = 30000.0


def main(seg, f0=300000, n=60000):
    seg = Path(seg)
    p = json.loads((seg / "interp_params.json").read_text())
    iss = np.loadtxt(seg / "inter_sample_shift.csv", delimiter=",")
    nraw = p["nchan_raw"]
    nchan = len(iss)

    import spikeinterface.full as si
    rec = si.read_binary(seg / "continuous.dat", sampling_frequency=FS, dtype="int16",
                         num_channels=nraw, is_filtered=False)
    rec = rec.select_channels(rec.channel_ids[:nchan])
    rec.set_property("inter_sample_shift", iss)
    # Work in float64: spikeinterface's default casts the result back to int16, and that
    # rounding (+/-0.5 counts, SD 0.289) would dominate and hide the real agreement. The
    # recording itself must be float too, or its margin taper fails on an int16 buffer.
    recf = si.astype(rec, dtype="float64")
    ps = si.phase_shift(recf, margin_ms=p["phase_shift_margin_ms"], inter_sample_shift=iss)
    py = ps.get_traces(start_frame=f0, end_frame=f0 + n).astype(np.float64)
    raw = rec.get_traces(start_frame=f0, end_frame=f0 + n).astype(np.float64)

    out = Path(tempfile.gettempdir()) / "ps_matlab.bin"
    cmd = [MATLAB, "-batch", "addpath('%s'); dump_phase_shift('%s', %d, %d, '%s')"
           % ((HERE.parent / "matlab").as_posix(), seg.as_posix(), f0, n, out.as_posix())]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not out.exists():
        raise SystemExit("MATLAB dump failed:\n%s\n%s" % (r.stdout[-2000:], r.stderr[-2000:]))
    mat = np.fromfile(out, dtype="float32").reshape(n, nchan).astype(np.float64)

    corr = py - raw                       # what the correction itself changes
    diff = mat - py                       # how far MATLAB is from spikeinterface
    print("window %d-%d, %d channels, max channel skew %.1f us"
          % (f0, f0 + n, nchan, 1e6 * (iss.max() - iss.min()) / FS))
    print("  raw signal           SD %8.3f counts" % raw.std())
    print("  correction (SI-raw)  SD %8.3f counts, max |.| %8.3f" % (corr.std(), np.abs(corr).max()))
    print("  MATLAB - SI          SD %8.5f counts, max |.| %8.5f" % (diff.std(), np.abs(diff).max()))
    print("  residual / correction  %.2e (SD ratio)" % (diff.std() / max(corr.std(), 1e-12)))
    worst = int(np.argmax(np.abs(diff).max(0)))
    print("  worst channel %d (shift %.3f samples): max |diff| %.5f counts"
          % (worst, iss[worst], np.abs(diff[:, worst]).max()))
    print("  (for scale: int16 rounding alone would give SD 0.289 counts)")
    ok = diff.std() < 1e-3 * corr.std() and np.abs(diff).max() < 0.01
    print("PASS - MATLAB reproduces spikeinterface.phase_shift" if ok else
          "FAIL - difference is larger than chunking error should allow")
    return 0 if ok else 1


if __name__ == "__main__":
    a = sys.argv[1:]
    sys.exit(main(a[0], *(int(x) for x in a[1:])))
