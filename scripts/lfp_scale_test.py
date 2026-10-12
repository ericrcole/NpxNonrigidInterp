"""Does DREDge-LFP report the right AMPLITUDE of a known shift? (msit004, same pipeline + config)

    python lfp_scale_test.py

Takes the real preprocessed LFP (bandpass 0.5-250, phase_shift, 250 Hz, averaged across columns)
over 900-1200 s, adds a KNOWN extra displacement by resampling each frame along depth
    extra(t) = 40 um * sin(2 pi t / 60 s) + 10 um * sin(2 pi 1.5 Hz t)
and re-runs estimate_motion with the nonrigid_800 configuration on both the unmodified and the
shifted data. The recovered extra motion is the difference of the two estimates, so the
recording's own motion cancels. Recovered/injected gain ~1 means DREDge-LFP measures shift
amplitude correctly and the ~0.5x (slow) / ~0.27x (cardiac) unit-vs-DREDge ratios are not a
scale artifact of the estimator path.
"""
import numpy as np
import spikeinterface.core as sc
from scipy.interpolate import interp1d
from scipy.signal import butter, sosfiltfilt
from spikeinterface.sortingcomponents.motion import estimate_motion

import dredge_recording as dr

FS = dr.LFP_FS
T0, T1 = 900.0, 1200.0


def main():
    rec, _ = dr.build("msit004", 805.0, 1745.0, 0.5)
    f0, f1 = int((T0 - 805.0) * FS), int((T1 - 805.0) * FS)
    X = rec.get_traces(start_frame=f0, end_frame=f1, return_scaled=False).astype(np.float32)  # (t, rows)
    loc = rec.get_channel_locations()
    y = loc[:, 1]
    order = np.argsort(y)
    X, y, loc = X[:, order], y[order], loc[order]
    t = np.arange(X.shape[0]) / FS
    extra = 40.0 * np.sin(2 * np.pi * t / 60.0) + 10.0 * np.sin(2 * np.pi * 1.5 * t)
    # content at depth y moves to y + extra: sample the original at y - extra
    Xs = np.empty_like(X)
    for i in range(X.shape[0]):
        Xs[i] = interp1d(y, X[i], kind="cubic", bounds_error=False, fill_value="extrapolate")(y - extra[i])
    cfg = dict(dr.CONFIGS["nonrigid_800"])
    est = []
    for name, data in (("orig", X), ("shifted", Xs)):
        r = sc.NumpyRecording(data, sampling_frequency=FS)
        r.set_channel_locations(loc)
        m = estimate_motion(r, method="dredge_lfp", progress_bar=False, **cfg)
        est.append(np.asarray(m.displacement[0]))
        print("%s: displacement %s, windows %s" % (name, est[-1].shape, np.round(m.spatial_bins_um, 0)), flush=True)
    d = est[1] - est[0]                          # (t, windows)
    n = min(len(d), len(extra))
    d, ex = d[:n], extra[:n]
    ex = ex - ex.mean()
    d = d - d.mean(0)
    slow = butter(3, 0.1 / (FS / 2), btype="low", output="sos")
    card = butter(3, [0.8 / (FS / 2), 5.0 / (FS / 2)], btype="band", output="sos")
    for label, sos in (("slow (<0.1 Hz)", slow), ("cardiac (0.8-5 Hz)", card)):
        e = sosfiltfilt(sos, ex)
        for w in range(d.shape[1]):
            r_ = sosfiltfilt(sos, d[:, w])
            s = slice(int(10 * FS), n - int(10 * FS))
            gain = np.dot(r_[s], e[s]) / np.dot(e[s], e[s])
            print("%-20s window %d: recovered gain %.2f  r %.2f  (sign + means DREDge reports content moving up)"
                  % (label, w, gain, np.corrcoef(r_[s], e[s])[0, 1]))


if __name__ == "__main__":
    main()
