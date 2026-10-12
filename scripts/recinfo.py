"""Per-recording paths and geometry, read from recordings.json. Shared by every script here.

Select the recording with the NPX_REC environment variable:

    set NPX_REC=msit004
    python template_tracking.py seg_check 1200 1330

Add a recording by copying an entry in ../recordings.json (see README.md, "Adding a recording").
Point NPX_RECORDINGS at a different file to keep your configs outside this folder.

Fields: ap (folder holding continuous.dat), lfp_base + lfp_stream (what SpikeInterface opens for
the LFP), motion (where motion estimates are written), chanmap (.mat with xcoords/ycoords, relative
paths resolve against the package root), nch (channels in the AP FILE, including sync),
nchan_neural (channels to interpolate), valid [t0, t1] in file-relative seconds, max_depth_um
(ignore sites above this, e.g. in saline; null = use all), motion_lowpass_hz, baseline_motion_tag,
slow_gain / cardiac_gain (per-recording gains, filled in after calibration).
nframes is derived from the file size unless the entry states it.
"""
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
_CFG = Path(os.environ.get("NPX_RECORDINGS", ROOT / "recordings.json"))
# every figure, cache and table the scripts write lands here (override with NPX_RESULTS)
RESULTS = Path(os.environ.get("NPX_RESULTS", ROOT / "results"))
RESULTS.mkdir(parents=True, exist_ok=True)
_CACHE = {}


def _load():
    if not _CACHE:
        if not _CFG.exists():
            raise SystemExit("no recordings file at %s (set NPX_RECORDINGS)" % _CFG)
        _CACHE.update(json.loads(_CFG.read_text()))
    return _CACHE


def key():
    k = os.environ.get("NPX_REC")
    if not k:
        names = [n for n in _load() if not n.startswith("_")]
        raise SystemExit("set NPX_REC to one of: %s" % ", ".join(names))
    return k


def _path(p):
    p = Path(p)
    return p if p.is_absolute() else (ROOT / p)


def get(k=None):
    k = k or key()
    cfg = _load()
    if k not in cfg:
        raise SystemExit("recording '%s' is not in %s" % (k, _CFG))
    r = dict(cfg[k])
    for f in ("ap", "motion", "chanmap", "lfp_base"):
        if r.get(f):
            r[f] = _path(r[f])
    if not r.get("nframes"):
        dat = r["ap"] / "continuous.dat"
        if not dat.exists():
            raise SystemExit("%s not found - check 'ap' for recording '%s'" % (dat, k))
        size = dat.stat().st_size
        if size % (2 * r["nch"]):
            raise SystemExit("%s is %d bytes, not a whole number of %d-channel int16 samples - "
                             "check 'nch' (sync channel included?)" % (dat, size, r["nch"]))
        r["nframes"] = size // (2 * r["nch"])
    r.setdefault("valid", [0.0, r["nframes"] / 30000.0])
    return r


def chanmap_yx(k=None):
    """Return (ycoords, xcoords) for the recording's probe, from either .mat format."""
    import numpy as np
    p = str(get(k)["chanmap"])
    try:
        import h5py
        with h5py.File(p, "r") as f:
            return np.array(f["ycoords"]).ravel(), np.array(f["xcoords"]).ravel()
    except (OSError, ImportError):
        from scipy.io import loadmat
        m = loadmat(p)
        return np.asarray(m["ycoords"]).ravel(), np.asarray(m["xcoords"]).ravel()


def row_matrix(k=None):
    """Collapse the probe's column pairs to one value per depth row.

    NP1.0-style probes put 2 sites at each y, offset in x. Cross-correlation over
    CHANNEL INDEX then mixes a depth shift with a column swap: an odd index shift lands
    the template on the other column, which carries a different amplitude for the same
    unit. Tracking in row space (uniform 20 um pitch) removes that confound. This is the
    same reduction average_across_direction performs before DREDge.

    Returns (y_rows, R) where x_rows = x @ R.
    """
    import numpy as np
    yv, _ = chanmap_yx(k)
    n = get(k)["nchan_neural"]
    yv = yv[:n]
    uy = np.unique(yv)
    R = np.zeros((n, len(uy)), dtype="float32")
    for j, u in enumerate(uy):
        idx = np.where(yv == u)[0]
        R[idx, j] = 1.0 / len(idx)
    return uy, R
