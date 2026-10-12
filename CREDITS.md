# Credits and provenance

This pipeline is a thin layer over other people's work. What is original here is the Python
estimation/verification tooling, the gain-calibration procedure, and fixes to the MATLAB
interpolation; everything underneath was written by the people below.

## The MATLAB interpolation

**InterpolationAfterDREDge** — Richard Hardstone and William Munoz
<https://github.com/williamunoz/InterpolationAfterDREDge> — released under **CC0 1.0 Universal**
(public domain dedication), which is what allows the modified copies here.

These four files in [matlab/](matlab/) are modified copies of files from that repository, and keep
its structure and function names:

- `shift_batch_on_disk2_dredge.m`
- `applyDredge.m`
- `RealignRecordingUsingDredgeAndKilosort2point5.m`
- `Kilosort25Configuration_Dredge.m`

The modifications are listed in [matlab/README.md](matlab/README.md): separating the file's channel
count from the interpolated channel count (and passing the sync channel through untouched),
streaming I/O instead of whole-file copies, removing a DC-offset mixing artifact, per-depth
(nonrigid) shifts with edge clamping, caching the kriging operator, and taking the displacement
reference from the whole recording. `run_interp.m` is new.

That repository in turn adapts **Kilosort 2.5** (`shift_batch_on_disk2.m`, `datashift2.m`) —
Marius Pachitariu and the MouseLand group, <https://github.com/MouseLand/Kilosort> — and the
original per-file credit lines ("Adapted from Kilosort2.5 ... by Richard Hardstone") are preserved.

## Motion estimation

**DREDge** — Charlie Windolf, Erdem Varol and colleagues,
<https://github.com/evarol/dredge> — the decentralized registration method used to estimate motion
from the LFP. It is called here through **SpikeInterface**
(<https://github.com/SpikeInterface/spikeinterface>, Buccino et al.), whose
`estimate_motion(method="dredge_lfp")` and Open Ephys readers this pipeline depends on.

Please cite DREDge, Kilosort and SpikeInterface in any publication that uses this pipeline, and
acknowledge InterpolationAfterDREDge for the interpolation code.

## What is new here

- Python tooling: `dredge_recording.py`, `derive_motion.py`, `make_segment.py`, `check_setup.py`,
  `check_segment.py`.
- Verification against spikes: `template_tracking.py`, `cardiac_calib.py`, `bracket_gain.py` and
  the analysis scripts, plus the ground-truth tests (`injection_test.py`, `synthetic_kriging.py`,
  `lfp_scale_test.py`).
- The per-recording gain calibration, and the findings in
  [FINDINGS_msit001_msit004.md](FINDINGS_msit001_msit004.md).

Developed by Eric Cole (erc84@cornell.edu) for the NPIX_MSIT human cortical Neuropixels recordings. This is an internal lab copy rather than a fork of InterpolationAfterDREDge; upstream
fixes are not automatically picked up, so check that repository for updates.
