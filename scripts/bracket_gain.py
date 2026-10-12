"""Find the slow-drift gain at which the correction neither over- nor under-shoots.

    python bracket_gain.py <segA> <gainA> <segB> <gainB>

e.g.  python bracket_gain.py seg_g100 1.0 seg_g50 0.5

Why a gain is needed at all
---------------------------
The interpolation applies ~100% of whatever displacement it is given (verified by registration and
by injecting known motion). Whether that REMOVES motion depends on whether DREDge's displacement
matches how the units actually move. On msit004 units moved only ~0.6x DREDge-LFP's slow estimate,
so correcting at face value left them moving backwards by almost as much as before.

Run template_tracking.py on two interpolated segments made with different slow gains, then pass the
two segment names here. For each unit, template_tracking regresses its tracked position on the
displacement that was APPLIED to that file; scaled back to units of the ORIGINAL (gain 1) estimate,
the after-correction slope is

    s_after(g) = tracker_gain * (s_unit - g)        s_unit = how much of DREDge the units really move

so the median after-slope is linear in g and crosses zero at g = s_unit. Two runs determine that
line; the printed crossing is the gain to use. Judge a run by the SIGN of the after-correction
slope/r, which does not depend on the tracker's gain:

    negative -> over-corrected (gain too high)    positive -> under-corrected    ~0 -> right

Then set "slow_gain" in recordings.json, derive the final motion with that gain, and re-run.
"""
import sys
from pathlib import Path

import numpy as np

import recinfo

OUT = recinfo.RESULTS


def median_slope(seg, gain):
    """Median after-correction slope and r for one run, in units of the gain-1 estimate."""
    p = OUT / ("%s%s_tracking.npz" % (recinfo.key() + "_", seg))
    if not p.exists():
        raise SystemExit("no tracking result at %s - run template_tracking.py %s <t0> <t1> first" % (p, seg))
    d = np.load(p)
    sa, r = d["slope_after"], d["r_after"]
    m = np.isfinite(sa) & np.isfinite(r)
    if m.sum() < 3:
        raise SystemExit("%s: only %d units tracked - widen the time range or loosen the gates" % (seg, m.sum()))
    print("%-16s gain %.2f: n=%2d  median slope after %+.2f (x applied) = %+.2f (x gain-1 estimate), median r %+.2f"
          % (seg, gain, m.sum(), np.median(sa[m]), gain * np.median(sa[m]), np.median(r[m])))
    return gain * float(np.median(sa[m])), float(np.median(r[m]))


def main(segA, gA, segB, gB):
    yA, rA = median_slope(segA, gA)
    yB, rB = median_slope(segB, gB)
    if abs(gA - gB) < 1e-9:
        raise SystemExit("the two runs must use different gains")
    if (yA > 0) == (yB > 0):
        print("\nBoth runs land on the same side (%s-corrected). The crossing is an extrapolation, "
              "so run a third gain beyond them before trusting it." % ("under" if yA > 0 else "over"))
    slope = (yB - yA) / (gB - gA)
    g = gA - yA / slope if abs(slope) > 1e-9 else np.nan
    print("\nbest slow gain ~ %.2f   (units move %.0f%% of DREDge's slow displacement)" % (g, 100 * g))
    print("Set \"slow_gain\": %.2f in recordings.json, then:" % round(g, 2))
    print("  python derive_motion.py <baseline_tag> <new_tag> \"0,0.1,%.2f\"" % round(g, 2))
    print("and re-run make_segment.py + run_interp.m + template_tracking.py to confirm slope/r ~ 0.")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3], float(sys.argv[4]))
