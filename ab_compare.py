#!/usr/bin/env python3
"""ab_compare.py -- beam-angle bias vs ridge height, on the same grids.

  python ab_compare.py 10=ab10-npz 20=refit-npz 40=ab40-npz

Uses the same perpendicular-slice energy-centroid tracker as centroid_track.py
(imported from it, so the two cannot drift apart). Only w/N = 0.8 fields on
grids present for EVERY ridge height are used, so each row is like for like.
"""
import os, sys, glob
import numpy as np
from scipy.interpolate import LinearNDInterpolator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from centroid_track import centres, centroid_track, fit_angle


def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, '*.npz'))):
        z = np.load(f)
        if abs(float(z['ratio']) - 0.8) > 1e-9:
            continue
        nz, nx = z['Z'].shape
        out[(nx, nz)] = z
    return out


def quality(z, C, lo, hi, b):
    """Flag fields the tracker cannot be trusted on: a localised spike in the
    rms field (max/mean far above the ~10-20 of a clean beam) or a centroid
    track that is not a line (linear-fit rms above 10 m)."""
    r = z['rms'][1:-1, 1:-1]
    r = r[np.isfinite(r)]
    mm = float(r.max()/r.mean()) if r.size and r.mean() > 0 else np.inf
    sel = (C[:, 0] >= lo*b) & (C[:, 0] <= hi*b)
    if sel.sum() < 6:
        return False, mm, np.inf
    x, zz = C[sel, 0], C[sel, 1]
    res = zz - np.polyval(np.polyfit(x, zz, 1), x)
    lr = float(np.sqrt((res**2).mean()))
    return (mm < 40.0 and lr < 10.0), mm, lr


def centroid_points(z):
    ratio = float(z['ratio']); D0 = float(np.abs(z['Z']).max())
    theta = np.arcsin(ratio); bounce = D0/np.tan(theta)
    x, zz, r = centres(z['X'], z['Z'], z['rms'])
    ok = np.isfinite(r)
    interp = LinearNDInterpolator(np.c_[x[ok], zz[ok]], r[ok])
    s_list = np.linspace(0.05, 0.90, 170)*bounce/np.cos(theta)
    C = np.vstack([centroid_track(interp, D0, theta, s_list, 150.0, 0.5, sd)
                   for sd in (+1, -1)])
    return C, np.degrees(theta), bounce, D0


def main():
    specs = [a.split('=', 1) for a in sys.argv[1:]]
    if not specs:
        raise SystemExit("usage: python ab_compare.py 10=ab10-npz 20=refit-npz 40=ab40-npz")
    data = {float(ab): load(d) for ab, d in specs}
    common = set.intersection(*[set(v) for v in data.values()])
    if not common:
        raise SystemExit("no grid is present for every ridge height")
    grids = sorted(common)
    print(f"\nw/N = 0.8, centroid tracker, grids common to all ridge heights: "
          + ", ".join(f"{nx-1}x{nz-1}" for nx, nz in grids))

    for lo, hi in ((0.10, 0.40), (0.10, 0.75)):
        print(f"\nBIAS over {lo:.2f}-{hi:.2f} of a bounce (deg)")
        print("%-8s" % "ab (m)" + "".join("%-12s" % f"{nx-1}x{nz-1}" for nx, nz in grids)
              + "limit p=1")
        lims = {}
        for ab in sorted(data):
            hs, es = [], []
            flags = []
            for g in grids:
                C, th, b, D0 = centroid_points(data[ab][g])
                ok, mm, lr = quality(data[ab][g], C, lo, hi, b)
                a = fit_angle(C, lo*b, hi*b)
                hs.append(D0/(g[1]-1))
                es.append(a - th if (a is not None and ok) else np.nan)
                if not ok:
                    flags.append(f"{g[0]-1}x{g[1]-1} (max/mean {mm:.0f}, fit rms {lr:.1f} m)")
            hs, es = np.array(hs), np.array(es)
            good = np.isfinite(es)
            if good.sum() >= 2:
                M = np.vstack([hs[good], np.ones(good.sum())]).T
                (A, e0), *_ = np.linalg.lstsq(M, es[good], rcond=None)
            else:
                e0 = np.nan
            lims[ab] = e0
            print("%-8g" % ab + "".join(("%+-12.2f" % e) if np.isfinite(e) else "REJECTED    "
                                        for e in es) + ("%+.2f" % e0 if np.isfinite(e0) else "n/a"))
            for f in flags:
                print("        rejected " + f)
        vals = np.array([v for v in lims.values() if np.isfinite(v)])
        if len(vals) >= 2:
            print("   spread of the limit across ridge heights: %.2f deg" % (vals.max() - vals.min()))

    print("\nA linear beam's angle cannot depend on ridge height, so any spread is")
    print("signal. Read its DIRECTION as well as its size:")
    print("   offset LARGER for a SMALLER ridge -> the barotropic tide is")
    print("      contaminating the rms diagnostic (the beam is weaker relative")
    print("      to the tide), and the fix is to track baroclinic rms;")
    print("   offset LARGER for a LARGER ridge  -> the offset scales with the")
    print("      topography itself -- grid distortion near the ridge, or the")
    print("      curved-bed treatment -- and should shrink as the ridge does.")
    print("Rejected fields are excluded from the limits; inspect them with view.py.")


if __name__ == '__main__':
    main()
