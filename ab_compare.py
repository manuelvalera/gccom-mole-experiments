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
            for g in grids:
                C, th, b, D0 = centroid_points(data[ab][g])
                a = fit_angle(C, lo*b, hi*b)
                hs.append(D0/(g[1]-1)); es.append(np.nan if a is None else a - th)
            hs, es = np.array(hs), np.array(es)
            good = np.isfinite(es)
            if good.sum() >= 2:
                M = np.vstack([hs[good], np.ones(good.sum())]).T
                (A, e0), *_ = np.linalg.lstsq(M, es[good], rcond=None)
            else:
                e0 = np.nan
            lims[ab] = e0
            print("%-8g" % ab + "".join("%+-12.2f" % e for e in es) + "%+.2f" % e0)
        vals = np.array([v for v in lims.values() if np.isfinite(v)])
        if len(vals) >= 2:
            print("   spread of the limit across ridge heights: %.2f deg" % (vals.max() - vals.min()))

    print("\nA linear beam's angle cannot depend on ridge height. A spread well")
    print("above ~0.2 deg means the barotropic tide is contaminating the rms")
    print("diagnostic; a spread near zero means the offset is in the beam itself.")


if __name__ == '__main__':
    main()
