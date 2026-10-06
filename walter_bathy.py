#!/usr/bin/env python3
"""walter_bathy.py -- the idealised bathymetries of Walter et al. (2012).

  python walter_bathy.py --case noncanonical -o walter_xi2.csv
  python walter_bathy.py --case canonical -o walter_xi02.csv
  python mry_setup.py --bathy walter_xi2.csv --out grids\\walter_xi2 --smooth 0

Their SUNTANS runs did not use the measured bathymetry. model_setup.m in the
archived setup (rectddatafromsuntansmodel.zip) builds two analytic profiles
over a 20 km transect, x = 0 at the shallow end:

  xi ~ 2 (non-canonical):  z = ( -35 - 47 tanh((x-700)/300)
                                 -40 - 40 tanh((x-2000)/4000) ) / 2
  xi ~ 0.2 (canonical):    z =   -35 - 47 tanh((x-4500)/6000)

The first reproduces the wall-like step of the real transect along the bore
path (15 to 50 m within ~500 m, slope up to 0.083 at 31 m depth); the second is
ten times gentler. Their stratification is a two-tanh fit to the 1 April 2010
MBARI C1 cast, which mry_N.txt already matches to within 0.1 degC.

Written as x, y, z triples on a regular grid three columns wide, the layout
mry_setup.py's read_transect expects, with y the along-transect distance. Pass
--smooth 0 to mry_setup.py: its default 200 m smoothing would flatten a step
whose own width is 300 m.
"""
import argparse

import numpy as np


def profile(case, x):
    if case == 'noncanonical':
        z1 = -35.0 - 47.0 * np.tanh((x - 700.0) / 300.0)
        z2 = -40.0 - 40.0 * np.tanh((x - 2000.0) / 4000.0)
        return 0.5 * (z1 + z2)
    if case == 'canonical':
        return -35.0 - 47.0 * np.tanh((x - 4500.0) / 6000.0)
    raise SystemExit(f'unknown case {case}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=['noncanonical', 'canonical'], required=True)
    ap.add_argument('--length', type=float, default=20000.0)
    ap.add_argument('--dx', type=float, default=10.0, help='along-transect spacing, m')
    ap.add_argument('-o', '--out', required=True)
    g = ap.parse_args()

    y = np.arange(0.0, g.length + 0.5 * g.dx, g.dx)
    z = profile(g.case, y)
    cols = np.array([0.0, 100.0, 200.0])
    rows = []
    for yy, zz in zip(y, z):
        for xx in cols:
            rows.append((xx, yy, zz))
    np.savetxt(g.out, np.array(rows), fmt='%.2f %.2f %.4f')

    d = -z
    s = np.abs(np.gradient(z, y))
    k15 = int(np.argmax(d >= 15.0))
    print(f"{g.out}: {g.case}, {len(y)} points over {g.length/1000:.1f} km, "
          f"depth {d.min():.1f} to {d.max():.1f} m")
    print(f"   steepest slope {s.max():.3f} at {d[np.argmax(s)]:.0f} m depth; "
          f"15 m isobath {y[k15]:.0f} m from the shallow end")


if __name__ == '__main__':
    main()
