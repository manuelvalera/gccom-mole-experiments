#!/usr/bin/env python3
"""time_compare.py -- beam angle and field statistics against run length.

  python time_compare.py [DIR]      (reads DIR/ab*_p*.npz, default time-npz)

Each run's rms is averaged over its last 10 periods, so nper = 20, 40, 60, 80
sample the solution at periods 10-20, 30-40, 50-60, 70-80. A converged,
stable solution gives the same angle in every row.
"""
import os, re, sys, glob
import numpy as np
from scipy.interpolate import LinearNDInterpolator
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from centroid_track import centres, centroid_track, fit_angle, argmax_track

rows = {}
DIR = sys.argv[1] if len(sys.argv) > 1 else 'time-npz'
for f in sorted(glob.glob(os.path.join(DIR, 'ab*_p*.npz'))):
    m = re.search(r'ab(\d+)_p(\d+)\.npz$', f)
    if not m:
        continue
    ab, npr = int(m.group(1)), int(m.group(2))
    z = np.load(f)
    theta = np.arcsin(float(z['ratio'])); th = np.degrees(theta)
    D0 = float(np.abs(z['Z']).max()); b = D0/np.tan(theta)
    x, zz, r = centres(z['X'], z['Z'], z['rms'])
    ok = np.isfinite(r)
    interp = LinearNDInterpolator(np.c_[x[ok], zz[ok]], r[ok])
    s = np.linspace(0.05, 0.90, 170)*b/np.cos(theta)
    C = np.vstack([centroid_track(interp, D0, theta, s, 150.0, 0.5, sd) for sd in (+1, -1)])
    A = argmax_track(z['X'], z['Z'], z['rms'])
    rr = r[ok]; mm = rr.max()/rr.mean()
    out = []
    for P in (C, A):
        for lo, hi in ((0.10, 0.40), (0.10, 0.75)):
            a = fit_angle(P, lo*b, hi*b)
            out.append(np.nan if a is None else a - th)
    nzz, nxx = z['Z'].shape
    rows.setdefault((ab, f'{nxx-1}x{nzz-1}'), []).append((npr, mm, float(rr.max()), out))

for (ab, grid), rs in sorted(rows.items()):
    print(f"\nab = {ab} m, {grid}   (bias vs theory, deg; rms averaged over the last 10 periods)")
    print("%-7s %-11s %-10s %-12s %-12s %-12s %-12s" % ("nper", "averages", "max/mean",
          "cent .10-.40", "cent .10-.75", "argm .10-.40", "argm .10-.75"))
    for npr, mm, mx, o in sorted(rs):
        print("%-7d %-11s %-10.1f %s" % (npr, f"{npr-10}-{npr}", mm,
              " ".join("%+-12.2f" % v for v in o)))
print("\nConverged and stable: every column flat down the table, max/mean flat.")
print("Spin-up: columns change early and then settle.")
print("Growth: max/mean keeps rising and the angles keep drifting.")
