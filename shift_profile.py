#!/usr/bin/env python3
"""shift_profile.py -- move the pycnocline up or down, keeping its shape.

  python shift_profile.py mry_N.txt --dz 10 -o mry_N_dn10.txt    (deeper)
  python shift_profile.py mry_N.txt --dz -10 -o mry_N_up10.txt   (shallower)

Whether a shoaling internal wave arrives as a depression or an elevation, and
so what shape its bores take, depends on where the pycnocline sits relative to
mid-depth. Walter et al. (2012) initialised with CTD casts from near the time of
their observations; the profile here is a fit of unknown date. Shifting it
tests whether the canonical/non-canonical shape is sensitive to that at all.

The table is (z, N) with z negative downward. A positive --dz moves every
feature DOWN by that many metres; values that would leave the column are
dropped and the vacated end is filled with the nearest remaining N.
"""
import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('profile')
    ap.add_argument('--dz', type=float, required=True,
                    help='metres to move the profile down (negative: up)')
    ap.add_argument('-o', '--out', required=True)
    g = ap.parse_args()

    p = np.loadtxt(g.profile)
    order = np.argsort(p[:, 0])
    z = p[order, 0]
    N = p[order, 1]
    ztop = float(z.max())
    zbot = float(z.min())
    shifted = np.interp(z + g.dz, z, N, left=N[0], right=N[-1])
    np.savetxt(g.out, np.column_stack([z, shifted]), fmt='%.4f %.6e')
    k0 = int(np.argmax(N))
    k1 = int(np.argmax(shifted))
    print(f"{g.out}: N max moved from z = {z[k0]:.1f} m to z = {z[k1]:.1f} m "
          f"(column {zbot:.1f} to {ztop:.1f} m)")


if __name__ == '__main__':
    main()
