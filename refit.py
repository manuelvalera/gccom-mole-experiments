#!/usr/bin/env python3
"""refit.py -- re-measure the beam angle from saved fields at many fit windows,
and ask whether the extrapolated limit depends on the window.

  python refit.py refit-npz

WHY THIS EXISTS

The seiche settled the discretization: at a floor-level alpha it converges at
second order to zero (p = 2.27, e0 = +5.9e-05 on the p=2 fit). The solver core
is fine.

The beam angle does not. Its series converges to a nonzero limit, and the
alpha sweep showed that is NOT the lateral boundary either -- 1e-4 and 1e-6
differ by 0.01-0.10 deg at every grid, while the series shape is unchanged:

    w/N=0.6   +1.86 +0.70 +0.09 -0.22   (alpha 1e-4)
              +1.95 +0.60 +0.02 -0.29   (alpha 1e-6)

So the offset lives in something the flat box cannot see: the ridge near
field, the sponge, or the fit window. This script tests the window, which is
the cheapest of the three because it is pure post-processing -- the same saved
field can be re-measured at any number of windows for no extra simulation.

THE TEST

If the extrapolated limit e0 is the SAME for every window, the beam really is
straight at some angle that is not the theory angle, and the offset is
physical -- finite ridge width, or the near field, or the sponge.

If e0 VARIES with the window, the beam is not straight over these distances,
and the "limit" is an artifact of where we chose to measure. That makes it a
statement about the measurement convention rather than about the model, and
the paper should quote the convention explicitly instead of a single number.

Windows are given as fractions of a bounce, so they mean the same thing at
every frequency. One bounce is D0 / tan(theta).
"""
import os, sys, glob
import numpy as np


def track(rms, XM, ZM, w0, w1):
    """Column-wise maximum of the rms field, then a straight-line fit over the
    window. Same convention as the solver's own tracker."""
    n, m = ZM.shape[0] - 1, ZM.shape[1] - 1
    Xc = np.zeros((n + 2, m + 2)); Zc = np.zeros((n + 2, m + 2))
    Xc[1:-1, 1:-1] = 0.25*(XM[:-1, :-1] + XM[1:, :-1] + XM[:-1, 1:] + XM[1:, 1:])
    Zc[1:-1, 1:-1] = 0.25*(ZM[:-1, :-1] + ZM[1:, :-1] + ZM[:-1, 1:] + ZM[1:, 1:])
    for A in (Xc, Zc):
        A[0, :] = A[1, :]; A[-1, :] = A[-2, :]
        A[:, 0] = A[:, 1]; A[:, -1] = A[:, -2]
    pts = []
    for j in range(1, m + 1):
        col = rms[1:-1, j]
        if not np.isfinite(col).any() or col.max() <= 0:
            continue
        k = int(np.nanargmax(col))
        xv = abs(Xc[k + 1, j])
        if w0 <= xv <= w1:
            pts.append((xv, Zc[k + 1, j]))
    if len(pts) < 8:
        return None, None, 0
    P = np.array(pts)
    A = np.vstack([P[:, 0], np.ones(len(P))]).T
    sl, ic = np.linalg.lstsq(A, P[:, 1], rcond=None)[0]
    res = P[:, 1] - A @ [sl, ic]
    return (np.degrees(np.arctan(abs(sl))),
            float(np.sqrt((res**2).mean())), len(P))


def extrapolate(h, e):
    """Fit e = A h^p + e0 for p = 1 and 2; return both e0 and the better fit."""
    h, e = np.asarray(h, float), np.asarray(e, float)
    out = {}
    for p in (1.0, 2.0):
        M = np.vstack([h**p, np.ones_like(h)]).T
        (A, e0), *_ = np.linalg.lstsq(M, e, rcond=None)
        out[p] = (e0, float(np.abs(e - (A*h**p + e0)).max()))
    return out


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'refit-npz'
    files = sorted(glob.glob(os.path.join(d, '*.npz')))
    if not files:
        raise SystemExit(f"no .npz in {d}/ -- run refit_fields.ps1 first")

    runs = []
    for f in files:
        z = np.load(f)
        ratio = float(z['ratio'])
        nz, nx = z['Z'].shape
        runs.append(dict(f=os.path.basename(f), rms=z['rms'], X=z['X'], Z=z['Z'],
                         ratio=ratio, nx=nx, nz=nz, dz=1000.0/(nz - 1)))
    print(f"\n{len(runs)} saved fields in {d}/\n")

    # windows as fractions of one bounce, so they are comparable across w/N
    fracs = [(0.10, 0.40), (0.10, 0.55), (0.10, 0.75),
             (0.10, 0.95), (0.30, 0.75), (0.50, 0.95)]

    for ratio in sorted({r['ratio'] for r in runs}):
        rs = sorted([r for r in runs if r['ratio'] == ratio],
                    key=lambda z: -z['dz'])
        if len(rs) < 3:
            print(f"w/N={ratio}: only {len(rs)} grids, need 3+ to extrapolate")
            continue
        theory = np.degrees(np.arcsin(ratio))
        bounce = 1000.0/np.tan(np.radians(theory))
        print("=" * 78)
        print(f"w/N = {ratio}   theory {theory:.2f} deg   one bounce = {bounce:.0f} m")
        print("=" * 78)
        print("%-14s %-28s %-11s %-11s" %
              ("window (bnc)", "bias by grid (coarse->fine)", "e0 (p=1)", "e0 (p=2)"))
        e0s = []
        for f0, f1 in fracs:
            w0, w1 = f0*bounce, f1*bounce
            row, hs, es, ok = [], [], [], True
            for r in rs:
                a, rms, npt = track(r['rms'], r['X'], r['Z'], w0, w1)
                if a is None or rms > 10.0:
                    row.append("  n/a "); ok = False; continue
                row.append("%+6.2f" % (a - theory))
                hs.append(r['dz']); es.append(a - theory)
            if not ok or len(hs) < 3:
                print("%-14s %-28s %-11s %-11s"
                      % (f"{f0:.2f}-{f1:.2f}", " ".join(row), "--", "--"))
                continue
            ex = extrapolate(hs, es)
            e0s.append(ex[2.0][0])
            print("%-14s %-28s %-11s %-11s"
                  % (f"{f0:.2f}-{f1:.2f}", " ".join(row),
                     "%+.2f" % ex[1.0][0], "%+.2f" % ex[2.0][0]))
        if len(e0s) >= 3:
            spread = max(e0s) - min(e0s)
            print(f"\n   e0 across windows: {min(e0s):+.2f} to {max(e0s):+.2f}"
                  f"   spread {spread:.2f} deg")
            if spread > 0.3:
                print("   -> the limit DEPENDS ON THE WINDOW. The beam is not"
                      " straight at a single")
                print("      angle over these distances, so there is no one"
                      " 'converged bias' to")
                print("      quote -- the measurement convention has to be"
                      " stated with the number.")
            else:
                print("   -> the limit is window-INDEPENDENT. The beam is"
                      " straight at an angle that")
                print("      is genuinely offset from theory; look to the ridge"
                      " near field or the")
                print("      sponge, not to the measurement.")
        print()


if __name__ == '__main__':
    main()
