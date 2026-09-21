#!/usr/bin/env python3
"""curvature.py -- is the beam straight?

  python curvature.py refit-npz

WHERE THIS COMES FROM

refit.py showed the extrapolated bias depends on the fit window: spread 1.22
deg at w/N=0.6 and 2.23 deg at 0.8, running from about zero at the innermost
window to -1 or -2 deg at the outermost. A straight beam cannot do that. The
natural reading is that the beam is CURVED, so "the beam angle" is only
defined once you say over what range it was measured.

This script tests that directly instead of inferring it, three ways.

1. LOCAL ANGLE vs DISTANCE. Track the maxima, then fit a straight line in
   successive bins of |x|. If the beam is straight the local angle is flat
   across bins; if it curves, the angle drifts monotonically. This is the
   plot the paper needs, because it shows the reader what the window choice
   is actually choosing between.

2. LINE vs QUADRATIC. Fit z = a + b|x| and z = a + b|x| + c|x|^2 over the
   same range and compare residuals. A curvature term that stays put as the
   grid refines is physical; one that shrinks is numerical.

3. THE PAPER'S OWN WINDOW. Garcia et al. fit over a fixed [200,500] m. At
   w/N=0.6 one bounce is 1333 m, so that is 0.15-0.375 of a bounce -- close
   to the innermost window here, where the extrapolated bias is near zero.
   At 0.8 a bounce is 750 m, so the same 200-500 m is 0.27-0.67, much further
   out. Reporting our numbers in the paper's own convention is the honest
   comparison against Fig 9, whatever we conclude about curvature.
"""
import os, sys, glob
import numpy as np


def maxima(rms, XM, ZM):
    """Column-wise maxima of the rms field, as (|x|, z) pairs."""
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
        pts.append((abs(Xc[k + 1, j]), Zc[k + 1, j]))
    P = np.array(sorted(pts))
    return P[:, 0], P[:, 1]


def angle_of(x, z):
    if len(x) < 6:
        return None, None
    A = np.vstack([x, np.ones(len(x))]).T
    sl, ic = np.linalg.lstsq(A, z, rcond=None)[0]
    res = z - A @ [sl, ic]
    return np.degrees(np.arctan(abs(sl))), float(np.sqrt((res**2).mean()))


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'refit-npz'
    files = sorted(glob.glob(os.path.join(d, '*.npz')))
    if not files:
        raise SystemExit(f"no .npz in {d}/")
    runs = []
    for f in files:
        z = np.load(f)
        nz, nx = z['Z'].shape
        # D0 from the grid itself, not assumed: the depth test varies it,
        # and a hardcoded 1000 would silently rescale every bounce fraction.
        D0 = float(np.abs(z['Z']).max())
        runs.append(dict(f=os.path.basename(f), ratio=float(z['ratio']),
                         rms=z['rms'], X=z['X'], Z=z['Z'], D0=D0,
                         nx=nx, nz=nz, dz=D0/(nz - 1)))

    for key in sorted({(r['ratio'], round(r['D0'])) for r in runs}):
        ratio, D0r = key
        rs = sorted([r for r in runs
                     if r['ratio'] == ratio and round(r['D0']) == D0r],
                    key=lambda z: -z['dz'])
        theory = np.degrees(np.arcsin(ratio))
        bounce = D0r/np.tan(np.radians(theory))
        print("=" * 76)
        print(f"w/N = {ratio}   D0 = {D0r:.0f} m   theory {theory:.2f} deg"
              f"   bounce {bounce:.0f} m")
        print("=" * 76)

        # ---- 1. local angle in bins of |x| --------------------------------
        edges = np.array([0.10, 0.25, 0.40, 0.55, 0.70, 0.85])*bounce
        print("\n1. LOCAL angle by distance band (bias vs theory, deg)")
        hdr = "  ".join("%4.2f-%4.2f" % (edges[i]/bounce, edges[i+1]/bounce)
                        for i in range(len(edges)-1))
        print("%-12s %s" % ("grid", hdr))
        for r in rs:
            x, z = maxima(r['rms'], r['X'], r['Z'])
            row = []
            for i in range(len(edges)-1):
                sel = (x >= edges[i]) & (x <= edges[i+1])
                a, _ = angle_of(x[sel], z[sel])
                row.append("   n/a   " if a is None else "%+8.2f " % (a - theory))
            print("%-12s %s" % (f"{r['nx']}x{r['nz']}", " ".join(row)))
        print("   flat across a row = straight beam; drifting = curved")

        # ---- 2. line vs quadratic ----------------------------------------
        print("\n2. LINE vs QUADRATIC over 0.10-0.75 of a bounce")
        print("%-12s %-11s %-11s %-13s %-11s" %
              ("grid", "lin rms", "quad rms", "curv c (1/m)", "c*L^2 (m)"))
        L = 0.65*bounce
        for r in rs:
            x, z = maxima(r['rms'], r['X'], r['Z'])
            sel = (x >= 0.10*bounce) & (x <= 0.75*bounce)
            xs, zs = x[sel], z[sel]
            if len(xs) < 10:
                print("%-12s insufficient points" % f"{r['nx']}x{r['nz']}"); continue
            A1 = np.vstack([xs, np.ones(len(xs))]).T
            c1 = np.linalg.lstsq(A1, zs, rcond=None)[0]
            r1 = np.sqrt(((zs - A1 @ c1)**2).mean())
            A2 = np.vstack([xs**2, xs, np.ones(len(xs))]).T
            c2 = np.linalg.lstsq(A2, zs, rcond=None)[0]
            r2 = np.sqrt(((zs - A2 @ c2)**2).mean())
            print("%-12s %-11.2f %-11.2f %-13.3e %-11.1f" %
                  (f"{r['nx']}x{r['nz']}", r1, r2, c2[0], c2[0]*L**2))
        print("   c*L^2 is the sag across the window, in metres.")
        print("   If it holds steady as the grid refines the curvature is real;")
        print("   if it shrinks it was numerical.")

        # ---- 3. the paper's own fixed window ------------------------------
        print("\n3. GARCIA ET AL.'S OWN WINDOW, fixed [200,500] m")
        print("   (= %.2f-%.2f of a bounce at this frequency)"
              % (200.0/bounce, 500.0/bounce))
        print("%-12s %-10s %-10s %-8s" % ("grid", "measured", "bias", "pts"))
        hs, es = [], []
        for r in rs:
            x, z = maxima(r['rms'], r['X'], r['Z'])
            sel = (x >= 200.0) & (x <= 500.0)
            a, rr = angle_of(x[sel], z[sel])
            if a is None:
                print("%-12s insufficient points" % f"{r['nx']}x{r['nz']}"); continue
            print("%-12s %-10.2f %+-10.2f %-8d"
                  % (f"{r['nx']}x{r['nz']}", a, a - theory, sel.sum()))
            hs.append(r['dz']); es.append(a - theory)
        if len(hs) >= 3:
            for p in (1.0, 2.0):
                M = np.vstack([np.array(hs)**p, np.ones(len(hs))]).T
                (A, e0), *_ = np.linalg.lstsq(M, np.array(es), rcond=None)
                resid = np.abs(np.array(es) - (A*np.array(hs)**p + e0)).max()
                print("   extrapolated e0 (p=%.0f) = %+.2f deg   resid %.2f"
                      % (p, e0, resid))
        print()


if __name__ == '__main__':
    main()
