#!/usr/bin/env python3
"""wiggle.py -- what is the oscillation on the tracked beam maxima?

  python wiggle.py a.npz [b.npz ...]

Section 9 showed the beam is straight (a line fits the maxima to ~3 m over
~900 m) but that short fit windows give slopes varying by 2-3 deg, and that
the variation GROWS under refinement.  That is the signature of an oscillation
whose amplitude does not fall with dx.

This refits the stored maxima, then measures the residual about that line:
its amplitude, and its dominant wavelength expressed BOTH in metres and in
grid cells.  The distinction is the whole point:

  * wavelength fixed in CELLS (~2 dx) and amplitude ~constant in cells
      -> a grid mode.  Refining makes it finer, not weaker, so the slope
         error over a fixed window does not improve.  Numerical.
  * wavelength fixed in METRES across grids
      -> real structure in the solution (beam sub-structure, interference
         with a reflected or secondary beam).  Physical.

Run it on the same frequency at two or three resolutions and compare.
"""
import os, sys
import numpy as np


def analyse(fn):
    d = np.load(fn)
    pts = d['pts'] if 'pts' in d else None
    if pts is None or len(pts) < 8:
        print(f"{os.path.basename(fn)}: no stored maxima"); return None
    X, Z = d['X'], d['Z']
    nx, nz = X.shape[1] - 1, X.shape[0] - 1
    Lx = float(X.max() - X.min()); D0 = float(-Z.min())
    dx, dz = Lx/nx, D0/nz

    x, z = pts[:, 0], pts[:, 1]
    o = np.argsort(x); x, z = x[o], z[o]
    A = np.vstack([x, np.ones(len(x))]).T
    sl, ic = np.linalg.lstsq(A, z, rcond=None)[0]
    res = z - (sl*x + ic)

    # dominant wavelength of the residual, via a periodogram on a uniform
    # resample (the maxima are one per column, so x is already near-uniform)
    xs = np.linspace(x.min(), x.max(), len(x))
    rs = np.interp(xs, x, res)
    rs = rs - rs.mean()
    F = np.abs(np.fft.rfft(rs * np.hanning(len(rs))))
    k = np.fft.rfftfreq(len(rs), d=(xs[1]-xs[0]))
    F[0] = 0.0
    lam = 1.0/k[np.argmax(F)] if k[np.argmax(F)] > 0 else np.inf
    # A peak at (or near) the fundamental means the residual is one gentle
    # arch across the window, not an oscillation, and the "wavelength" is
    # meaningless -- it just reports the span back.  Say so instead of
    # dressing it up as a measurement.
    degenerate = lam > 0.6*(xs[-1] - xs[0])

    span = x.max() - x.min()
    print(f"{os.path.basename(fn):22s} w/N={float(d['ratio']):.1f}  "
          f"{nx}x{nz}  dx={dx:.1f} dz={dz:.1f} m")
    print(f"   fit {np.degrees(np.arctan(abs(sl))):.2f} deg   "
          f"span {span:.0f} m, {len(x)} pts")
    print(f"   residual: rms {res.std():.2f} m  peak {np.abs(res).max():.2f} m"
          f"   = {res.std()/dz:.2f} dz")
    if degenerate:
        print(f"   residual is ONE ARCH across the window, not an oscillation "
              f"(peak at {lam:.0f} m vs span {xs[-1]-xs[0]:.0f} m)")
        print(f"   -> no wavelength to report; the metres-vs-cells test below "
              f"does not apply")
    else:
        print(f"   dominant wavelength {lam:.0f} m = {lam/dx:.1f} dx")
    # slope error a window of length L would inherit from this oscillation
    for L in (200, 300, 500):
        if L < span:
            err = np.degrees(np.arctan(2*np.sqrt(2)*res.std()/L))
            print(f"   -> a {L} m window would see about {err:.2f} deg of slope error")
    return dict(dx=dx, dz=dz, lam=lam, rms=res.std(), degen=degenerate)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    out = [analyse(f) for f in sys.argv[1:]]
    out = [o for o in out if o]
    if len(out) > 1:
        print("\nacross grids:")
        if any(o['degen'] for o in out):
            print("   (at least one residual is a single arch -- the")
            print("    wavelength comparison below is not meaningful; read the")
            print("    amplitude column instead)")
        print("   %-10s %-12s %-12s %-10s" % ("dx (m)", "lambda (m)", "lambda/dx", "rms/dz"))
        for o in out:
            print("   %-10.1f %-12.0f %-12.1f %-10.2f"
                  % (o['dx'], o['lam'], o['lam']/o['dx'], o['rms']/o['dz']))
        # amplitude scaling is the reliable diagnostic
        r = np.array([o['rms'] for o in out]); h = np.array([o['dx'] for o in out])
        if len(r) > 1 and r[0] > 0 and h[0] != h[-1]:
            p = np.log(r[-1]/r[0])/np.log(h[-1]/h[0])
            print(f"\n   residual amplitude scales as dx^{p:.2f}"
                  + ("  -> shrinking with refinement: ordinary discretization"
                     if p > 0.5 else
                     "  -> NOT shrinking: something grid-locked"))
        if any(o['degen'] for o in out):
            return
        lam_m = np.array([o['lam'] for o in out])
        lam_c = np.array([o['lam']/o['dx'] for o in out])
        vm = lam_m.std()/lam_m.mean()
        vc = lam_c.std()/lam_c.mean()
        print(f"\n   scatter of lambda in metres {vm:.2f}, in cells {vc:.2f}")
        if vc < vm:
            print("   -> wavelength is fixed in CELLS: a grid mode, not physical.")
            print("      Refining will not reduce the slope error over a fixed")
            print("      window; only lengthening the window will.")
        else:
            print("   -> wavelength is fixed in METRES: real structure in the")
            print("      solution, and refining should eventually resolve it.")


if __name__ == '__main__':
    main()
