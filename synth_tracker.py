"""Feed the column-max + lstsq tracker a beam of KNOWN angle and see what it
reports back.

HANDOFF sec 5 records "measurement bias in the angular peak-find" as a tested
dead end -- but that was cone3d's spherical-shell peak-find, tested against a
synthetic cone.  iwbcurv uses a different measurement entirely: per-column
vertical argmax of the rms field, then a least-squares line through those
(|x|, z) maxima.  That one has never been tested against a known answer.

This builds an analytic beam at a prescribed angle on the REAL curvilinear
grid, with the real Xcen/Zcen, and runs the identical fit block.  Any
difference between input and output is measurement, not physics.

Variants:
  clean   the four beam arms only (the V and its mirror)
  refl    plus one surface and one bed reflection per arm, which is what the
          field actually contains after a few periods
"""
import os, sys
import numpy as np
sys.path.insert(0, '.')
from oct2py import Oct2Py

MOLE = (sys.argv[1] if len(sys.argv) > 1 else os.environ.get('MOLE_SRC', ''))
if not MOLE or not os.path.isfile(os.path.join(MOLE, 'gridGen.m')):
    raise SystemExit(
        "set MOLE_SRC (or pass the path as argv[1]) to MOLE's src/matlab_octave.\n"
        "  Octave's addpath does NOT error on a missing directory, so a wrong\n"
        "  path shows up later as 'gridGen undefined' -- checking it here instead.\n"
        "  got: %r" % MOLE)
Lx, D0 = 6000.0, 1000.0
nx, nz = 256, 101


def build_grid():
    open('grids/iwbridge/geom_over.m', 'w').write(
        "function [Lx, D0, ab, Lb] = geom_over()\n"
        "    Lx = %r; D0 = %r; ab = 20.0; Lb = 30.0;\nend\n" % (Lx, D0))
    oc = Oct2Py(); oc.eval("warning('off','all'); more off;")
    oc.addpath(MOLE); oc.addpath(os.path.abspath('grids/iwbridge'))
    here = os.getcwd(); os.chdir('grids')
    try:
        oc.eval("global BT BB; BT=1.6; BB=3.2;")
        oc.eval("cd('%s'); [X,Z]=gridGen('TFI','iwbridge',%d,%d,false);"
                % (os.path.abspath('.').replace(os.sep, '/'), nx, nz))
    finally:
        os.chdir(here)
    X = oc.pull('X'); Z = oc.pull('Z'); oc.exit()
    return X.T.copy(), Z.T.copy()


XM, ZM = build_grid()
m, n = nx - 1, nz - 1

# cell centres, exactly as iwbcurv builds them
Xcen = np.zeros((n + 2, m + 2)); Zcen = np.zeros((n + 2, m + 2))
Xin = 0.25*(XM[:-1, :-1] + XM[1:, :-1] + XM[:-1, 1:] + XM[1:, 1:])
Zin = 0.25*(ZM[:-1, :-1] + ZM[1:, :-1] + ZM[:-1, 1:] + ZM[1:, 1:])
Xcen[1:-1, 1:-1] = Xin; Zcen[1:-1, 1:-1] = Zin
for Ac in (Xcen, Zcen):
    Ac[0, :] = Ac[1, :]; Ac[-1, :] = Ac[-2, :]
    Ac[:, 0] = Ac[:, 1]; Ac[:, -1] = Ac[:, -2]


def synth(theta_deg, sigma=40.0, nrefl=0, ramp=0.6, decay=3000.0, skew=0.0):
    """Beam leaving the ridge crest, folded off the surface/bed nrefl times.

    NOTE: the first version of this file summed four arms each of peak
    amplitude 1.0, so every column had several tied maxima and argmax picked
    arbitrarily -- it reported -9 to -33 deg of "bias" that was entirely an
    artefact of the synthetic field.  The `seg` mask below keeps exactly one
    arm alive at each x, which is what the physical field looks like inside
    the first bounce (where the fit window sits).
    """
    sl = np.tan(np.radians(theta_deg))
    ax = np.abs(Xcen)
    F = np.zeros_like(Xcen)
    for b in range(nrefl + 1):
        fold = np.mod(sl*ax, 2*D0)
        zr = -D0 + np.where(fold <= D0, fold, 2*D0 - fold)
        seg = ((sl*ax) // D0).astype(int)
        d = Zcen - zr
        sig = sigma * (1.0 + skew*np.sign(d))
        amp = (ramp**b) * np.exp(-ax/decay)
        F = np.where(seg == b, np.maximum(F, amp*np.exp(-d**2/(2*sig**2))), F)
    return F


def fit(rms, win):
    """The identical block from iwbcurv."""
    pts = []
    for j in range(1, m + 1):
        col = rms[1:-1, j]
        if not np.isfinite(col).any() or col.max() <= 0:
            continue
        k = int(np.nanargmax(col)) + 1
        xv, zv = Xcen[k, j], Zcen[k, j]
        if not (win[0] <= abs(xv) <= win[1]):
            continue
        pts.append((abs(xv), zv))
    if len(pts) < 4:
        return None
    P = np.array(pts)
    A = np.vstack([P[:, 0], np.ones(len(P))]).T
    slope, icpt = np.linalg.lstsq(A, P[:, 1], rcond=None)[0]
    res = P[:, 1] - A @ [slope, icpt]
    return (np.degrees(np.arctan(abs(slope))), len(P),
            np.sqrt(np.mean(res**2)))




def window(theta_deg):
    bounce = D0 / np.tan(np.radians(theta_deg))
    return [max(0.15*bounce, 3*Lx/nx), min(0.70*bounce, 0.45*Lx)]


print("synthetic beam through the iwbcurv column-max + lstsq tracker")
print("grid: real TFI curvilinear, %dx%d, Lx=%.0f m, D0=%.0f m\n" % (nx, nz, Lx, D0))
print("%-7s %-9s %8s %10s %8s %6s %8s"
      % ('ratio', 'case', 'input', 'reported', 'bias', 'pts', 'fit rms'))
for ratio in (0.2, 0.4, 0.6, 0.8):
    th = np.degrees(np.arctan(np.sqrt(ratio**2/(1-ratio**2))))
    for case, nr in (('primary', 0), ('+2 refl', 2)):
        a, npts, rr = fit(synth(th, nrefl=nr), window(th))
        print("%-7s %-9s %8.2f %10.2f %+8.2f %6d %8.1f"
              % (ratio, case, th, a, a-th, npts, rr))

print("\nbeam width (ratio 0.6, theory 36.87, dz = %.1f m)\n" % (D0/(nz-1)))
print("%8s %7s %10s %8s" % ('sigma_m', 'cells', 'reported', 'bias'))
for sig in (200, 100, 60, 40, 25, 15, 10, 6):
    a, npts, rr = fit(synth(36.87, sigma=sig), window(36.87))
    print("%8.0f %7.1f %10.2f %+8.2f" % (sig, sig/(D0/(nz-1)), a, a-36.87))

print("\ncross-beam asymmetry (ratio 0.6, sigma = 25 m)\n")
print("%8s %10s %8s" % ('skew', 'reported', 'bias'))
for sk in (0.0, 0.1, 0.2, 0.4):
    a, npts, rr = fit(synth(36.87, sigma=25, skew=sk), window(36.87))
    print("%8.2f %10.2f %+8.2f" % (sk, a, a-36.87))

worst = max(abs(fit(synth(np.degrees(np.arctan(np.sqrt(r**2/(1-r**2)))))
                    , window(np.degrees(np.arctan(np.sqrt(r**2/(1-r**2))))))[0]
                - np.degrees(np.arctan(np.sqrt(r**2/(1-r**2)))))
            for r in (0.2, 0.4, 0.6, 0.8))
print("\nworst bias across the four frequencies: %+.2f deg" % worst)
assert worst < 0.15, "tracker bias exceeds 0.15 deg -- measurement is NOT clean"
print("tracker is clean (< 0.15 deg)")
