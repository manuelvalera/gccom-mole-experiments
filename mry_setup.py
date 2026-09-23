#!/usr/bin/env python3
"""mry_setup.py -- turn the Monterey data into a model configuration.

  python mry_setup.py --bathy narrow_bathy_100m.csv --out grids/mryshelf

Reads the two data files and writes what the solver needs:

  grids/mryshelf/*.m   boundary curves for gridGen, with the measured depth
                       profile embedded (no file I/O from inside Octave)
  mry_N.txt            N(z) table, from the temperature fit
  mry_setup.png        what the configuration looks like, and whether the
                       internal tide is sub- or supercritical along it

WHAT THE DATA ARE

narrow_bathy_100m.csv is a 10.1 km transect in UTM zone 10N (EPSG 26910),
three columns 100 m apart, so effectively a 1-D depth profile, 6-88 m deep.
MRY_profile.txt is a seven-Gaussian fit of temperature against depth over the
same 0-110 m range. Together they describe a shallow shelf slab, NOT the deep
canyon -- the canyon is >1000 m and is not in these files.

WHY THE CONFIGURATION MATTERS MORE THAN THE GEOMETRY

At 36.8 N the M2 tide has f/omega = 0.62, so rotation cuts the beam slope by
22% and cannot be left out. N varies from 2.1e-3 near the surface to 8.9e-3
near the bottom, a factor of four, so a single scalar N is not usable either.
Both are reported below so the numbers are on record before any run.
"""
import argparse
import os

import numpy as np

# MRY_profile.txt: T(z) = sum a_i exp(-((z-b_i)/c_i)^2), z in m (negative down)
FIT_A = [12.69, 0.3666, 2.163, -0.2957, 1.343, 4.714, 0.5697]
FIT_B = [-2.649, -43.06, -56.65, -29.76, -70.12, -88.98, -25.25]
FIT_C = [85.04, 8.703, 14.45, 5.569, 11.93, 17.79, 11.71]

G = 9.81
ALPHA_T = 2.0e-4          # thermal expansion, 1/degC
OMEGA_M2 = 2*np.pi/44712.0
LAT = 36.8


def temperature(z):
    t = np.zeros_like(np.asarray(z, dtype=float))
    for a, b, c in zip(FIT_A, FIT_B, FIT_C):
        t = t + a*np.exp(-((z - b)/c)**2)
    return t


def buoyancy_frequency(z, floor=1e-4):
    """N(z) from the temperature fit with a linear equation of state.

    A floor is applied because the fit is only constrained over 0 to -110 m and
    a fitted profile can go weakly unstable where it is unconstrained; the floor
    keeps N real without inventing stratification.
    """
    t = temperature(z)
    n2 = G*ALPHA_T*np.gradient(t, z)
    return np.sqrt(np.maximum(n2, floor**2))


def read_transect(path):
    d = np.loadtxt(path)
    ux = np.unique(np.round(d[:, 0], 1))
    uy = np.unique(np.round(d[:, 1], 1))
    if len(ux)*len(uy) != len(d):
        raise SystemExit('bathymetry is not on a regular grid')
    z = d[:, 2].reshape(len(uy), len(ux))
    h = np.abs(np.nanmean(z, axis=1))          # average the three columns
    s = uy - uy.min()                          # along-transect distance, m
    order = np.argsort(s)
    return s[order], h[order]


def write_grid(outdir, s, h, smooth_m):
    """Boundary curves for gridGen, with the depth profile embedded."""
    os.makedirs(outdir, exist_ok=True)
    Lx = float(s.max() - s.min())
    D0 = float(h.max())

    # light smoothing: the solver's bed constraint differentiates this profile,
    # and 100 m data on a ~20 m grid would otherwise be a staircase
    if smooth_m > 0:
        w = max(int(round(smooth_m/np.median(np.diff(s)))), 1)
        k = np.ones(w)/w
        hs = np.convolve(np.pad(h, (w, w), mode='edge'), k, mode='same')[w:-w]
    else:
        hs = h

    xs = s - s.mean()                          # centre the domain on x = 0
    prof = '\n'.join(f'        {a:.4f} {b:.4f};' for a, b in zip(xs, hs))

    with open(os.path.join(outdir, 'geom.m'), 'w') as fh:
        fh.write("function [Lx, D0, ab, Lb] = geom()\n"
                 "% Monterey shelf transect. ab and Lb are unused here (the bed comes\n"
                 "% from measured bathymetry) but are kept so the signature matches the\n"
                 "% other grids and the driver's geom_over mechanism still works.\n"
                 "    if exist('geom_over') == 2\n"
                 "        [Lx, D0, ab, Lb] = geom_over();\n"
                 "    else\n"
                 f"        Lx = {Lx:.1f};\n"
                 f"        D0 = {D0:.2f};\n"
                 "        ab = 0;\n"
                 "        Lb = 1;\n"
                 "    end\n"
                 "end\n")

    with open(os.path.join(outdir, 'mry_depth.m'), 'w') as fh:
        fh.write("function d = mry_depth(x)\n"
                 "% Measured depth (positive, metres) against along-transect distance,\n"
                 "% centred on x = 0. Embedded rather than read from disk so gridGen\n"
                 "% works whatever the current directory is.\n"
                 "    persistent P\n"
                 "    if isempty(P)\n"
                 "        P = [\n" + prof + "\n        ];\n"
                 "    end\n"
                 "    d = interp1(P(:, 1), P(:, 2), x, 'linear', 'extrap');\n"
                 "end\n")

    with open(os.path.join(outdir, 'bottom.m'), 'w') as fh:
        fh.write("function XY = bottom(s)\n"
                 "    [Lx, D0, ab, Lb] = geom();\n"
                 "    x = -Lx/2 + Lx*stretch(s, 1, betaOf(s, 'bot'));\n"
                 "    XY = [x, -mry_depth(x)];\n"
                 "end\n")

    for name, body in (
        ('top.m', "    [Lx, D0, ab, Lb] = geom();\n"
                  "    x = -Lx/2 + Lx*stretch(s, 1, betaOf(s, 'top'));\n"
                  "    XY = [x, 0];\n"),
        ('left.m', "    [Lx, D0, ab, Lb] = geom();\n"
                   "    global BULGE\n"
                   "    if isempty(BULGE)\n"
                   "        BULGE = 0.0;\n"
                   "    end\n"
                   "    d = mry_depth(-Lx/2);\n"
                   "    XY = [-Lx/2 + BULGE*D0*sin(pi*s), -d*(1 - s)];\n"),
        ('right.m', "    [Lx, D0, ab, Lb] = geom();\n"
                    "    global BULGE\n"
                    "    if isempty(BULGE)\n"
                    "        BULGE = 0.0;\n"
                    "    end\n"
                    "    d = mry_depth(Lx/2);\n"
                    "    XY = [Lx/2 - BULGE*D0*sin(pi*s), -d*(1 - s)];\n"),
    ):
        with open(os.path.join(outdir, name), 'w') as fh:
            fh.write(f"function XY = {name[:-2]}(s)\n" + body + "end\n")

    for name in ('betaOf.m', 'stretch.m'):
        src = os.path.join('grids', 'iwbridge', name)
        if os.path.isfile(src):
            with open(src) as a, open(os.path.join(outdir, name), 'w') as b:
                b.write(a.read())
    return Lx, D0, xs, hs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bathy', default='narrow_bathy_100m.csv')
    ap.add_argument('--out', default=os.path.join('grids', 'mryshelf'))
    ap.add_argument('--trim', type=float, default=0.0,
                    help='drop the part of the transect shallower than this depth '
                         '(m). The shallow end is where the cells are worst: the '
                         'domain is 115:1 overall and its last kilometre is far '
                         'worse. Trimming at 20 m removes the most distorted '
                         'cells at the cost of about a kilometre of transect.')
    ap.add_argument('--smooth', type=float, default=200.0,
                    help='along-transect smoothing length in m (0 = none)')
    ap.add_argument('--ntable', default='mry_N.txt')
    ap.add_argument('--figure', default='mry_setup.png')
    g = ap.parse_args()

    s, h = read_transect(g.bathy)
    if g.trim > 0:
        keep = h >= g.trim
        if keep.sum() < 10:
            raise SystemExit(f'--trim {g.trim} leaves too little transect')
        s, h = s[keep], h[keep]
        print(f"trimmed to depths >= {g.trim:.0f} m: {len(h)} points, "
              f"{(s.max()-s.min())/1000:.1f} km")
    Lx, D0, xs, hs = write_grid(g.out, s, h, g.smooth)

    zt = np.linspace(-max(110.0, D0), 0.0, 1101)
    N = buoyancy_frequency(zt)
    np.savetxt(g.ntable, np.c_[zt, N], fmt='%12.4f %14.6e',
               header='z (m, negative down)   N (1/s)   from MRY_profile.txt, linear EOS')

    f = 2*7.292e-5*np.sin(np.radians(LAT))
    slope = np.sqrt(max(OMEGA_M2**2 - f**2, 0.0)/np.maximum(N**2 - OMEGA_M2**2, 1e-14))
    slope_norot = np.sqrt(OMEGA_M2**2/np.maximum(N**2 - OMEGA_M2**2, 1e-14))
    topo = np.abs(np.gradient(hs, xs))
    slope_at_bed = np.interp(-hs, zt, slope)
    crit = topo/np.maximum(slope_at_bed, 1e-12)

    print(f"\ntransect      {Lx/1000:.1f} km long, {hs.min():.1f} to {hs.max():.1f} m deep")
    print(f"grid files    {g.out}")
    print(f"N(z)          {N.min():.2e} to {N.max():.2e} 1/s -> {g.ntable}")
    print(f"M2            omega {OMEGA_M2:.3e}, f {f:.3e}, f/omega {f/OMEGA_M2:.2f}")
    print(f"beam slope    {np.percentile(slope, 5):.4f} to {np.percentile(slope, 95):.4f} "
          f"(without rotation: {np.percentile(slope_norot, 5):.4f} to "
          f"{np.percentile(slope_norot, 95):.4f})")
    print(f"criticality   median {np.median(crit):.2f}, max {crit.max():.2f}, "
          f"{100*np.mean(crit > 1):.0f}% of the transect supercritical")
    print(f"bounce        {D0/np.interp(-D0, zt, slope)/1000:.1f} km at the deepest point "
          f"-> {Lx/(D0/np.interp(-D0, zt, slope)):.1f} bounces across the domain\n")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(11, 6.5))
    ax[0, 0].plot(xs/1000, -hs, color='saddlebrown')
    ax[0, 0].fill_between(xs/1000, -hs, -D0*1.05, color='saddlebrown', alpha=0.3)
    ax[0, 0].set_xlabel('along transect (km)')
    ax[0, 0].set_ylabel('z (m)')
    ax[0, 0].set_title(f'Monterey shelf transect ({Lx/1000:.1f} km)')
    ax[0, 1].semilogx(N, zt, color='tab:blue')
    ax[0, 1].set_xlabel('N (1/s)')
    ax[0, 1].set_ylabel('z (m)')
    ax[0, 1].set_title('stratification from the temperature fit')
    ax[1, 0].plot(xs/1000, topo, color='k', label='bathymetric slope')
    ax[1, 0].plot(xs/1000, slope_at_bed, color='tab:blue', label='M2 beam slope')
    ax[1, 0].set_yscale('log')
    ax[1, 0].set_xlabel('along transect (km)')
    ax[1, 0].set_title('slopes')
    ax[1, 0].legend(fontsize=8)
    ax[1, 1].plot(xs/1000, crit, color='tab:red')
    ax[1, 1].axhline(1, color='k', ls='--', lw=1)
    ax[1, 1].set_yscale('log')
    ax[1, 1].set_xlabel('along transect (km)')
    ax[1, 1].set_title('criticality (above 1 = supercritical, reflects)')
    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(g.figure, dpi=130)
    print(f"figure        {g.figure}")


if __name__ == '__main__':
    main()
