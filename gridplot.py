#!/usr/bin/env python3
"""gridplot.py -- draw the curvilinear grid, the ridge, and the side bulge.

  python gridplot.py --mole $env:MOLE_SRC                  # coarse, generated
  python gridplot.py --mole ... --nx 64 --nz 33 --Lx 6000
  python gridplot.py --npz paper-logs/bounce0.6.npz --skip 8

Three panels:
  (a) the whole domain.  The eta-lines bow because betaOf() uses a different
      stretching at the top and the bottom, and the side walls bulge by
      0.15 D0, so the xi-lines are NOT vertical -- that is what makes this a
      genuinely curvilinear grid rather than a sigma grid.
  (b) the ridge, at true aspect.  ab = 20 m on D0 = 1000 m is a 2% bump; with
      nz = 101 it is TWO CELLS tall, which is why raising --ab was the only
      handle on how well the source is resolved.
  (c) a side wall, showing the bulge.

Generating the grid (--mole) is usually clearer than reading a production npz:
at 512x201 the lines merge into a solid block, so --npz needs --skip.
"""
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def from_octave(mole, grids, nx, nz, Lx, D0, ab, lb, bt, bb):
    from oct2py import Oct2Py
    gd = os.path.join(grids, 'iwbridge')
    with open(os.path.join(gd, 'geom_over.m'), 'w') as fh:
        fh.write("function [Lx, D0, ab, Lb] = geom_over()\n"
                 f"    Lx = {Lx}; D0 = {D0}; ab = {ab}; Lb = {lb};\nend\n")
    oc = Oct2Py(); oc.eval("warning('off','all'); more off;")
    oc.addpath(mole); oc.addpath(os.path.abspath(gd))
    here = os.getcwd(); os.chdir(grids)
    try:
        oc.eval(f"global BT BB; BT={bt}; BB={bb};")
        oc.eval(f"cd('{os.path.abspath('.').replace(os.sep,'/')}'); "
                f"[X,Z]=gridGen('TFI','iwbridge',{nx},{nz},false);")
    finally:
        os.chdir(here)
    X = oc.pull('X'); Z = oc.pull('Z'); oc.exit()
    return X.T.copy(), Z.T.copy()


def wander(XM):
    """Max horizontal drift of a xi-line from bed to surface, in m."""
    return float(np.abs(XM[-1, :] - XM[0, :]).max())


def draw(ax, XM, ZM, skip, lw=0.5):
    for i in range(0, XM.shape[0], skip):
        ax.plot(XM[i, :], ZM[i, :], color='#2c6fb0', lw=lw)
    if (XM.shape[0]-1) % skip:
        ax.plot(XM[-1, :], ZM[-1, :], color='#2c6fb0', lw=lw)
    for j in range(0, XM.shape[1], skip):
        ax.plot(XM[:, j], ZM[:, j], color='#c0392b', lw=lw)
    if (XM.shape[1]-1) % skip:
        ax.plot(XM[:, -1], ZM[:, -1], color='#c0392b', lw=lw)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mole', default=os.environ.get('MOLE_SRC'))
    p.add_argument('--grids', default='grids')
    p.add_argument('--npz', default=None)
    p.add_argument('--nx', type=int, default=64)
    p.add_argument('--nz', type=int, default=33)
    p.add_argument('--Lx', type=float, default=6000.0)
    p.add_argument('--D0', type=float, default=1000.0)
    p.add_argument('--ab', type=float, default=20.0)
    p.add_argument('--lb', type=float, default=30.0)
    p.add_argument('--bt', type=float, default=1.6)
    p.add_argument('--bb', type=float, default=3.2)
    p.add_argument('--skip', type=int, default=1, help='draw every Nth line')
    p.add_argument('-o', '--out', default='grid.png')
    g = p.parse_args()

    if g.npz:
        d = np.load(g.npz)
        XM, ZM = d['X'], d['Z']
        src = os.path.basename(g.npz)
    else:
        if not g.mole or not os.path.isfile(os.path.join(g.mole, 'gridGen.m')):
            raise SystemExit("set --mole (or MOLE_SRC) to MOLE's src/matlab_octave")
        XM, ZM = from_octave(g.mole, g.grids, g.nx, g.nz,
                             g.Lx, g.D0, g.ab, g.lb, g.bt, g.bb)
        src = f"generated {g.nx}x{g.nz}"

    Lx = float(XM.max() - XM.min()); D0 = float(-ZM.min())
    nx, nz = XM.shape[1]-1, XM.shape[0]-1
    wd = wander(XM)
    crest = float(ZM[0, :].max() + D0)
    print(f"{src}: {nx}x{nz}  Lx={Lx:.0f} m  D0={D0:.0f} m")
    print(f"   xi-line wander bed->surface: {wd:.1f} m")
    print(f"   ridge crest height: {crest:.1f} m = {crest/(D0/nz):.1f} cells "
          f"(dz = {D0/nz:.1f} m)")

    fig = plt.figure(figsize=(13, 7), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0, :])
    draw(ax, XM, ZM, g.skip)
    ax.fill_between(XM[0, :], ZM[0, :], -D0*1.05, color='#6b5b45', zorder=3)
    ax.plot(XM[0, :], ZM[0, :], 'k', lw=1.4, zorder=4)
    ax.set_xlim(XM.min(), XM.max()); ax.set_ylim(-D0*1.05, 0.03*D0)
    ax.set_xlabel('$x$ (m)'); ax.set_ylabel('$z$ (m)')
    ax.set_title(f"(a) whole domain — {nx}x{nz}, xi-line wander {wd:.0f} m "
                 f"(0 would be a sigma grid)", fontsize=10)

    axr = fig.add_subplot(gs[1, 0])
    draw(axr, XM, ZM, 1, lw=0.7)
    axr.fill_between(XM[0, :], ZM[0, :], -D0*1.05, color='#6b5b45', zorder=3)
    axr.plot(XM[0, :], ZM[0, :], 'k', lw=1.8, zorder=4)
    zoom = max(6*g.lb, 200.0)
    axr.set_xlim(-zoom, zoom)
    axr.set_ylim(-D0 - 0.004*D0, -D0 + 6*crest)
    axr.set_aspect('equal', adjustable='box')
    axr.set_xlabel('$x$ (m)'); axr.set_ylabel('$z$ (m)')
    axr.set_title(f"(b) the ridge, true aspect — crest {crest:.0f} m "
                  f"= {crest/(D0/nz):.1f} cells tall", fontsize=10)

    axw = fig.add_subplot(gs[1, 1])
    draw(axw, XM, ZM, 1, lw=0.7)
    bulge = float(XM[:, 0].max() - XM[:, 0].min())
    axw.plot(XM[:, 0], ZM[:, 0], 'k', lw=2.0, zorder=5)
    axw.set_xlim(XM.min() - 0.01*Lx, XM.min() + 0.10*Lx)
    axw.set_ylim(-D0, 0)
    axw.set_xlabel('$x$ (m)'); axw.set_ylabel('$z$ (m)')
    axw.set_title(f"(c) left wall — bulges {bulge:.0f} m "
                  f"(0.15 D0, scaled by DEPTH not Lx)", fontsize=10)

    fig.suptitle("Fully curvilinear TFI grid: eta-lines bow because the top and "
                 "bottom stretching differ (betaOf), the walls bulge,\n"
                 "and the ridge is a 2% bump only a couple of cells tall  —  "
                 "blue = eta-lines, red = xi-lines", fontsize=11)
    fig.savefig(g.out, dpi=150, bbox_inches='tight')
    print(f"wrote {g.out}")


if __name__ == '__main__':
    main()
