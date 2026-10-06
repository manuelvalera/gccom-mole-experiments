#!/usr/bin/env python3
"""grid_view.py -- draw the grid a run actually used.

  python grid_view.py fr_walter_xi2_a10\\grid.npz -o grid.png
  python grid_view.py fr_walter_xi2_a10\\grid.npz --zoom 7800 10000 --skew -o step.png
  python grid_view.py .gridcache\\grid_<hash>.npz --lsl 0.20 --lslr 0.005 --isobath 15

Reads any .npz holding the node coordinates X and Z: the grid.npz every run
with --frames writes, or an entry in the solver's grid cache. The orientation
is detected from the coordinates themselves, so either layout works.

Panels: the whole domain, and optionally a zoomed window. --skew colours each
cell by how far its corner angle departs from 90 degrees -- the quantity that
went to 60-80 degrees with the default stretching on the shelf and stays near
the bottom slope's own angle on a sigma grid. --lsl/--lslr shade the offshore
(left) and shoreward (right) sponges as fractions of the domain length, and
--isobath marks where a mooring on that depth contour would sit.
"""
import argparse

import numpy as np


def load_grid(path):
    d = np.load(path)
    X = np.asarray(d['X'], float)
    Z = np.asarray(d['Z'], float)
    # rows should run from bed to surface and columns along x; whichever axis
    # X varies along is the horizontal one
    if np.ptp(X, axis=0).mean() > np.ptp(X, axis=1).mean():
        X = X.T
        Z = Z.T
    if Z[0, :].mean() > Z[-1, :].mean():
        X = X[::-1, :]
        Z = Z[::-1, :]
    return X, Z


def cell_skew(X, Z):
    """Deviation of each cell's lower-left corner angle from 90 degrees."""
    ex = X[:-1, 1:] - X[:-1, :-1]
    ez = Z[:-1, 1:] - Z[:-1, :-1]
    fx = X[1:, :-1] - X[:-1, :-1]
    fz = Z[1:, :-1] - Z[:-1, :-1]
    cosang = (ex * fx + ez * fz) / (np.hypot(ex, ez) * np.hypot(fx, fz))
    ang = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
    return np.abs(ang - 90.0)


def draw(ax, X, Z, xlim, sx, sz, skew, g, title):
    from matplotlib import colors
    x0, x1 = xlim
    cols = np.where((X[0, :] >= x0) & (X[0, :] <= x1))[0]
    c0 = max(int(cols.min()) - 1, 0)
    c1 = min(int(cols.max()) + 1, X.shape[1] - 1)
    if skew is not None:
        pc = ax.pcolormesh(X[:, c0:c1 + 1], Z[:, c0:c1 + 1], skew[:, c0:c1],
                           cmap='magma_r', norm=colors.Normalize(0, max(10.0, float(skew.max()))),
                           shading='flat')
        ax.figure.colorbar(pc, ax=ax, fraction=0.03, pad=0.01,
                           label='degrees from orthogonal')
    for j in range(0, Z.shape[0], sz):
        ax.plot(X[j, c0:c1 + 1], Z[j, c0:c1 + 1], color='0.45', lw=0.4)
    for i in range(c0, c1 + 1, sx):
        ax.plot(X[:, i], Z[:, i], color='0.45', lw=0.4)
    zb = Z[0, :]
    ax.fill_between(X[0, :], zb, zb.min() - 10.0, color='saddlebrown', zorder=3)
    ax.plot(X[0, :], zb, 'k', lw=1.3, zorder=4)
    Lx = float(X.max() - X.min())
    if g.lsl > 0:
        ax.axvspan(X.min(), X.min() + g.lsl * Lx, color='tab:blue', alpha=0.18, zorder=2)
    if g.lslr > 0:
        ax.axvspan(X.max() - g.lslr * Lx, X.max(), color='tab:orange', alpha=0.35, zorder=2)
    if g.isobath:
        h = -zb
        k = np.where(np.diff(np.sign(h - g.isobath)) != 0)[0]
        if k.size:
            xm = float(X[0, k[-1]])
            ax.plot(xm, -g.isobath + 1.5, 'v', color='red', ms=9, zorder=5)
            ax.annotate(f'{g.isobath:g} m', (xm, -g.isobath + 3.0), color='red',
                        ha='center', fontsize=8, zorder=5)
    ax.set_xlim(x0, x1)
    ax.set_ylim(float(Z[0, c0:c1 + 1].min()) - 3.0, 2.0)
    ax.set_ylabel('z (m)')
    ax.set_title(title, fontsize=10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('grid', help='an .npz with node coordinates X and Z')
    ap.add_argument('-o', '--out', default='grid.png')
    ap.add_argument('--zoom', type=float, nargs=2, default=None, metavar=('X0', 'X1'),
                    help='add a second panel over this x range, m')
    ap.add_argument('--every', type=int, nargs=2, default=None, metavar=('SX', 'SZ'),
                    help='draw every SX-th column and SZ-th row (default: about 50 of each)')
    ap.add_argument('--skew', action='store_true', help='colour cells by skewness')
    ap.add_argument('--lsl', type=float, default=0.0, help='offshore sponge, fraction of Lx')
    ap.add_argument('--lslr', type=float, default=0.0, help='shoreward sponge, fraction of Lx')
    ap.add_argument('--isobath', type=float, default=None, help='mark a mooring depth, m')
    ap.add_argument('--dpi', type=int, default=120)
    g = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    X, Z = load_grid(g.grid)
    nz, nx = Z.shape[0] - 1, Z.shape[1] - 1
    sk = cell_skew(X, Z)
    h = -Z[0, :]
    dz = (Z[-1, :] - Z[0, :]) / nz
    worst = np.unravel_index(int(np.argmax(sk)), sk.shape)
    print(f"{g.grid}: {nx} x {nz} cells, x {X.min():.0f} to {X.max():.0f} m, "
          f"depth {h.min():.1f} to {h.max():.1f} m")
    print(f"   offshore (deep) end on the {'left' if h[0] > h[-1] else 'RIGHT'}")
    print(f"   skewness: median {np.median(sk):.1f} deg, worst {sk.max():.1f} deg at "
          f"x = {X[worst]:.0f} m, z = {Z[worst]:.1f} m")
    print(f"   cell thickness {dz.min()*100:.1f} cm (shallowest column) to {dz.max():.2f} m; "
          f"cell width {np.median(np.diff(X[0, :])):.1f} m")

    sx, sz = g.every if g.every else (max(nx // 50, 1), max(nz // 15, 1))
    npanel = 2 if g.zoom else 1
    fig, axes = plt.subplots(npanel, 1, figsize=(11, 4.2 * npanel), squeeze=False)
    draw(axes[0, 0], X, Z, (float(X.min()), float(X.max())), sx, sz,
         sk if g.skew else None, g,
         f'{nx} x {nz} grid, lines drawn every {sx} columns and {sz} rows')
    if g.zoom:
        zsx = max(sx // 8, 1)
        zsz = max(sz // 2, 1)
        draw(axes[1, 0], X, Z, tuple(g.zoom), zsx, zsz, sk if g.skew else None, g,
             f'x {g.zoom[0]:.0f} to {g.zoom[1]:.0f} m, lines every {zsx} columns and {zsz} rows')
    axes[-1, 0].set_xlabel('x (m)')
    fig.tight_layout()
    fig.savefig(g.out, dpi=g.dpi)
    plt.close(fig)
    print(f"   -> {g.out}")


if __name__ == '__main__':
    main()
