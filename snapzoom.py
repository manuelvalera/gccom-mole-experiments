#!/usr/bin/env python3
"""snapzoom.py -- look at the field in a small window over the last few frames.

  python iwbcurv.py ... --frames fr_fail --framerate 600 --framefrom 2.85
  python snapzoom.py fr_fail --x 8800 9500 --z -45 -10 --last 6 -o zoom.png

Maxima and cell counts say THAT something grows and WHERE. They cannot say what
it looks like, and that is the question that separates physics from numerics:

  a coherent overturn a few metres across, dense fluid rolling over light,
  is a real feature the model may be failing to mix;

  a checkerboard at the grid scale -- alternating signs cell to cell -- is a
  numerical mode, and no physical parameterisation will cure it.

Two rows per frame: total buoyancy anomaly (background removed, so overturns
show as closed contours) and vertical velocity, both on the true curvilinear
grid, with the cell edges drawn so the grid scale is visible.
"""
import argparse
import glob
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('--x', type=float, nargs=2, required=True, help='x window, m')
    ap.add_argument('--z', type=float, nargs=2, required=True, help='z window, m')
    ap.add_argument('--last', type=int, default=6, help='how many final frames to show')
    ap.add_argument('--every', type=int, default=1, help='take every Nth of those')
    ap.add_argument('--edges', action='store_true', help='draw the cell edges')
    ap.add_argument('-o', '--out', default='zoom.png')
    g = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import colors

    grid = np.load(os.path.join(g.dir, 'grid.npz'))
    X, Z = grid['X'], grid['Z']
    n, m = Z.shape[0] - 1, Z.shape[1] - 1
    files = sorted(glob.glob(os.path.join(g.dir, 'frame_*.npz')))
    files = files[-g.last::max(g.every, 1)]
    if not files:
        raise SystemExit(f'no frames in {g.dir}/')

    # columns inside the x window, so the plot does not render 4000 columns
    xc_row = 0.5 * (X[0, :-1] + X[0, 1:])
    cols = np.where((xc_row >= g.x[0]) & (xc_row <= g.x[1]))[0]
    if cols.size == 0:
        raise SystemExit('the x window contains no cells')
    c0, c1 = int(cols.min()), int(cols.max()) + 1
    Xw = X[:, c0:c1 + 1]
    Zw = Z[:, c0:c1 + 1]

    fig, axes = plt.subplots(2, len(files), figsize=(3.6 * len(files), 6.4),
                             squeeze=False, sharey=True)
    for k, f in enumerate(files):
        fr = np.load(f)
        b = fr['b'].reshape(n + 2, m + 2)[1:-1, 1:-1][:, c0:c1]
        w = fr['w'].reshape(n + 2, m + 2)[1:-1, 1:-1][:, c0:c1]
        for row, (fld, cmap, lab) in enumerate(((b, 'RdBu_r', 'b anomaly'),
                                                (w, 'PuOr_r', 'w'))):
            ax = axes[row, k]
            lim = float(np.percentile(np.abs(fld), 99.5)) or 1e-12
            pc = ax.pcolormesh(Xw, Zw, fld, cmap=cmap,
                               norm=colors.Normalize(-lim, lim), shading='flat',
                               edgecolors='k' if g.edges else 'face',
                               linewidth=0.05 if g.edges else 0)
            ax.set_xlim(*g.x)
            ax.set_ylim(*g.z)
            if row == 0:
                ax.set_title(f"t/T = {float(fr['tT']):.4f}", fontsize=9)
            if k == 0:
                ax.set_ylabel(f'{lab}\nz (m)')
            if row == 1:
                ax.set_xlabel('x (m)')
            fig.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
    br = 0 if Z[0, :].mean() < Z[-1, :].mean() else -1
    for ax in axes.ravel():
        ax.plot(X[br, :], Z[br, :], 'k', lw=1)
    fig.suptitle(f'{os.path.basename(os.path.normpath(g.dir))}: last {len(files)} frames, '
                 f'x {g.x[0]:.0f}-{g.x[1]:.0f} m', fontsize=10)
    fig.tight_layout()
    fig.savefig(g.out, dpi=110)
    plt.close(fig)
    print(f'{g.out}: {len(files)} frames, {c1 - c0} columns in the window')


if __name__ == '__main__':
    main()
