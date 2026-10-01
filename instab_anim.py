#!/usr/bin/env python3
"""instab_anim.py -- watch where the trouble starts.

  python iwbcurv.py ... --frames fr_fail --framerate 120 --framefrom 3.3
  python instab_anim.py fr_fail -o fail.gif --x 7000 10100 --z -60 0
  python instab_anim.py fr_fail -o fail.gif --sheet fail_sheet.png --last 6
  python instab_anim.py fr_fail -o fail.gif --n2=-1e-5     (note the = sign)

Every run-up failure so far has been read through maxima and cell counts. This
draws the field and marks, frame by frame,

  red dots     cells where the TOTAL stratification is unstable, N^2 < 0:
               dense water over light, an overturn;
  yellow x     cells where the local advective CFL exceeds --cfl (default
               0.5), from the contravariant velocities and the grid metrics;
  white +      the location of max|u|.

Whether the markers first appear at the bed, at the surface, in the thin
shallow cells, in the sponge, or at the nose of the surge -- and which of the
two kinds appears first -- separates a physical overturn from a numerical
limit.

The frames must come from a solver version that stores dt and N2c in
grid.npz; older frame directories are refused with a message.
"""
import argparse
import glob
import os

import numpy as np


def centre_metrics(X, Z):
    """xi_x, xi_z, eta_x, eta_z at cell centres, index-space derivatives."""
    x_xi = np.gradient(X, axis=1)
    x_et = np.gradient(X, axis=0)
    z_xi = np.gradient(Z, axis=1)
    z_et = np.gradient(Z, axis=0)
    J = x_xi * z_et - x_et * z_xi
    avg = lambda A: 0.25 * (A[:-1, :-1] + A[1:, :-1] + A[:-1, 1:] + A[1:, 1:])
    Jc = avg(J)
    return avg(z_et) / Jc, -avg(x_et) / Jc, -avg(z_xi) / Jc, avg(x_xi) / Jc


def diagnose(fr, n, m, met, N2bg, dt):
    """Total N^2 and local CFL at interior cell centres."""
    xix, xiz, etx, etz = met
    B = fr['b'].reshape(n + 2, m + 2)
    dbx = 0.5 * (B[1:-1, 2:] - B[1:-1, :-2])
    dbe = 0.5 * (B[2:, 1:-1] - B[:-2, 1:-1])
    N2 = N2bg + xiz * dbx + etz * dbe
    u = fr['u'].reshape(n + 2, m + 2)[1:-1, 1:-1]
    w = fr['w'].reshape(n + 2, m + 2)[1:-1, 1:-1]
    cfl = dt * np.maximum(np.abs(u * xix + w * xiz), np.abs(u * etx + w * etz))
    return N2, cfl, B[1:-1, 1:-1], u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('-o', '--out', default='instab.gif')
    ap.add_argument('--x', type=float, nargs=2, default=None, help='x window, m')
    ap.add_argument('--z', type=float, nargs=2, default=None, help='z window, m')
    ap.add_argument('--field', default='b', choices=['b', 'w', 'u'])
    ap.add_argument('--cfl', type=float, default=0.5, help='CFL marker threshold')
    ap.add_argument('--n2', type=float, default=-1e-6,
                    help='mark cells with total N^2 below this (1/s^2). Slightly below '
                         'zero, so round-off-level noise does not paint the whole domain. '
                         'Negative values need the = form: --n2=-1e-5')
    ap.add_argument('--fps', type=float, default=8.0)
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--sheet', default=None,
                    help='also write a PNG of the last --last frames side by side')
    ap.add_argument('--last', type=int, default=6)
    ap.add_argument('--dpi', type=int, default=100)
    g = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation, colors

    grid = np.load(os.path.join(g.dir, 'grid.npz'))
    if 'dt' not in grid.files or 'N2c' not in grid.files:
        raise SystemExit('grid.npz has no dt/N2c: rerun with the current iwbcurv.py '
                         'to write frames this script can read')
    X, Z = grid['X'], grid['Z']
    n, m = Z.shape[0] - 1, Z.shape[1] - 1
    dt = float(grid['dt'])
    N2bg = grid['N2c'].reshape(n + 2, m + 2)[1:-1, 1:-1]
    met = centre_metrics(X, Z)
    xc = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:] + X[1:, 1:])
    zc = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])
    files = sorted(glob.glob(os.path.join(g.dir, 'frame_*.npz')))[::max(g.stride, 1)]
    if not files:
        raise SystemExit(f'no frames in {g.dir}/')
    # a failure dump ends with frames that are already NaN; they show nothing and
    # swamp the colour scale, so stop at the last fully finite frame
    good = [f for f in files if np.isfinite(np.load(f)['b']).all()]
    if len(good) < len(files):
        print(f"skipping {len(files) - len(good)} non-finite frame(s) at the end")
    files = good
    if not files:
        raise SystemExit('every frame is non-finite')
    xwin = g.x if g.x else (float(X.min()), float(X.max()))
    zwin = g.z if g.z else (float(Z.min()), 0.0)
    br = 0 if Z[0, :].mean() < Z[-1, :].mean() else -1

    # one colour scale for the whole animation, from the last frame
    last = np.load(files[-1])
    _, _, Bl, _ = diagnose(last, n, m, met, N2bg, dt)
    fld_of = {'b': lambda fr: fr['b'].reshape(n + 2, m + 2)[1:-1, 1:-1],
              'w': lambda fr: fr['w'].reshape(n + 2, m + 2)[1:-1, 1:-1],
              'u': lambda fr: fr['u'].reshape(n + 2, m + 2)[1:-1, 1:-1]}
    sel = (xc >= xwin[0]) & (xc <= xwin[1]) & (zc >= zwin[0]) & (zc <= zwin[1])
    ref = fld_of[g.field](last)[sel]
    ref = ref[np.isfinite(ref)]
    lim = float(np.percentile(np.abs(ref), 99.0)) if ref.size else 1.0
    lim = lim or 1e-12

    def draw(ax, fr):
        N2, cfl, _, u = diagnose(fr, n, m, met, N2bg, dt)
        f = fld_of[g.field](fr)
        ax.clear()
        ax.pcolormesh(X, Z, np.nan_to_num(f), cmap='RdBu_r',
                      norm=colors.Normalize(-lim, lim), shading='flat')
        ax.fill_between(X[br, :], Z[br, :], Z[br, :].min() - 5, color='saddlebrown',
                        zorder=3)
        unst = (N2 < g.n2) & sel
        hot = (cfl > g.cfl) & sel
        ax.scatter(xc[unst], zc[unst], s=6, c='red', marker='o', lw=0, zorder=4,
                   label=f'N^2<{g.n2:g} ({int(unst.sum())})')
        ax.scatter(xc[hot], zc[hot], s=18, c='yellow', marker='x', lw=0.8, zorder=5,
                   label=f'CFL>{g.cfl:g} ({int(hot.sum())})')
        uu = np.where(np.isfinite(u), np.abs(u), -1.0)
        k = np.unravel_index(np.argmax(uu), uu.shape)
        ax.plot(xc[k], zc[k], '+', color='white', ms=14, mew=2, zorder=6)
        ax.set_xlim(*xwin)
        ax.set_ylim(*zwin)
        cflmax = float(np.nanmax(cfl[sel])) if sel.any() else 0.0
        ax.set_title(f"t/T = {float(fr['tT']):.3f}   N^2<{g.n2:g}: {int(unst.sum())} cells   "
                     f"max CFL {cflmax:.2f}   white + = max|u|", fontsize=9)
        ax.legend(loc='lower left', fontsize=7, framealpha=0.7)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('z (m)')

    fig, ax = plt.subplots(figsize=(11, 4.4))
    anim = animation.FuncAnimation(fig, lambda i: draw(ax, np.load(files[i])),
                                   frames=len(files), blit=False)
    anim.save(g.out, writer=animation.PillowWriter(fps=g.fps), dpi=g.dpi)
    plt.close(fig)
    print(f"{g.out}: {len(files)} frames, field {g.field}, "
          f"window x {xwin[0]:.0f}..{xwin[1]:.0f} z {zwin[0]:.0f}..{zwin[1]:.0f}")

    if g.sheet:
        pick = files[-g.last:]
        fig, axes = plt.subplots(len(pick), 1, figsize=(11, 3.0 * len(pick)))
        axes = np.atleast_1d(axes)
        for a_, f_ in zip(axes, pick):
            draw(a_, np.load(f_))
        fig.tight_layout()
        fig.savefig(g.sheet, dpi=g.dpi)
        plt.close(fig)
        print(f"{g.sheet}: last {len(pick)} frames stacked")

    # text summary: when and where each kind of marker first appears
    print("\nfirst appearance in the window:")
    first_unst = first_hot = None
    for f_ in files:
        fr = np.load(f_)
        N2, cfl, _, _ = diagnose(fr, n, m, met, N2bg, dt)
        if first_unst is None and ((N2 < -1e-5) & sel).any():
            idx = np.argwhere((N2 < -1e-5) & sel)[0]
            first_unst = (float(fr['tT']), float(xc[tuple(idx)]), float(zc[tuple(idx)]))
        if first_hot is None and ((cfl > g.cfl) & sel).any():
            idx = np.argwhere((cfl > g.cfl) & sel)[0]
            first_hot = (float(fr['tT']), float(xc[tuple(idx)]), float(zc[tuple(idx)]))
    for name, v in (('N^2 < -1e-5', first_unst), (f'CFL > {g.cfl:g}', first_hot)):
        if v:
            print(f"   {name:12s} t/T = {v[0]:.3f} at x = {v[1]:.0f} m, z = {v[2]:.1f} m")
        else:
            print(f"   {name:12s} never")


if __name__ == '__main__':
    main()
