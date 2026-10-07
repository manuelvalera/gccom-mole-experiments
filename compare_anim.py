#!/usr/bin/env python3
"""compare_anim.py -- two runs, same moments, with their 15 m mooring records.

  python compare_anim.py fr_walter_xi2_anim fr_walter_xi02_anim `
      --labels "steep step, xi ~ 2" "gentle slope, xi ~ 0.2" `
      --moorings mry-fig10/walter_xi2_full_mooring.npz mry-fig10/walter_xi02_full_mooring.npz `
      --x 6500 10000 --z -60 0 -o docs/figures/bores.gif

Top panels: temperature near the shore in each run, with isotherms, so cold
water can be seen running up each slope and draining back. Bottom panel: the
2 m above bed record of each mooring with a cursor at the current time, so the
canonical and non-canonical shapes measured by events.py are seen as motion.

Temperature is the background profile plus the perturbation, converted with a
linear equation of state: T = T0 + (b_background + b) / (g alpha). The
background buoyancy is integrated up each column from the stored N^2, so the
frames must come from a solver version that writes N2c into grid.npz.

Frames are paired by time, nearest t/T, so the two runs need not have been
written at exactly the same steps.
"""
import argparse
import glob
import os

import numpy as np


def load_run(d):
    grid = np.load(os.path.join(d, 'grid.npz'))
    if 'N2c' not in grid.files:
        raise SystemExit(f'{d}/grid.npz has no N2c: rerun with the current iwbcurv.py')
    X = grid['X']
    Z = grid['Z']
    n = Z.shape[0] - 1
    m = Z.shape[1] - 1
    zc = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])
    N2 = grid['N2c'].reshape(n + 2, m + 2)[1:-1, 1:-1]
    # Background buoyancy is a function of absolute depth, not of height above
    # each column's own bed: integrate N^2 up the deepest column once, then
    # read every cell off that profile by its z. (Integrating column by column
    # would put the same temperature on every bed, 5 m deep or 81.)
    jd = int(np.argmin(Z[0, :]))
    zref = zc[:, jd]
    n2ref = N2[:, jd]
    bref = np.zeros_like(zref)
    bref[1:] = np.cumsum(0.5 * (n2ref[1:] + n2ref[:-1]) * np.diff(zref))
    bbg = np.interp(zc, zref, bref)
    files = sorted(glob.glob(os.path.join(d, 'frame_*.npz')))
    if not files:
        raise SystemExit(f'no frames in {d}/')
    times = np.array([float(np.load(f)['tT']) for f in files])
    return dict(X=X, Z=Z, n=n, m=m, bbg=bbg, zref=zref, bref=bref, zc=zc,
                files=files, times=times)


def temperature(run, f, T0, alpha):
    b = np.load(f)['b'].reshape(run['n'] + 2, run['m'] + 2)[1:-1, 1:-1]
    return T0 + (run['bbg'] + b) / (9.80665 * alpha)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs=2, help='two frame directories')
    ap.add_argument('--labels', nargs=2, default=None)
    ap.add_argument('--moorings', nargs=2, default=None,
                    help='the two runs\' _mooring.npz records, for the bottom panel')
    ap.add_argument('--x', type=float, nargs=2, default=(6500.0, 10000.0))
    ap.add_argument('--z', type=float, nargs=2, default=(-60.0, 0.0))
    ap.add_argument('--T0', type=float, default=11.0)
    ap.add_argument('--alpha', type=float, default=2.0e-4)
    ap.add_argument('--fps', type=float, default=10.0)
    ap.add_argument('--stride', type=int, default=1)
    ap.add_argument('--dpi', type=int, default=90)
    ap.add_argument('-o', '--out', default='bores.gif')
    g = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation

    runs = [load_run(d) for d in g.runs]
    labels = g.labels if g.labels else [os.path.basename(os.path.normpath(d)) for d in g.runs]
    t0 = max(r['times'].min() for r in runs)
    t1 = min(r['times'].max() for r in runs)
    master = runs[0]['times']
    master = master[(master >= t0) & (master <= t1)][::max(g.stride, 1)]
    if master.size == 0:
        raise SystemExit('the two runs have no overlapping times')

    # one colour scale for the whole animation, from the cells actually on
    # screen in the first pair -- the deep water off-screen would otherwise
    # stretch the scale and saturate the window
    vals = []
    for r in runs:
        xc_ = 0.5 * (r['X'][:-1, :-1] + r['X'][:-1, 1:])
        win = ((xc_ >= g.x[0]) & (xc_ <= g.x[1]) & (r['zc'] >= g.z[0]) & (r['zc'] <= g.z[1]))
        for tT in (master[0], master[len(master) // 2], master[-1]):
            T = temperature(r, r['files'][int(np.argmin(np.abs(r['times'] - tT)))],
                            g.T0, g.alpha)
            vals.append(T[win])
    vals = np.concatenate(vals)
    lo = float(np.percentile(vals, 1))
    hi = float(np.percentile(vals, 99))
    levels = np.linspace(lo, hi, 9)

    moor = []
    if g.moorings:
        for p, r in zip(g.moorings, runs):
            d = np.load(p)
            tt = d['t'] / float(d['T'])
            # full temperature at the sensor, on the same background as the
            # panels: the record itself holds only the perturbation
            zs = -float(d['isobath']) + float(d['mab'][0])
            bg = float(np.interp(zs, r['zref'], r['bref']))
            temp = g.T0 + (bg + d['b'][:, 0]) / (9.80665 * g.alpha)
            ok = np.isfinite(temp)
            moor.append((tt[ok], temp[ok], float(d['x'])))

    nrow = 3 if moor else 2
    fig = plt.figure(figsize=(10, 3.0 * nrow))
    gs = fig.add_gridspec(nrow, 1, height_ratios=[1, 1, 0.8][:nrow])
    axes = [fig.add_subplot(gs[k]) for k in range(nrow)]

    def draw(k):
        tT = float(master[k])
        for ax, r, lab in zip(axes[:2], runs, labels):
            j = int(np.argmin(np.abs(r['times'] - tT)))
            T = temperature(r, r['files'][j], g.T0, g.alpha)
            ax.clear()
            sel = (r['X'][0, :] >= g.x[0] - 200) & (r['X'][0, :] <= g.x[1] + 200)
            c0 = max(int(np.argmax(sel)) - 1, 0)
            c1 = min(len(sel) - 1 - int(np.argmax(sel[::-1])) + 1, r['m'])
            Xw = r['X'][:, c0:c1 + 1]
            Zw = r['Z'][:, c0:c1 + 1]
            Tw = T[:, c0:c1]
            ax.pcolormesh(Xw, Zw, Tw, cmap='RdYlBu_r', vmin=lo, vmax=hi, shading='flat')
            xc = 0.5 * (Xw[:-1, :-1] + Xw[:-1, 1:])
            zc = 0.25 * (Zw[:-1, :-1] + Zw[1:, :-1] + Zw[:-1, 1:] + Zw[1:, 1:])
            ax.contour(xc, zc, Tw, levels=levels, colors='k', linewidths=0.5)
            xb = r['X'][0, :]
            zb = r['Z'][0, :]
            ax.fill_between(xb, zb, zb.min() - 10, color='saddlebrown', zorder=3)
            if moor:
                xm = moor[axes.index(ax)][2]
                ax.plot([xm, xm], [float(np.interp(xm, xb, zb)), 0.0], color='k', lw=0.8, ls=':', zorder=4)
            ax.set_xlim(*g.x)
            ax.set_ylim(*g.z)
            ax.set_ylabel('z (m)')
            ax.set_title(f'{lab}    t/T = {tT:.2f}', fontsize=10)
        axes[1].set_xlabel('x (m)')
        if moor:
            ax = axes[2]
            ax.clear()
            for (tt, temp, xm), lab, col in zip(moor, labels, ('tab:red', 'tab:blue')):
                ax.plot(tt, temp, color=col, lw=0.9, label=lab)
            ax.axvline(tT, color='k', lw=1)
            ax.set_xlim(t0, t1)
            ax.set_xlabel('t / T')
            ax.set_ylabel('T at 2 mab (degC)')
            ax.legend(fontsize=8, loc='lower left')
            ax.set_title('15 m mooring, dotted line in the panels above', fontsize=9)
        fig.tight_layout()

    anim = animation.FuncAnimation(fig, draw, frames=len(master), blit=False)
    if g.out.lower().endswith('.mp4'):
        writer = animation.FFMpegWriter(fps=g.fps)
    else:
        writer = animation.PillowWriter(fps=g.fps)
    anim.save(g.out, writer=writer, dpi=g.dpi)
    plt.close(fig)
    print(f'{g.out}: {len(master)} frames, t/T {t0:.2f} to {t1:.2f}, '
          f'temperature {lo:.2f} to {hi:.2f} degC')


if __name__ == '__main__':
    main()
