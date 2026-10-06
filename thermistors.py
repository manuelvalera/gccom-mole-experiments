#!/usr/bin/env python3
"""thermistors.py -- virtual thermistors at an isobath, as in Walter et al. (2012).

  python iwbcurv.py ... --mode1 3.7 --frames fr_bore --framerate 16
  python thermistors.py fr_bore --isobath 15 --mab 2 4 6 -o thermistors.png
  python thermistors.py run_mooring.npz -o thermistors.png   (from --moor)

The bottom row of their Figure 10 is temperature at 2, 4 and 6 m above bed at
the 15 m isobath -- the same quantity their moorings measured, which is what
lets the model be compared with the observations rather than with another
model. The signature they report for the non-canonical case is a gradual
decline over hours, then an abrupt warm front as the wave drains back downslope.

Buoyancy is converted to temperature with

    b = -g (rho - rho0) / rho0 = g alpha (T - T0),

so T = T0 + b / (g alpha), with alpha the thermal expansion coefficient. The
offset is arbitrary in a Boussinesq model with no absolute reference; what is
comparable is the SHAPE and the RANGE, so the series are also shown as
anomalies from their own means.
"""
import argparse
import glob
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('--isobath', type=float, default=15.0,
                    help='water depth at which to place the mooring, in metres')
    ap.add_argument('--mab', type=float, nargs='+', default=[2.0, 4.0, 6.0],
                    help='heights above bed, in metres')
    ap.add_argument('--alpha', type=float, default=2.0e-4,
                    help='thermal expansion coefficient, 1/degC')
    ap.add_argument('--T0', type=float, default=11.0,
                    help='reference temperature for the conversion, degC')
    ap.add_argument('-o', '--out', default='thermistors.png')
    ap.add_argument('--sponge', type=float, default=0.05,
                    help='shoreward sponge width as a fraction of the domain, used '
                         'only to warn when the mooring falls inside it')
    ap.add_argument('--hours', action='store_true',
                    help='x axis in hours rather than tidal periods')
    g = ap.parse_args()

    if g.dir.endswith('_mooring.npz'):
        # a mooring recorded every step by --moor: no frames to interpolate
        d = np.load(g.dir)
        T = float(d['T'])
        times = d['t']
        temps = g.T0 + d['b'] / (9.80665 * g.alpha)
        # a record that ends in a failure carries NaN at the tail; keep what
        # came before it rather than letting one bad sample void the ranges
        bad = (~np.isfinite(temps)).any(axis=1) | (np.abs(temps - temps[0]) > 20.0).any(axis=1)
        if bad.any():
            k = int(np.argmax(bad))
            print(f"record cut at t/T = {times[k] / float(d['T']):.2f}: the run was "
                  f"diverging from there ({len(times) - k} samples dropped)")
            times = times[:k]
            temps = temps[:k]
        g.mab = list(d['mab'])
        g.isobath = float(d['isobath'])
        x_moor = float(d['x'])
        dt = float(np.median(np.diff(times)))
        print(f"mooring record: {len(times)} samples, every {dt:.1f} s "
              f"({dt/60:.1f} min), at the {g.isobath:g} m isobath, x = {x_moor:.0f} m")
        plot(g, T, times, temps, x_moor)
        return
    grid = np.load(os.path.join(g.dir, 'grid.npz'))
    files = sorted(glob.glob(os.path.join(g.dir, 'frame_*.npz')))
    if not files:
        raise SystemExit(f'no frames in {g.dir}/')
    if 'b' not in np.load(files[0]):
        raise SystemExit('these frames carry no buoyancy; rerun with a solver '
                         'version that saves b')

    X, Z = grid['X'], grid['Z']
    T = float(grid['T'])
    n, m = Z.shape[0] - 1, Z.shape[1] - 1
    xc = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:] + X[1:, 1:])
    zc = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])

    # the bed, and the column where it crosses the requested isobath
    br = 0 if Z[0, :].mean() < Z[-1, :].mean() else -1
    bed_x = X[br, :]
    bed_z = Z[br, :]
    depth = -bed_z
    # walk from deep to shallow and take the first crossing
    order = np.argsort(bed_x)
    bx, bd = bed_x[order], depth[order]
    cross = np.where(np.diff(np.sign(bd - g.isobath)) != 0)[0]
    if cross.size == 0:
        raise SystemExit(f'the {g.isobath:g} m isobath is not in this domain '
                         f'(depths {bd.min():.1f} to {bd.max():.1f} m)')
    x_moor = float(np.interp(g.isobath, [bd[cross[-1] + 1], bd[cross[-1]]],
                             [bx[cross[-1] + 1], bx[cross[-1]]]))
    col = int(np.argmin(np.abs(xc[0, :] - x_moor)))
    bed_here = float(np.interp(x_moor, bx, bd))
    print(f"mooring at x = {x_moor:.0f} m, local depth {bed_here:.1f} m "
          f"(column {col} of {m})")
    # A mooring inside the sponge measures the sponge. The shoreward sponge is
    # usually a small fraction of the domain, so warn when the site is close to
    # that edge and suggest a deeper isobath, which sits further offshore.
    span = float(xc[0, :].max() - xc[0, :].min())
    from_edge = float(xc[0, :].max() - x_moor)
    if from_edge < g.sponge * span:
        print(f"   *** {from_edge:.0f} m from the shoreward edge, i.e. inside a "
              f"sponge of {g.sponge:g} Lx ({g.sponge*span:.0f} m) ***")
        print(f"   *** the series below is damped; use a deeper isobath, which "
              f"sits further offshore ***")

    # cell centres in that column, and the rows closest to each height above bed
    zcol = zc[:, col]
    rows = []
    for h in g.mab:
        ztarget = -bed_here + h
        r = int(np.argmin(np.abs(zcol - ztarget)))
        rows.append(r)
        print(f"   {h:g} mab -> z = {zcol[r]:.2f} m "
              f"({zcol[r] + bed_here:.2f} m above the bed)")

    times, series = [], []
    for f in files:
        fr = np.load(f)
        b = fr['b'].reshape(n + 2, m + 2)[1:-1, 1:-1]
        times.append(float(fr['t']))
        series.append([b[r, col] for r in rows])
    times = np.array(times)
    series = np.array(series)
    # b -> temperature
    temps = g.T0 + series / (9.80665 * g.alpha)
    plot(g, T, times, temps, x_moor)


def plot(g, T, times, temps, x_moor):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    tx = times / 3600.0 if g.hours else times / T
    xlabel = 'time (hours)' if g.hours else 't / T'
    fig, ax = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    colors = ['#8c510a', '#2166ac', '#b2182b', '#4d4d4d', '#1b7837']
    for i, h in enumerate(g.mab):
        ax[0].plot(tx, temps[:, i], color=colors[i % len(colors)],
                   label=f'{h:g} mab')
        ax[1].plot(tx, temps[:, i] - temps[:, i].mean(),
                   color=colors[i % len(colors)], label=f'{h:g} mab')
    ax[0].set_ylabel('T (degC)')
    ax[0].legend(loc='upper right', fontsize=9)
    ax[0].set_title(f'virtual thermistors at the {g.isobath:g} m isobath '
                    f'(x = {x_moor:.0f} m)\n'
                    f'compare Walter et al. (2012) Figure 10, bottom row')
    ax[1].set_ylabel('T anomaly (degC)')
    ax[1].set_xlabel(xlabel)
    ax[1].axhline(0.0, color='k', lw=0.5)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(g.out, dpi=120)
    plt.close(fig)

    rng = temps.max(axis=0) - temps.min(axis=0)
    print(f"\n{g.out}: {len(times)} samples, {tx[0]:.2f} to {tx[-1]:.2f} {xlabel}")
    for i, h in enumerate(g.mab):
        print(f"   {h:g} mab: range {rng[i]:.3f} degC")
    print("   Walter et al. observed bores of about 0.5 degC and warm-front "
          "relaxations of 1 degC or more")


if __name__ == '__main__':
    main()
