#!/usr/bin/env python3
"""ctd_profile.py -- an N(z) table from a dated CTD cast, for --nprofile.

  python ctd_profile.py C1_CTD.mat --date 2010-05-18 -o mry_N_c1_20100518.txt
  python ctd_profile.py C1_CTD.mat --list-month 5          (what casts exist)

Reads MBARI's C1 station casts (a MATLAB struct CTD.cast with mtime, depth,
temp, sal, dens) and computes the buoyancy frequency from the measured density,

    N^2 = (g / rho0) d(rho)/d(depth),

after light smoothing, since a 1 m CTD profile differenced directly is noisy
and goes spuriously negative. N is floored at --nmin like the existing tables.

Walter et al. (2012) initialised from MBARI CTD casts near their 1-16 May 2010
observations. C1 sits near the canyon head rather than at Hopkins, but its
18 May 2010 cast is two days after their window and its 9-11 degC range matches
their 15 m thermistors; the profile used until now was fitted to late-summer
conditions.

Also reports the effective thermal expansion coefficient of the cast,
alpha = -(1/rho0) d(rho)/dT across the thermocline, for thermistors.py --alpha.
"""
import argparse
from datetime import datetime, timedelta

import numpy as np
import scipy.io as sio


def dn2dt(x):
    return datetime.fromordinal(int(x)) + timedelta(days=float(x) % 1) - timedelta(days=366)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mat')
    ap.add_argument('--date', default=None, help='YYYY-MM-DD of the cast to use')
    ap.add_argument('--list-month', dest='month', type=int, default=None)
    ap.add_argument('--zmax', type=float, default=110.0, help='depth of the table, m')
    ap.add_argument('--smooth', type=float, default=3.0, help='running mean, m')
    ap.add_argument('--nmin', type=float, default=1e-4, help='floor on N, 1/s')
    ap.add_argument('--rho0', type=float, default=1025.0)
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--fig', default=None, help='optional comparison figure')
    ap.add_argument('--compare', default=None, help='existing N table to plot against')
    g = ap.parse_args()

    d = sio.loadmat(g.mat, squeeze_me=True, struct_as_record=False)
    casts = list(np.atleast_1d(d['CTD'].cast))
    dates = [dn2dt(c.mtime) for c in casts]

    if g.month:
        for c, dt in zip(casts, dates):
            if dt.month == g.month:
                print(f"   {dt:%Y-%m-%d %H:%M}  {int(c.depth.min())}-{int(c.depth.max())} m  "
                      f"T {np.nanmin(c.temp):.2f}-{np.nanmax(c.temp):.2f} degC")
        return
    if not g.date or not g.out:
        raise SystemExit('give --date and -o, or --list-month')

    hits = [i for i, dt in enumerate(dates) if dt.strftime('%Y-%m-%d') == g.date]
    if not hits:
        raise SystemExit(f'no cast on {g.date}')
    c = casts[hits[0]]
    dep = c.depth.astype(float)
    rho = np.asarray(c.dens, float)
    tem = np.asarray(c.temp, float)
    ok = np.isfinite(dep) & np.isfinite(rho) & np.isfinite(tem)
    dep, rho, tem = dep[ok], rho[ok], tem[ok]
    order = np.argsort(dep)
    dep, rho, tem = dep[order], rho[order], tem[order]

    # resample to 0.1 m, extending the shallowest reading to the surface
    zz = np.arange(0.0, g.zmax + 1e-9, 0.1)
    rho_i = np.interp(zz, dep, rho)
    tem_i = np.interp(zz, dep, tem)
    k = max(int(round(g.smooth / 0.1)), 1)
    pad = np.pad(rho_i, (k // 2, k - 1 - k // 2), mode='edge')
    rho_s = np.convolve(pad, np.ones(k) / k, mode='valid')
    n2 = 9.80665 / g.rho0 * np.gradient(rho_s, zz)
    N = np.sqrt(np.clip(n2, g.nmin**2, None))

    with open(g.out, 'w') as fh:
        fh.write(f"# z (m, negative down)   N (1/s)   from {g.mat} cast {g.date}, "
                 f"density smoothed over {g.smooth:g} m\n")
        for z_, n_ in zip(-zz[::-1], N[::-1]):
            fh.write(f"{z_:12.4f}   {n_:.6e}\n")

    # summary and the effective expansion coefficient across the thermocline
    kmax = int(np.argmax(N))
    sel = (zz >= 3) & (zz <= 25)
    alpha = -np.polyfit(tem_i[sel], rho_i[sel], 1)[0] / g.rho0
    print(f"{g.out}: cast {dates[hits[0]]:%Y-%m-%d %H:%M}, "
          f"{int(dep.min())}-{int(dep.max())} m")
    print(f"   N max {N.max():.3e} 1/s at {zz[kmax]:.1f} m depth; "
          f"N at 5/15/30/60 m: " + ", ".join(f"{np.interp(x, zz, N):.2e}" for x in (5, 15, 30, 60)))
    print(f"   effective alpha (3-25 m) = {alpha:.2e} 1/degC  -> thermistors.py --alpha {alpha:.2e}")

    if g.fig:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(8, 5), sharey=True)
        ax[0].plot(tem_i, -zz, label=f'C1 {g.date}')
        ax[0].set_xlabel('T (degC)')
        ax[0].set_ylabel('z (m)')
        ax[1].semilogx(N, -zz, label=f'C1 {g.date}')
        if g.compare:
            p = np.loadtxt(g.compare)
            ax[1].semilogx(p[:, 1], p[:, 0], label=g.compare)
        ax[1].set_xlabel('N (1/s)')
        for a in ax:
            a.grid(alpha=0.3)
            a.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(g.fig, dpi=110)
        print(f"   figure -> {g.fig}")


if __name__ == '__main__':
    main()
