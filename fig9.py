#!/usr/bin/env python3
"""fig9.py -- reproduce Garcia et al. (2019) Fig 9 from the run logs.

  python fig9.py paper-logs [-o fig9_repro.png]

Fig 9 plots the measured internal wave beam angle against omega/N, with the
hydrostatic (Eq. 31) and nonhydrostatic (Eq. 30) curves.  This draws the same
axes, our points, and GCCOM's published points read off their figure.

WHICH OF OUR RUNS IS PLOTTED, and why it matters.  Sections 8 and 11 showed
the measured angle depends strongly on how far OUT the fit window reaches --
2.12 deg of spread at w/N=0.6 across windows ending at 0.30 to 0.75 of a
bounce -- and only weakly on where it starts (0.28 deg).  So a single point
per frequency is only meaningful once the window convention is fixed.  This
selects, per frequency, the usable fit with the LONGEST window at the finest
grid, which in practice means the one reaching ~0.75 of a bounce.  The window
actually used is printed and annotated, so the convention travels with the
figure instead of being hidden in it.

Pass --win-substr to force a particular window instead.
"""
import argparse, os, sys, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_paper import scan, GCCOM, theory


def pick(rows, ratio, win_substr=None):
    """Usable fit, longest window, finest grid, for one frequency."""
    c = [r for r in rows if round(r['ratio'], 1) == ratio
         and r['rms'] < 10 and r['mm'] > 5]
    if win_substr:
        c = [r for r in c if win_substr in r['win']]
    if not c:
        return None
    def span(r):
        a, b = (float(v) for v in r['win'].split(','))
        return b - a
    best_grid = max(r['nx']*r['nz'] for r in c)
    c = [r for r in c if r['nx']*r['nz'] == best_grid]
    return max(c, key=span)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('logdir', nargs='?', default='paper-logs')
    ap.add_argument('-o', '--out', default='fig9_repro.png')
    ap.add_argument('--win-substr', default=None)
    g = ap.parse_args()

    rows = [r for r in (scan(f) for f in
                        sorted(glob.glob(os.path.join(g.logdir, '*.log')))) if r]
    if not rows:
        raise SystemExit(f"no completed runs in {g.logdir}/")

    ratios = [0.2, 0.4, 0.6, 0.8]
    sel = {r: pick(rows, r, g.win_substr) for r in ratios}
    missing = [r for r in ratios if sel[r] is None]
    if missing:
        print(f"warning: no usable run for omega/N = {missing}")

    print("plotted points:")
    for r in ratios:
        s = sel[r]
        if not s:
            continue
        print(f"   w/N={r}  {s['meas']:.2f} deg  (theory {theory(r):.2f}, "
              f"{s['meas']-theory(r):+.2f})  Lx={s['Lx']:.0f}km "
              f"{s['nx']}x{s['nz']} win[{s['win']}] rms {s['rms']:.1f} m")

    x = np.linspace(0.001, 0.995, 400)
    fig, ax = plt.subplots(figsize=(6.2, 5.4), constrained_layout=True)
    ax.plot(x, np.degrees(np.arcsin(x)), 'k-', lw=1.6,
            label='nonhydrostatic, Eq. 30')
    ax.plot(x, np.degrees(np.arctan(x)), 'k--', lw=1.4,
            label='hydrostatic, Eq. 31')
    ax.plot(list(GCCOM.keys()), list(GCCOM.values()), 'ks', ms=9,
            mfc='none', mew=1.6, ls='none',
            label='GCCOM, Garcia et al. Fig 9 (digitized)')
    px = [r for r in ratios if sel[r]]
    py = [sel[r]['meas'] for r in px]
    ax.plot(px, py, 'o', ms=9, color='#c0392b', ls='none',
            label='this solver, fully curvilinear')

    ax.set_xlim(0, 1); ax.set_ylim(0, 90)
    ax.set_xlabel(r'$\omega/N$'); ax.set_ylabel(r'Beam Angle  $\phi$  (deg)')
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 91, 10))
    ax.grid(alpha=0.25)
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)

    wins = {r: sel[r]['win'] for r in px}
    note = "fit windows (m): " + ", ".join(f"{r}:[{wins[r]}]" for r in px)
    mb_o = np.mean([abs(sel[r]['meas']-theory(r)) for r in px])
    mb_g = np.mean([abs(GCCOM[r]-theory(r)) for r in px])
    ax.set_title("Internal wave beam angle vs forcing frequency\n"
                 f"mean |bias|: this solver {mb_o:.2f}$\\degree$, "
                 f"GCCOM {mb_g:.2f}$\\degree$", fontsize=10)
    ax.text(0.02, -0.13, note, transform=ax.transAxes, fontsize=6.5,
            va='top', color='0.35')

    fig.savefig(g.out, dpi=160, bbox_inches='tight')
    print(f"\nwrote {g.out}")


if __name__ == '__main__':
    main()
