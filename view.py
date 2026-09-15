#!/usr/bin/env python3
"""
view.py -- plot iwbcurv / ridge2d output on the actual curvilinear grid.

ASCII loses everything that matters here: the grid is curvilinear so cells are
not axis-aligned, the beam is a thin feature at a specific angle, and the
question is always "is the bright stuff ON the theory characteristic or not".

Draws, per field:
  * rms speed on the true (X, Z) cells, so the plot is in physical space
  * the seabed
  * theory characteristics from the ridge crest at +/- asin(w/N), both the
    nonhydrostatic (Eq.30) and hydrostatic (Eq.31) slopes
  * the tracked column maxima and the least-squares fit actually used
  * the fit window, so you can see whether it sits where the beam is

USAGE
  python view.py g0.6.npz
  python view.py g0.2.npz g0.4.npz g0.6.npz g0.8.npz      # panel per file
  python view.py g0.6.npz --grid          # overlay the grid lines
  python view.py g0.6.npz --log           # log colour scale
  python view.py g*.npz -o sweep.png
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Geometry is NOT hardcoded any more.  It used to be
#     D0, AB, LB, LX = 1000.0, 20.0, 30.0, 3000.0
# which silently mis-drew every run that was not the 3 km default: the seabed
# came from the analytic ridge at ab=20 regardless of --ab, and set_xlim
# cropped a 6 or 12 km domain to the middle 3 km.  Everything below is taken
# from the grid stored in the npz, so the bed drawn is the bed used.


def panel(ax, fn, show_grid=False, logscale=False, force_rms=False,
          clim=None, force_total=False, xlim=None):
    d = np.load(fn)
    rms = d['rms']                      # (n+2, m+2) centres incl. ghost ring
    XM, ZM = d['X'], d['Z']             # (n+1, m+1) nodes, (eta, xi)
    ratio = float(d['ratio'])
    pts = d['pts'] if 'pts' in d else None

    snap = ('snap_u' in d) and (not force_rms)
    if snap:
        # instantaneous u / u0 at the u_bc = 0 phase, as Garcia et al. Fig 8
        F = d['snap_u'][1:-1, 1:-1] / float(d['u0'])
        if not force_total:
            # Plot the BAROCLINIC part, u - <u>_z.  The barotropic tide is
            # depth-uniform and carries no beam; here it is ~0.5 u0 at the
            # snapshot phase (the sponge-driven interior lags u_bc, so it is
            # near its peak when u_bc = 0) and swamps the beam on any sensible
            # colour scale.  Removing the depth mean isolates the wave field.
            F = F - F.mean(axis=0, keepdims=True)
        # Scale from the INTERIOR only.  Thin boundary layers at the surface
        # and the bed (the penalty layer) are much stronger than the beam and
        # will otherwise set the range and wash everything out.  Garcia et al.
        # Fig 8 fix the scale at u/u0 = +/- 0.15.
        nb = max(2, F.shape[0] // 12)
        lim = clim if clim else np.percentile(np.abs(F[nb:-nb, :]), 99.0)
        if force_total:
            # The snapshot is taken at the u_bc = 0 phase, but the interior
            # barotropic flow lags the boundary forcing and does NOT pass
            # through zero with it -- measured at -0.27 u0, near-uniform in x.
            # Against a +/-0.15 limit that pins 99-100% of the field and the
            # beam disappears into a flat colour.  Say so rather than draw it.
            bt = float(np.median(d['snap_u'][1:-1, 1:-1].mean(axis=0)
                                 / float(d['u0'])))
            sat = float(np.mean(np.abs(F) > lim))
            if sat > 0.5:
                print(f"  WARNING {os.path.basename(fn)}: --total saturates "
                      f"({100*sat:.0f}% of cells outside +/-{lim:g}).")
                print(f"          depth-mean (barotropic) u/u0 = {bt:+.3f} at "
                      f"this phase, so the whole field is offset.")
                print(f"          Drop --total (baroclinic is the default), or "
                      f"raise --clim above {abs(bt)+lim:.2f}.")
        pc = ax.pcolormesh(XM, ZM, F, cmap='RdBu_r', vmin=-lim, vmax=lim,
                           shading='auto')
        Fq = np.abs(F)
    else:
        F = rms[1:-1, 1:-1]
        Fq = F
    v = (F if snap else F / F.max())
    if (not snap) and logscale:
        v = np.log10(np.maximum(v, 1e-4))
        vmin, vmax = -4, 0
    elif not snap:
        vmin, vmax = 0, 1

    if not snap:
        pc = ax.pcolormesh(XM, ZM, v, cmap='magma', vmin=vmin, vmax=vmax,
                           shading='auto')

    if show_grid:
        for i in range(0, XM.shape[0], 4):
            ax.plot(XM[i, :], ZM[i, :], color='w', lw=0.25, alpha=0.35)
        for j in range(0, XM.shape[1], 4):
            ax.plot(XM[:, j], ZM[:, j], color='w', lw=0.25, alpha=0.35)

    # the actual bed and domain, straight off the grid
    xs, bed = XM[0, :], ZM[0, :]
    LX = float(XM.max() - XM.min())
    D0 = float(-ZM.min())
    ax.fill_between(xs, bed, -D0*1.06, color='#6b5b45', zorder=3)
    ax.plot(xs, bed, color='k', lw=1.4, zorder=4)

    # theory characteristics from the ridge crest
    phi_nh = np.degrees(np.arctan(np.sqrt(ratio**2/(1-ratio**2))))
    phi_h = np.degrees(np.arctan(ratio))
    zc = float(bed[np.argmin(np.abs(xs))])          # true crest height
    for phi, col, ls, lab in ((phi_nh, 'cyan', '-', f'nonhydro {phi_nh:.1f}$\\degree$'),
                              (phi_h, 'lime', '--', f'hydro {phi_h:.1f}$\\degree$')):
        s = np.tan(np.radians(phi))
        for sx in (1, -1):
            xr = np.linspace(0, sx*LX/2, 200)
            zr = zc + np.abs(xr)*s
            k = zr <= 0
            ax.plot(xr[k], zr[k], color=col, ls=ls, lw=1.5, zorder=5,
                    label=lab if sx == 1 else None)

    if pts is not None and len(pts):
        ax.plot(pts[:, 0], pts[:, 1], 'o', ms=2.5, color='white',
                mec='k', mew=0.3, zorder=6, label='tracked maxima')
        ax.plot(-pts[:, 0], pts[:, 1], 'o', ms=2.5, color='white',
                mec='k', mew=0.3, zorder=6)
        A = np.vstack([pts[:, 0], np.ones(len(pts))]).T
        sl, ic = np.linalg.lstsq(A, pts[:, 1], rcond=None)[0]
        meas = np.degrees(np.arctan(abs(sl)))
        xr = np.linspace(pts[:, 0].min(), pts[:, 0].max(), 20)
        for sx in (1, -1):
            ax.plot(sx*xr, sl*xr + ic, color='orange', lw=2.0, zorder=7,
                    label=f'fit {meas:.1f}$\\degree$' if sx == 1 else None)
        for sx in (1, -1):
            ax.axvspan(sx*pts[:, 0].min(), sx*pts[:, 0].max(),
                       color='w', alpha=0.06, zorder=2)
        err = meas - phi_nh
        ttl = (f"$\\omega/N$={ratio}   fit {meas:.2f}$\\degree$ vs "
               f"{phi_nh:.2f}$\\degree$  ({err:+.2f})")
    else:
        ttl = f"$\\omega/N$={ratio}"

    kind = ((f"baroclinic $u'/u_0$ at t={float(d['snap_tT']):.1f}T"
             if not force_total else
             f"total $u/u_0$ at t={float(d['snap_tT']):.1f}T")
            if snap else f"rms, max/mean {Fq.max()/Fq.mean():.1f}")
    geo = (f"Lx={LX/1000:.0f} km  D0={D0:.0f} m  "
           f"{ZM.shape[1]-1}x{ZM.shape[0]-1}  crest {zc+D0:.0f} m")
    ax.set_title(f"{ttl}\n{kind}   {geo}   {os.path.basename(fn)}", fontsize=9)
    if xlim:
        ax.set_xlim(-xlim, xlim)
    else:
        ax.set_xlim(float(XM.min()), float(XM.max()))
    ax.set_ylim(-D0*1.06, 0.02*D0)
    ax.set_xlabel('$x$ (m)'); ax.set_ylabel('$z$ (m)')
    ax.legend(loc='upper right', fontsize=7, framealpha=0.75)
    return pc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('files', nargs='+')
    ap.add_argument('--grid', action='store_true', help='overlay grid lines')
    ap.add_argument('--log', action='store_true', help='log colour scale')
    ap.add_argument('--clim', type=float, default=0.15,
                    help='fixed |u/u0| colour limit (Garcia et al. use 0.15); 0 = auto')
    ap.add_argument('--total', action='store_true',
                    help='plot total u instead of the baroclinic part')
    ap.add_argument('--rms', action='store_true',
                    help='plot the rms average instead of the snapshot')
    ap.add_argument('--xlim', type=float, default=None,
                    help='half-width in m to crop the x axis to, for comparing '
                         'against a figure drawn on a different domain (Garcia '
                         'et al. Fig 8 is +/-1500 m).  Cropping only; the run '
                         'itself is unchanged.')
    ap.add_argument('-o', '--out', default=None)
    g = ap.parse_args()

    nf = len(g.files)
    ncol = 1 if nf == 1 else 2
    nrow = int(np.ceil(nf / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(9*ncol, 3.6*nrow),
                             squeeze=False, constrained_layout=True)
    pc = None
    for ax, fn in zip(axes.ravel(), g.files):
        pc = panel(ax, fn, g.grid, g.log, g.rms, g.clim or None, g.total, g.xlim)
    for ax in axes.ravel()[nf:]:
        ax.axis('off')
    lab = ('rms speed / max' if g.rms else '$u/u_0$')
    fig.colorbar(pc, ax=axes.ravel().tolist(), shrink=0.7, label=lab)
    fig.suptitle('Internal wave beam over a Gaussian ridge, fully curvilinear grid\n'
                 'cyan = nonhydrostatic theory, green = hydrostatic, '
                 'orange = least-squares fit to the white tracked maxima',
                 fontsize=11)
    out = g.out or (os.path.splitext(g.files[0])[0] + '.png' if nf == 1
                    else 'iwb_fields.png')
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")


if __name__ == '__main__':
    main()
