#!/usr/bin/env python3
"""animate.py -- turn saved snapshots into an animation.

  python iwbcurv.py ... --frames frames_08 --framerate 12
  python animate.py frames_08 -o beam.gif
  python animate.py frames_08 -o beam.mp4 --field w --fps 20

Shows the field evolving on the real curvilinear grid, from the first step to
whenever the run stopped, so the beams can be watched forming out of the
start-up transient.

  --field speed   sqrt(u^2 + w^2), the default: shows the beams most clearly
  --field w       vertical velocity, signed: shows the phase propagating
  --field u       horizontal velocity, signed: dominated by the barotropic tide

MP4 needs ffmpeg on PATH; GIF needs nothing beyond pillow. If ffmpeg is
missing, ask for a .gif.

The colour scale is fixed across the animation, computed from a late frame, so
brightening means the field is actually growing rather than the scale moving.
With --clim-from-frame the reference frame can be chosen; with --log the scale
is logarithmic, which suits a beam that spans orders of magnitude.
"""
import argparse
import glob
import os

import numpy as np


def load_frames(d):
    grid = np.load(os.path.join(d, 'grid.npz'))
    files = sorted(glob.glob(os.path.join(d, 'frame_*.npz')))
    if not files:
        raise SystemExit(f'no frames in {d}/')
    return grid, files


def field_of(fr, which, shape, baroclinic=True):
    """Field at cell centres, reshaped to (nz+2, nx+2).

    The barotropic tide fills the whole domain at a few times 1e-3 m/s and
    swamps the beams, so the depth mean of u is removed by default: what is
    left is the baroclinic signal the beams are made of. --with-barotropic
    keeps it.
    """
    u = fr['u'].reshape(shape)
    w = fr['w'].reshape(shape)
    if baroclinic and which != 'b':
        u = u - u.mean(axis=0, keepdims=True)         # rows are depth; x is fastest
    if which == 'b':
        # buoyancy: the field that makes a lock release or a bore visible at all
        return fr['b'].reshape(shape) if 'b' in fr else np.zeros(shape)
    if which == 'speed':
        return np.sqrt(u**2 + w**2)
    if which == 'w':
        return w
    return u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('-o', '--out', default='animation.gif')
    ap.add_argument('--field', default='speed', choices=['speed', 'w', 'u', 'b'])
    ap.add_argument('--fps', type=float, default=12.0)
    ap.add_argument('--log', action='store_true',
                    help='logarithmic colour scale (speed only)')
    ap.add_argument('--clim', type=float, default=None,
                    help='fix the colour limit by hand instead of from a frame')
    ap.add_argument('--clim-from-frame', dest='climfrom', type=float, default=0.9,
                    help='fraction through the run of the frame that sets the scale')
    ap.add_argument('--dpi', type=int, default=110)
    ap.add_argument('--stride', type=int, default=1, help='use every Nth frame')
    ap.add_argument('--with-barotropic', dest='baro', action='store_true',
                    help='keep the depth-mean flow instead of removing it')
    ap.add_argument('--xlim', type=float, default=None)
    g = ap.parse_args()

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import animation, colors

    grid, files = load_frames(g.dir)
    files = files[::max(g.stride, 1)]
    X = grid['X']
    Z = grid['Z']
    D0 = float(grid['D0'])
    n, m = Z.shape[0] - 1, Z.shape[1] - 1

    shape = (n + 2, m + 2)
    ref = np.load(files[min(int(g.climfrom*len(files)), len(files) - 1)])
    refv = field_of(ref, g.field, shape, not g.baro)[1:-1, 1:-1]
    scale = g.clim if g.clim else float(np.percentile(np.abs(refv), 99.8))
    signed = g.field in ('w', 'u', 'b')

    fig, ax = plt.subplots(figsize=(11, 4.2))
    first = field_of(np.load(files[0]), g.field, shape, not g.baro)[1:-1, 1:-1]
    if signed:
        norm = colors.Normalize(-scale, scale)
        cmap = 'RdBu_r'
    elif g.log:
        norm = colors.LogNorm(max(scale*1e-3, 1e-12), scale)
        cmap = 'magma'
    else:
        norm = colors.Normalize(0, scale)
        cmap = 'magma'
    mesh = ax.pcolormesh(X, Z, first, norm=norm, cmap=cmap, shading='auto')

    # the bed, from the grid itself
    br = 0 if Z[0, :].mean() < Z[-1, :].mean() else -1
    ax.fill_between(X[br, :], Z[br, :], Z[br, :].min() - 0.05*D0,
                    color='saddlebrown', alpha=0.85, zorder=3)
    ax.set_xlabel('x (m)')
    ax.set_ylabel('z (m)')
    if g.xlim:
        ax.set_xlim(-g.xlim, g.xlim)
    fig.colorbar(mesh, ax=ax,
                 label=('buoyancy (m/s^2)' if g.field == 'b' else f'{g.field} (m/s)'))
    title = ax.set_title('')
    fig.tight_layout()

    def draw(i):
        fr = np.load(files[i])
        v = field_of(fr, g.field, shape, not g.baro)[1:-1, 1:-1]
        mesh.set_array(v.ravel())
        _unit = 'm/s^2' if g.field == 'b' else 'm/s'
        _pre = '' if (g.baro or g.field == 'b') else 'baroclinic '
        title.set_text(f"t / T = {float(fr['tT']):6.2f}      {_pre}{g.field}, "
                       f"scale fixed at {scale:.2e} {_unit}")
        return mesh, title

    anim = animation.FuncAnimation(fig, draw, frames=len(files), blit=False)
    if g.out.lower().endswith('.gif'):
        anim.save(g.out, writer=animation.PillowWriter(fps=g.fps), dpi=g.dpi)
    else:
        try:
            anim.save(g.out, writer=animation.FFMpegWriter(fps=g.fps, bitrate=2400),
                      dpi=g.dpi)
        except (FileNotFoundError, RuntimeError) as e:
            raise SystemExit(f"could not write {g.out} ({e.__class__.__name__}); "
                             f"ffmpeg is needed for video -- ask for a .gif instead")
    plt.close(fig)
    print(f"{g.out}: {len(files)} frames, "
          f"t/T {float(np.load(files[0])['tT']):.2f} to "
          f"{float(np.load(files[-1])['tT']):.2f}, {g.fps:g} fps")


if __name__ == '__main__':
    main()
