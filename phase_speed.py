#!/usr/bin/env python3
"""phase_speed.py -- does the imposed wave travel at the speed it should?

  python iwbcurv.py ... --mode1 3.0 --frames fr_mode --framerate 12
  python phase_speed.py fr_mode --omega 1.405e-4

The structure test shows the wave is generated with the right vertical shape.
This checks the other half: that it propagates at the mode speed. Two
independent measurements, neither using the value the solver computed.

  WAVELENGTH. At the depth where the mode peaks, the horizontal profile of w is
  a sinusoid in x. Its dominant Fourier component gives k, hence lambda = 2 pi / k
  and c = omega / k.

  PHASE TRACKING. The phase of that component advances linearly in time; its
  rate is omega, and dividing by k gives c again -- this one is sensitive to the
  direction of travel, which the wavelength alone is not.

Both are taken in the interior, clear of the sponges where the wave is imposed
and absorbed.
"""
import argparse
import glob
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('--omega', type=float, required=True, help='forcing frequency, rad/s')
    ap.add_argument('--skip', type=float, default=0.35,
                    help='fraction of the domain to drop at each end (sponges)')
    ap.add_argument('--from-frame', dest='ff', type=float, default=0.4,
                    help='fraction of the run to skip before measuring')
    g = ap.parse_args()

    grid = np.load(os.path.join(g.dir, 'grid.npz'))
    files = sorted(glob.glob(os.path.join(g.dir, 'frame_*.npz')))
    if len(files) < 4:
        raise SystemExit(f'need at least 4 frames in {g.dir}/')
    files = files[int(g.ff * len(files)):]

    X, Z = grid['X'], grid['Z']
    n, m = Z.shape[0] - 1, Z.shape[1] - 1
    xc = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:] + X[1:, 1:])
    zc = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])

    # depth where the wave is strongest, from the first frame used
    w0 = np.load(files[0])['w'].reshape(n + 2, m + 2)[1:-1, 1:-1]
    row = int(np.argmax(np.abs(w0).mean(axis=1)))
    xs = xc[row, :]
    keep = (xs > xs.min() + g.skip * np.ptp(xs)) & (xs < xs.max() - g.skip * np.ptp(xs))
    xs = xs[keep]
    dx = float(np.mean(np.diff(xs)))
    win = np.hanning(len(xs))

    # k from the gradient of the analytic-signal phase, not from FFT bins: with a
    # window only a wavelength or two wide the dominant bin is just the window
    # length, which is how this first read 5.957 km for a 6.83 km wave.
    from scipy.signal import hilbert
    ks, phases, times, amps = [], [], [], []
    for f in files:
        fr = np.load(f)
        w = fr['w'].reshape(n + 2, m + 2)[1:-1, 1:-1][row, :][keep]
        a = hilbert(w - w.mean())
        phx = np.unwrap(np.angle(a))
        mid = slice(len(xs) // 6, -len(xs) // 6 or None)
        ks.append(abs(np.polyfit(xs[mid], phx[mid], 1)[0]))
        phases.append(np.angle(a).mean() if False else float(phx[len(xs) // 2]))
        amps.append(float(np.abs(a).mean()))
        times.append(float(fr['t']))
    k = float(np.median(ks))
    lam = 2 * np.pi / k
    times = np.array(times)
    ph = np.unwrap(np.array(phases))
    rate = np.polyfit(times, ph, 1)[0]

    print(f"\nmeasured at z = {zc[row, 0]:.1f} m, over x = {xs.min():.0f} to {xs.max():.0f} m")
    print(f"  wavelength   {lam/1000:.3f} km   (k = {k:.3e} 1/m, "
          f"spread across frames {100*np.std(ks)/k:.1f}%)")
    print(f"  phase speed from wavelength   c = omega/k = {g.omega/k:.4f} m/s")
    print(f"  phase advance {rate:+.4e} rad/s vs omega {g.omega:.4e} "
          f"({abs(rate)/g.omega:.3f} of it, sign = direction of travel)")
    print(f"  phase speed from tracking     c = {abs(rate)/k:.4f} m/s")
    print(f"  amplitude drift over the window: "
          f"{100*(amps[-1]/amps[0] - 1):+.1f}%  (0 = no growth or decay)")


if __name__ == '__main__':
    main()
