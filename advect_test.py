#!/usr/bin/env python3
"""advect_test.py -- does the scalar advection do what it claims?

  python advect_test.py --mole $env:MOLE_SRC

Three tests, each with an answer known in advance, run on the solver's own
curvilinear grid and mimetic divergence -- not on a toy uniform mesh, because
the point is to exercise the operators the model actually uses.

1. CONSERVATION. Flux form with the mimetic divergence should conserve the
   integral of b to round-off, whatever the scheme and however ugly the field.

2. ORDER OF ACCURACY. A smooth bump translated by a uniform flow has an exact
   solution. upwind1 should converge at first order and upwind2 near second;
   this is the difference that decided the fate of the beams in 2021, so it is
   worth measuring rather than assuming.

3. MONOTONICITY. A step profile must not develop new extrema under minmod.
   upwind2 will overshoot, which is why a limiter exists.
"""
import argparse
import importlib.util
import os
import sys

import numpy as np


def load_solver(path):
    spec = importlib.util.spec_from_file_location('iwbcurv', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['iwbcurv'] = mod
    spec.loader.exec_module(mod)
    return mod


def uniform_grid(Lx, Lz, m, n):
    """(eta, xi) arrays, the layout the solver hands to MOLE.

    The solver pushes the grid with eta as rows on purpose; pushing it the
    other way round swaps which block of the face vector MOLE treats as which
    direction, and the divergence of a constant field then comes out nonzero.
    """
    x = np.linspace(-Lx / 2, Lx / 2, m + 1)
    z = np.linspace(-Lz, 0.0, n + 1)
    X, Z = np.meshgrid(x, z)
    return X, Z


def centres(X, Z):
    n, m = Z.shape[0] - 1, Z.shape[1] - 1
    xc = np.zeros((n + 2, m + 2))
    zc = np.zeros((n + 2, m + 2))
    xc[1:-1, 1:-1] = 0.25 * (X[:-1, :-1] + X[1:, :-1] + X[:-1, 1:] + X[1:, 1:])
    zc[1:-1, 1:-1] = 0.25 * (Z[:-1, :-1] + Z[1:, :-1] + Z[:-1, 1:] + Z[1:, 1:])
    xc[0, :] = xc[1, :]
    xc[-1, :] = xc[-2, :]
    xc[:, 0] = xc[:, 1] - (xc[:, 2] - xc[:, 1])
    xc[:, -1] = xc[:, -2] + (xc[:, -2] - xc[:, -3])
    zc[0, :] = zc[1, :] - (zc[2, :] - zc[1, :])
    zc[-1, :] = zc[-2, :] + (zc[-2, :] - zc[-3, :])
    zc[:, 0] = zc[:, 1]
    zc[:, -1] = zc[:, -2]
    return xc, zc


def run_case(mod, oc_div, X, Z, m, n, U, scheme, periods, nsteps, bump):
    """Translate a field by a uniform flow and return the final state."""
    D = oc_div
    nc = (n + 2) * (m + 2)
    xc, _ = centres(X, Z)
    Lx = X.max() - X.min()
    b = bump(xc)
    u = np.full((n, m + 1), U).ravel()
    w = np.zeros((n + 1) * m)
    dt = periods * Lx / U / nsteps
    for _ in range(nsteps):
        adv = mod.advect_centred(b, u, w, D, n, m, nc, scheme, np)
        B = b.reshape(n + 2, m + 2).copy()
        A = adv.reshape(n + 2, m + 2)
        # ONLY the interior is advanced: the boundary rows of the mimetic
        # divergence are not cell-balance equations, and feeding their output
        # back through the face interpolation diverges within a few steps.
        B[1:-1, 1:-1] -= dt * A[1:-1, 1:-1]
        B[:, 0] = B[:, m]                    # periodic in x
        B[:, m + 1] = B[:, 1]
        B[0, :] = B[1, :]                    # no flux through the eta walls
        B[-1, :] = B[-2, :]
        b = B.ravel()
    return b, xc, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mole', default=os.environ.get('MOLE_SRC'))
    ap.add_argument('--solver', default='iwbcurv.py')
    g = ap.parse_args()
    if not g.mole:
        raise SystemExit('set MOLE_SRC or pass --mole')

    mod = load_solver(g.solver)
    from oct2py import Oct2Py
    oc = Oct2Py()
    oc.addpath(os.path.abspath(g.mole))

    Lx, Lz, U = 1.0, 0.25, 1.0
    print('\n1. CONSERVATION (flux form with the mimetic divergence)')
    for scheme in ('upwind1', 'upwind2', 'minmod'):
        m, n = 64, 16
        X, Z = uniform_grid(Lx, Lz, m, n)
        oc.push('Xm', X.copy())
        oc.push('Zm', Z.copy())
        oc.eval(f"D = div2DCurv(2, Xm, Zm);", verbose=False)
        D = mod.pull(oc, 'D')
        rng = np.random.default_rng(0)
        b0 = rng.standard_normal((n + 2) * (m + 2))
        u = rng.standard_normal((n, m + 1)).ravel()
        w = rng.standard_normal((n + 1) * m).ravel()
        # a constant field carried by any flow: div(u b) must reduce to b div(u),
        # which the projection makes zero -- the sharpest consistency check
        bconst = np.ones((n + 2) * (m + 2))
        adv = mod.advect_centred(bconst, u, w, D, n, m, (n + 2) * (m + 2), scheme, np)
        interior = np.zeros((n + 2, m + 2), bool)
        interior[1:-1, 1:-1] = True
        divu = D @ np.concatenate([u, w])
        err = float(np.abs(adv[interior.ravel()] - divu[interior.ravel()]).max())
        print(f"   {scheme:8s} max|div(u*1) - div(u)| over interior = {err:.2e}")

    print('\n2. ORDER OF ACCURACY (smooth bump, one full traverse)')
    bump = lambda xc: np.exp(-((xc - 0.0) / 0.12) ** 2).ravel()
    for scheme in ('upwind1', 'upwind2', 'minmod'):
        errs, hs = [], []
        for m in (32, 64, 128):
            n = 8
            X, Z = uniform_grid(Lx, Lz, m, n)
            oc.push('Xm', X.copy())
            oc.push('Zm', Z.copy())
            oc.eval("D = div2DCurv(2, Xm, Zm);", verbose=False)
            D = mod.pull(oc, 'D')
            steps = 20 * m
            b, xc, dt = run_case(mod, D, X, Z, m, n, U, scheme, 1.0, steps, bump)
            exact = bump(xc)
            sel = np.zeros((n + 2, m + 2), bool)
            sel[1:-1, 1:-1] = True
            e = np.abs(b - exact)[sel.ravel()].mean()
            errs.append(e)
            hs.append(Lx / m)
        p = [np.log(errs[i] / errs[i + 1]) / np.log(hs[i] / hs[i + 1]) for i in range(2)]
        print(f"   {scheme:8s} errors {errs[0]:.3e} {errs[1]:.3e} {errs[2]:.3e}   "
              f"orders {p[0]:.2f}, {p[1]:.2f}")

    print('\n3. MONOTONICITY (step profile, quarter traverse)')
    step = lambda xc: (xc.ravel() < 0.0).astype(float)
    m, n = 128, 8
    X, Z = uniform_grid(Lx, Lz, m, n)
    oc.push('Xm', X.copy())
    oc.push('Zm', Z.copy())
    oc.eval("D = div2DCurv(2, Xm, Zm);", verbose=False)
    D = mod.pull(oc, 'D')
    for scheme in ('upwind1', 'upwind2', 'minmod'):
        b, xc, dt = run_case(mod, D, X, Z, m, n, U, scheme, 0.25, 5 * m, step)
        sel = np.zeros((n + 2, m + 2), bool)
        sel[1:-1, 1:-1] = True
        v = b[sel.ravel()]
        print(f"   {scheme:8s} range [{v.min():+.4f}, {v.max():+.4f}]   "
              f"overshoot {max(v.max() - 1.0, -v.min(), 0.0):.2e}  "
              + ('OK' if max(v.max() - 1.0, -v.min()) < 1e-9 else 'oscillates'))
    print('\n4. MOMENTUM ADVECTION (exact cases)')
    m, n = 48, 24
    X, Z = uniform_grid(1.0, 1.0, m, n)
    met = mod.logical_metrics(X, Z, np)
    xu = 0.5 * (X[:-1, :] + X[1:, :])
    zu = 0.5 * (Z[:-1, :] + Z[1:, :])
    xw = 0.5 * (X[:, :-1] + X[:, 1:])
    zw = 0.5 * (Z[:, :-1] + Z[:, 1:])

    # (a) uniform flow acting on a uniform field: advection must vanish
    u = np.ones((n, m + 1)).ravel()
    w = 0.5 * np.ones((n + 1) * m)
    au, aw = mod.advect_momentum(u, w, met, n, m, 'upwind2', np)
    print(f"   uniform field, uniform flow: max|adv| = "
          f"{max(np.abs(au).max(), np.abs(aw).max()):.2e}   (exact answer 0)")

    # (b) u = x with flow (1, 0): u du/dx = x, known everywhere
    u = xu.ravel().copy()
    w = np.zeros((n + 1) * m)
    au, aw = mod.advect_momentum(u, w, met, n, m, 'upwind2', np)
    err = np.abs(au.reshape(n, m + 1)[:, 2:-2] - xu[:, 2:-2]).max()
    print(f"   u = x, w = 0, expect u du/dx = x: max error = {err:.2e}")

    # (c) solid-body rotation (u, w) = (-z, x): the advection of each component
    #     is the centripetal term, u.grad(u) = -x and u.grad(w) = -z
    u = (-zu).ravel().copy()
    w = xw.ravel().copy()
    au, aw = mod.advect_momentum(u, w, met, n, m, 'upwind2', np)
    eu = np.abs(au.reshape(n, m + 1)[2:-2, 2:-2] - (-xu[2:-2, 2:-2])).max()
    ew = np.abs(aw.reshape(n + 1, m)[2:-2, 2:-2] - (-zw[2:-2, 2:-2])).max()
    print(f"   solid-body rotation, expect (-x, -z): max error = "
          f"{max(eu, ew):.2e}")

    # (d) the same on a stretched grid, refined: this is where a scheme that
    #     advects with the Cartesian velocity instead of the contravariant one
    #     stops converging. The stretching has to be SMOOTH -- a mapping with a
    #     cusp (|t|^1.6, say) collapses the cell spacing at one point and no
    #     finite-difference scheme converges there.
    print("   solid-body rotation on a smoothly stretched grid:")
    prev = None
    for mm, nn in ((48, 24), (96, 48), (192, 96)):
        tt = np.linspace(-1, 1, mm + 1)
        xs = 0.5 * np.tanh(1.6 * tt) / np.tanh(1.6)
        ss = np.linspace(0, 1, nn + 1)
        zs = -(0.6 * ss + 0.4 * ss ** 2)
        Xs, Zs = np.meshgrid(xs, zs)
        mets = mod.logical_metrics(Xs, Zs, np)
        xus = 0.5 * (Xs[:-1, :] + Xs[1:, :])
        zus = 0.5 * (Zs[:-1, :] + Zs[1:, :])
        xws = 0.5 * (Xs[:, :-1] + Xs[:, 1:])
        zws = 0.5 * (Zs[:, :-1] + Zs[:, 1:])
        au, aw = mod.advect_momentum((-zus).ravel(), xws.ravel(), mets, nn, mm,
                                     'upwind2', np)
        e = max(np.abs(au.reshape(nn, mm + 1)[2:-2, 2:-2] + xus[2:-2, 2:-2]).max(),
                np.abs(aw.reshape(nn + 1, mm)[2:-2, 2:-2] + zws[2:-2, 2:-2]).max())
        order = '' if prev is None else f"   order {np.log(prev / e) / np.log(2):.2f}"
        print(f"     {mm:4d}x{nn:<4d} max error {e:.3e}{order}")
        prev = e
    oc.exit()


if __name__ == '__main__':
    main()
