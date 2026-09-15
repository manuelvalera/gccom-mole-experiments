#!/usr/bin/env python3
"""bedorder.py -- what order is the bed constraint operator?

  python bedorder.py --mole $env:MOLE_SRC

The beam angle converges at ~first order (fitted p = 1.0-1.3 on a four-grid
series) even though the interior scheme is second order and `--order 2/4/6`
changes nothing.  This tests the bed constraint in isolation, with no
timestepping.

Construct an analytic flow that is EXACTLY tangent to the bed:

    psi(x,z) = (z - b(x)) * cos(k x) * exp(z/H)      psi = 0 on z = b(x)
    u =  d(psi)/dz        w = -d(psi)/dx

psi vanishes along the bed, so the velocity there is exactly tangential and
the true constraint residual is zero.  Sample u and w at the face centres the
solver uses, apply C, and see how the residual falls with h.

    residual ~ h^2  -> the constraint is second order; look elsewhere
    residual ~ h^1  -> the constraint is first order and is the likely cause

SUSPECT, stated before running: the u in C's `avg(u)` sits at iu2[0,:], which
is half a cell ABOVE the bed, while w sits on the bed face itself.  Equating
them across that offset is an O(dz) error.  If that is right, the residual
scales like h and shrinks when the two are evaluated at the same height.
"""
import argparse, os
import numpy as np


def build(mole, grids, nx, nz, Lx, D0, ab, lb, bt, bb):
    from oct2py import Oct2Py
    gd = os.path.join(grids, 'iwbridge')
    with open(os.path.join(gd, 'geom_over.m'), 'w') as fh:
        fh.write("function [Lx, D0, ab, Lb] = geom_over()\n"
                 f"    Lx = {Lx}; D0 = {D0}; ab = {ab}; Lb = {lb};\nend\n")
    oc = Oct2Py(); oc.eval("warning('off','all'); more off;")
    oc.addpath(mole); oc.addpath(os.path.abspath(gd))
    here = os.getcwd(); os.chdir(grids)
    try:
        oc.eval(f"global BT BB; BT={bt}; BB={bb};")
        oc.eval(f"cd('{os.path.abspath('.').replace(os.sep,'/')}'); "
                f"[X,Z]=gridGen('TFI','iwbridge',{nx},{nz},false);")
    finally:
        os.chdir(here)
    X = oc.pull('X'); Z = oc.pull('Z'); oc.exit()
    return X.T.copy(), Z.T.copy()


def bed(x, D0, ab, lb):
    return -(D0 - ab*np.exp(-x**2/(2*lb**2)))


def dbed(x, D0, ab, lb):
    return ab*np.exp(-x**2/(2*lb**2))*(-x/lb**2)


def flow(x, z, D0, ab, lb, k, H):
    """u, w from psi = (z - b(x)) cos(kx) exp(z/H); psi = 0 on the bed."""
    b = bed(x, D0, ab, lb); bp = dbed(x, D0, ab, lb)
    c, s, e = np.cos(k*x), np.sin(k*x), np.exp(z/H)
    # psi   = (z-b) c e
    # dpsi/dz = c e + (z-b) c e / H
    u = c*e + (z-b)*c*e/H
    # dpsi/dx = -bp c e + (z-b)(-k s) e
    w = -(-bp*c*e + (z-b)*(-k*s)*e)
    return u, w


def residual(XM, ZM, D0, ab, lb, k, H, shift_u_to_bed=False):
    n, m = ZM.shape[0]-1, ZM.shape[1]-1
    Xf = 0.5*(XM[:, :-1] + XM[:, 1:]); Zf = 0.5*(ZM[:, :-1] + ZM[:, 1:])
    Xu = 0.5*(XM[:-1, :] + XM[1:, :]); Zu = 0.5*(ZM[:-1, :] + ZM[1:, :])

    zx_ = np.gradient(ZM[0, :]); xx_ = np.gradient(XM[0, :])
    zc_ = 0.5*(zx_[:-1] + zx_[1:]); xc_ = 0.5*(xx_[:-1] + xx_[1:])
    rat = zc_/np.where(np.abs(xc_) < 1e-12, 1e-12, xc_)

    xw, zw = Xf[0, :], Zf[0, :]                  # bed w-faces, ON the bed
    xu, zu = Xu[0, :], Zu[0, :]                  # first u row, HALF A CELL UP
    if shift_u_to_bed:
        zu = np.interp(xu, XM[0, :], ZM[0, :])   # same height as w

    _, w = flow(xw, zw, D0, ab, lb, k, H)
    u, _ = flow(xu, zu, D0, ab, lb, k, H)
    r = w - rat*0.5*(u[:-1] + u[1:])
    scale = np.abs(flow(xw, zw, D0, ab, lb, k, H)[1]).max()
    off = float(np.mean(np.abs(Zu[0, :] - np.interp(Xu[0, :], XM[0, :], ZM[0, :]))))
    return float(np.abs(r).max()), scale, off


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--mole', default=os.environ.get('MOLE_SRC'))
    p.add_argument('--grids', default='grids')
    p.add_argument('--Lx', type=float, default=6000.0)
    p.add_argument('--D0', type=float, default=1000.0)
    p.add_argument('--ab', type=float, default=20.0)
    p.add_argument('--lb', type=float, default=30.0)
    p.add_argument('--bt', type=float, default=1.6)
    p.add_argument('--bb', type=float, default=3.2)
    g = p.parse_args()
    if not g.mole or not os.path.isfile(os.path.join(g.mole, 'gridGen.m')):
        raise SystemExit("set --mole / MOLE_SRC")

    k, H = 2*np.pi/1500.0, 400.0
    grids = [(256, 101), (384, 151), (512, 201), (640, 251)]

    for shift, label in ((False, "AS BUILT: u half a cell above the bed"),
                         (True,  "CONTROL: u sampled ON the bed")):
        print(f"\n{label}")
        print("%-12s %-10s %-14s %-14s %-8s" %
              ("grid", "dx (m)", "u-offset (m)", "max residual", "order"))
        prev = None
        for nx, nz in grids:
            XM, ZM = build(g.mole, g.grids, nx, nz, g.Lx, g.D0, g.ab, g.lb,
                           g.bt, g.bb)
            r, sc, off = residual(XM, ZM, g.D0, g.ab, g.lb, k, H, shift)
            h = g.Lx/nx
            o = ""
            if prev:
                o = "%.2f" % (np.log(r/prev[1])/np.log(h/prev[0]))
            print("%-12s %-10.2f %-14.3f %-14.3e %-8s"
                  % (f"{nx}x{nz}", h, off, r/sc, o))
            prev = (h, r)


if __name__ == '__main__':
    main()
