#!/usr/bin/env python3
"""
cone3d.py -- 3-D internal-wave cone test against MOLE's operators.

A point source oscillating at omega in uniform-N fluid radiates a CONE at
theta = asin(omega/N) from the horizontal.  The 2-D St Andrew's cross is a
slice through it.  Exact answer, no fitted parameters.

  cart   Cartesian grid, grad3D / div3D            -> validates the solver
  curv   curvilinear grid, grad3DCurv / div3DCurv  -> exercises GI13

--patch puts the fixed GI13 ahead of MOLE's on the Octave path.  Running
`curv` both ways is the direct test of the patch on physics.

SETUP
    pip install oct2py numpy scipy
    # Windows: point oct2py at Octave if it is not on PATH
    set OCTAVE_EXECUTABLE=C:\\Users\\<you>\\AppData\\Local\\Programs\\GNU Octave\\Octave-11.3.0\\mingw64\\bin\\octave-cli.exe

USAGE
    python cone3d.py cart --mole C:/path/to/mole/src/matlab_octave --a 48
    python cone3d.py curv --mole ... --a 48
    python cone3d.py curv --mole ... --a 48 --patch

Everything needed is printed to stdout.  Paste the whole block back.
"""
import argparse, os, sys, tempfile, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from oct2py import Oct2Py

PATCHED_GI13 = r"""
function I = GI13(M, m, n, o, type)
  switch type
    case 'Gn',  Ax=P(m); Ay=Q(n); Az=speye(o);
    case 'Gc',  Ax=P(m); Ay=speye(n); Az=Q(o);
    case 'Ge',  Ax=Q(m); Ay=P(n); Az=speye(o);
    case 'Gcy', Ax=speye(m); Ay=P(n); Az=Q(o);
    case 'Gee', Ax=Q(m); Ay=speye(n); Az=P(o);
    case 'Gnn', Ax=speye(m); Ay=Q(n); Az=P(o);
    otherwise, error('GI13: unknown type %s', type);
  end
  I = kron(Az, kron(Ay, Ax)) * M;
end
function A = Q(N)
  A = spdiags(0.5*ones(N,2), [0 1], N, N+1);
end
function A = P(N)
  A = spdiags(0.5*ones(N+1,2), [-1 0], N+1, N);
  A(1,1)=1.5; A(1,2)=-0.5; A(N+1,N)=1.5; A(N+1,N-1)=-0.5;
end
"""


def pull(oc, nm):
    oc.eval(f"[ii,jj,vv]=find({nm}); sz=size({nm});")
    return sp.csr_matrix(
        (oc.pull('vv').ravel(),
         (oc.pull('ii').ravel().astype(int) - 1,
          oc.pull('jj').ravel().astype(int) - 1)),
        shape=tuple(oc.pull('sz').ravel().astype(int)))


def build(mole, mode, a, amp, patched):
    oc = Oct2Py(); oc.eval("more off; warning('off','all');")
    oc.addpath(mole)
    if patched:
        d = os.path.join(tempfile.gettempdir(), 'cone3d_patch')
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'GI13.m'), 'w').write(PATCHED_GI13)
        oc.addpath(d)
    which = oc.feval("which", "GI13").strip()

    N = a + 1
    oc.eval(f"k=2; a={a}; N={N}; amp={amp}; h=1/a;")
    if mode == 'cart':
        oc.eval("G=grad3D(k,a,h,a,h,a,h); D=div3D(k,a,h,a,h,a,h);")
    else:
        oc.eval("""
        [p,q,r]=meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
        X=p+amp*sin(2*pi*q); Y=q+amp*sin(2*pi*p);
        Z=r+0.5*amp*sin(2*pi*p).*sin(2*pi*q);
        G=grad3DCurv(k,X,Y,Z); D=div3DCurv(k,X,Y,Z);
        """)
    oc.eval("B=robinBC3D(k,a,h,a,h,a,h,0,1);")
    oc.eval("Ic=interpol3D(a,a,a,0.5,0.5,0.5); Idf=interpolD3D(a,a,a,0.5,0.5,0.5);")
    ops = [pull(oc, s) for s in ('G', 'D', 'B', 'Ic', 'Idf')]

    s = np.linspace(0, 1, N); ctr = s[:-1] + 0.5 / a
    c = np.concatenate(([0], ctr, [1]))
    if mode == 'cart':
        Pw, Qw, Rw = np.meshgrid(ctr, ctr, s, indexing='ij')
        Xw, Yw, Zw = Pw, Qw, Rw
        Pc, Qc, Rc = np.meshgrid(c, c, c, indexing='ij')
        Xc, Yc, Zc = Pc, Qc, Rc
    else:
        Pw, Qw, Rw = np.meshgrid(ctr, ctr, s, indexing='ij')
        Xw = Pw + amp * np.sin(2 * np.pi * Qw)
        Yw = Qw + amp * np.sin(2 * np.pi * Pw)
        Zw = Rw + 0.5 * amp * np.sin(2 * np.pi * Pw) * np.sin(2 * np.pi * Qw)
        Pc, Qc, Rc = np.meshgrid(c, c, c, indexing='ij')
        Xc = Pc + amp * np.sin(2 * np.pi * Qc)
        Yc = Qc + amp * np.sin(2 * np.pi * Pc)
        Zc = Rc + 0.5 * amp * np.sin(2 * np.pi * Pc) * np.sin(2 * np.pi * Qc)
    fl = lambda A: A.transpose(2, 1, 0).ravel()      # x fastest
    oc.exit()
    return ops, (fl(Xw), fl(Yw), fl(Zw)), (fl(Xc), fl(Yc), fl(Zc)), which


def run(mole, mode, a, amp, ratio, nper, patched, sig, spg):
    (G, D, B, Ic, Idf), Wp, Cp, which = build(mole, mode, a, amp, patched)
    m = n = o = a
    nu, nv, nw = (m + 1) * n * o, m * (n + 1) * o, m * n * (o + 1)
    nc = (m + 2) ** 3
    print(f"  GI13 : {which}")
    print(f"  grid : {a}^3   faces {nu+nv+nw}   centres {nc}")

    L = (D @ G + B)
    asym = abs(L - L.T).max() / abs(L).max()
    L = L.tolil(); ip = nc // 2; L[ip, :] = 0; L[ip, ip] = 1; L = L.tocsc()
    print(f"  L    : nonsymmetry {asym:.2e}  (D*G+BC is not SPD -- AMG will not work)")
    t = time.time(); lu = spl.splu(L)
    print(f"  splu : {time.time()-t:.1f}s"); sys.stdout.flush()

    Nb = 1.0; om = ratio * Nb
    xw, yw, zw = Wp; xc, yc, zc = Cp
    Fw = 1e-3 * np.exp(-(((xw-.5)**2 + (yw-.5)**2 + (zw-.5)**2) / sig**2))
    dW = np.minimum.reduce([xw, 1-xw, yw, 1-yw, zw, 1-zw])
    dC = np.minimum.reduce([xc, 1-xc, yc, 1-yc, zc, 1-zc])
    Sw = (np.maximum(0, spg - dW) / spg)**2
    Sc = (np.maximum(0, spg - dC) / spg)**2
    tau = 0.5 / Nb

    T = 2*np.pi/om; dt = min(0.02/Nb, T/300); nt = int(nper*T/dt)
    u = np.zeros(nu); v = np.zeros(nv); w = np.zeros(nw); b = np.zeros(nc)
    Wc = np.zeros(nc, complex)
    iu = np.arange(nu).reshape(o, n, m+1); um = np.ones(nu)
    um[iu[:, :, 0].ravel()] = 0; um[iu[:, :, -1].ravel()] = 0
    iv = np.arange(nv).reshape(o, n+1, m); vm = np.ones(nv)
    vm[iv[:, 0, :].ravel()] = 0; vm[iv[:, -1, :].ravel()] = 0
    iw = np.arange(nw).reshape(o+1, n, m); wm = np.ones(nw)
    wm[iw[0].ravel()] = 0; wm[iw[-1].ravel()] = 0
    z0 = np.zeros(nu+nv)

    t0 = time.time()
    for it in range(1, nt+1):
        tt = it*dt
        bz = (Ic @ b)[nu+nv:]
        us = u*um; vs = v*vm
        ws = (w + dt*(bz + Fw*np.sin(om*tt)))*wm
        rhs = (D @ np.concatenate([us, vs, ws]))/dt; rhs[ip] = 0
        p = lu.solve(rhs); gp = G @ p
        u = (us - dt*gp[:nu])*um
        v = (vs - dt*gp[nu:nu+nv])*vm
        w = (ws - dt*gp[nu+nv:])*wm
        w -= dt*Sw*w/tau
        wc = Idf @ np.concatenate([z0, w])
        b = b - dt*Nb**2*wc - dt*Sc*b/tau
        if tt/T > nper-2:
            Wc += wc*np.exp(-1j*om*tt)*dt
        if it % max(1, nt//5) == 0:
            print(f"    t/T={tt/T:5.2f}  max|w|={np.abs(w).max():.3e}  "
                  f"({time.time()-t0:.0f}s)"); sys.stdout.flush()
    return np.abs(Wc), Cp, a


def report(A, C, a, ratio):
    nc = a + 2
    xc, yc, zc = [q.reshape(nc, nc, nc) for q in C]
    F = A.reshape(nc, nc, nc)
    j = nc // 2
    S, xs, zs = F[:, j, :], xc[:, j, :], zc[:, j, :]
    tgt = np.degrees(np.arcsin(ratio))

    print(f"\n  x-z slice through the source  (theory {tgt:.1f} deg)")
    step = max(1, nc // 22)
    Cm = S[::step, ::step] / S.max(); ch = ' .:-=+*#%@'
    for row in Cm[::-1]:
        print('    ' + ''.join(ch[min(9, int(v*9.999))] for v in row))

    R = np.hypot(xs-.5, zs-.5)
    TH = np.degrees(np.arctan2(np.abs(zs-.5), np.abs(xs-.5)))
    sel = (R > 0.18) & (R < 0.36)
    prof = np.array([S[sel & (np.abs(TH-q) < 3)].mean()
                     if (sel & (np.abs(TH-q) < 3)).sum() else np.nan
                     for q in range(0, 91, 5)])
    mx = np.nanmax(prof)
    print(f"\n  angular profile, deg from horizontal   (peak should sit at {tgt:.1f})")
    for q, vv in zip(range(0, 91, 5), prof):
        bar = '#'*int(45*vv/mx) if np.isfinite(vv) else ''
        mark = '  <== theory' if abs(q-tgt) < 2.5 else ''
        print(f"    {q:3d} {bar}{mark}")

    # 3-D angular peak, excluding the forced vertical column
    R3 = np.sqrt((xc-.5)**2 + (yc-.5)**2 + (zc-.5)**2)
    hz = np.hypot(xc-.5, yc-.5)
    A3 = np.degrees(np.arctan2(np.abs(zc-.5), hz))
    est = []
    for rr in np.linspace(0.20, 0.36, 10):
        s = np.abs(R3-rr) < 0.012
        if s.sum() < 80: continue
        aa, vv = A3[s], F[s]
        k = (aa > 6) & (aa < 80)
        if k.sum() < 40: continue
        g = np.arange(6, 80, 0.5)
        sm = np.array([vv[k][np.abs(aa[k]-q) < 2.5].mean()
                       if (np.abs(aa[k]-q) < 2.5).sum() else np.nan for q in g])
        if not np.all(np.isnan(sm)): est.append(g[np.nanargmax(sm)])
    meas = float(np.median(est)) if est else float('nan')
    print(f"\n  CONE HALF-ANGLE  measured {meas:.2f}   theory {tgt:.2f}   "
          f"error {meas-tgt:+.2f} deg ({100*(meas-tgt)/tgt:+.1f}%)")
    return meas


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['cart', 'curv'])
    ap.add_argument('--mole', required=True)
    ap.add_argument('--a', type=int, default=48)
    ap.add_argument('--ratio', type=float, default=0.6)
    ap.add_argument('--amp', type=float, default=0.10)
    ap.add_argument('--nper', type=int, default=10)
    ap.add_argument('--sig', type=float, default=0.10)
    ap.add_argument('--spg', type=float, default=0.09)
    ap.add_argument('--patch', action='store_true')
    ap.add_argument('--out', default=None)
    g = ap.parse_args()
    print(f"\n=== cone3d  {g.mode}  a={g.a}  omega/N={g.ratio}  amp={g.amp}  "
          f"nper={g.nper}  {'PATCHED' if g.patch else 'STOCK'} GI13 ===")
    A, C, a = run(g.mole, g.mode, g.a, g.amp, g.ratio, g.nper, g.patch, g.sig, g.spg)
    report(A, C, a, g.ratio)
    if g.out:
        np.savez_compressed(g.out, A=A, xc=C[0], yc=C[1], zc=C[2],
                            ratio=g.ratio, a=a)
        print(f"  saved {g.out}")
    print()
