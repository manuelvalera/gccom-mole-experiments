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

# MKL reads thread settings at import time, so this must precede numpy.
# On hybrid Intel CPUs (P-cores + E-cores) the best value is the number of
# PHYSICAL P-CORES: MKL splits work statically and E-cores become stragglers.
# Measured on an i7-14700F (8 P + 12 E): 6-8 threads = 1.77x, 12+ regresses.
_nt = os.environ.get('CONE3D_THREADS')
if _nt is None:
    _nt = '8'
for _v in ('MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, _nt)
os.environ.setdefault('MKL_DYNAMIC', 'FALSE')

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


def run(mole, mode, a, amp, ratio, nper, patched, sigs, spg, dtfac=0.02, spp=300,
        taufac=0.5):
    (G, D, B, Ic, Idf), Wp, Cp, which = build(mole, mode, a, amp, patched)
    m = n = o = a
    nu, nv, nw = (m + 1) * n * o, m * (n + 1) * o, m * n * (o + 1)
    nc = (m + 2) ** 3
    print(f"  GI13 : {which}")
    print(f"  grid : {a}^3   faces {nu+nv+nw}   centres {nc}")

    L = (D @ G + B).tocsc()
    asym = abs(L - L.T).max() / abs(L).max()

    # L is singular (pure Neumann): nullity 1.  Pinning one node removes it ONLY
    # if the LEFT null vector u is nonzero there -- pinning at the geometric
    # centre silently fails.  Find u as the residual of a least-squares fit to a
    # random vector: r = b - L x is orthogonal to range(L), so L^T r = 0.
    t = time.time()
    bb = np.random.default_rng(7).standard_normal(nc)
    u = bb - L @ spl.lsqr(L, bb, atol=1e-10, btol=1e-10, iter_lim=3000)[0]
    ip = int(np.argmax(np.abs(u)))
    print(f"  L    : nonsymmetry {asym:.2e}   pin node {ip} "
          f"(|u|={np.abs(u[ip])/np.linalg.norm(u):.2e}, centre would be "
          f"{np.abs(u[nc//2])/np.linalg.norm(u):.2e})   [{time.time()-t:.1f}s]")

    L = L.tolil(); L[ip, :] = 0; L[ip, ip] = 1; L = L.tocsc()

    # MKL PARDISO is ~7x faster to factorise than SuperLU; fall back if absent.
    solve = None
    try:
        import pypardiso
        t = time.time()
        _ps = pypardiso.PyPardisoSolver(); _Lc = L.tocsr(); _ps.factorize(_Lc)
        print(f"  PARDISO factorise : {time.time()-t:.1f}s")
        solve = lambda r: _ps.solve(_Lc, r)
    except Exception as e:
        print(f"  (pypardiso unavailable: {str(e)[:50]} -- using SuperLU)")
    if solve is None:
        t = time.time(); lu = spl.splu(L)
        print(f"  splu : {time.time()-t:.1f}s")
        solve = lu.solve
    print("  (factorisation done once, reused for every source width)")

    # verify the projection actually removes divergence before trusting anything
    fk = np.random.default_rng(1).standard_normal(G.shape[0])
    rr = D @ fk; rr[ip] = 0
    pk = solve(rr)
    dv = np.linalg.norm(D @ (fk - G @ pk)) / np.linalg.norm(D @ fk)
    print(f"  check: divergence after projection {dv:.3e}  "
          f"{'OK' if dv < 1e-8 else '*** PROJECTION IS NOT WORKING ***'}")
    sys.stdout.flush()

    Nb = 1.0; om = ratio * Nb
    xw, yw, zw = Wp; xc, yc, zc = Cp
    dW = np.minimum.reduce([xw, 1-xw, yw, 1-yw, zw, 1-zw])
    dC = np.minimum.reduce([xc, 1-xc, yc, 1-yc, zc, 1-zc])
    Sw = (np.maximum(0, spg - dW) / spg) ** 2
    Sc = (np.maximum(0, spg - dC) / spg) ** 2
    tau = taufac / Nb

    T = 2*np.pi/om; dt = min(dtfac/Nb, T/spp); nt = int(nper*T/dt)
    print(f"  dt   : {dt:.4g}  (dt*N={dt*Nb:.3g}, stability limit 2)  {nt} steps for {nper} periods, {T/dt:.0f}/period")
    iu = np.arange(nu).reshape(o, n, m+1); um = np.ones(nu)
    um[iu[:, :, 0].ravel()] = 0; um[iu[:, :, -1].ravel()] = 0
    iv = np.arange(nv).reshape(o, n+1, m); vm = np.ones(nv)
    vm[iv[:, 0, :].ravel()] = 0; vm[iv[:, -1, :].ravel()] = 0
    iw = np.arange(nw).reshape(o+1, n, m); wm = np.ones(nw)
    wm[iw[0].ravel()] = 0; wm[iw[-1].ravel()] = 0
    Ic_z  = Ic[nu+nv:, :].tocsr()      # centres -> z-faces only
    Idf_w = Idf[:, nu+nv:].tocsr()     # z-faces -> centres only

    results = []
    for sig in sigs:
        cells = sig / (1.0 / a)
        travel = (0.5 - spg) / sig
        print(f"\n  --- sig={sig}: {cells:.1f} cells wide, "
              f"{travel:.1f} source-widths of travel before the sponge ---")
        print(f"  scheme: SSPRK104 + incremental projection   sponge tau={taufac}/N, "
              f"thickness {spg}")
        sys.stdout.flush()
        Fw = 1e-3 * np.exp(-(((xw-.5)**2 + (yw-.5)**2 + (zw-.5)**2) / sig**2))
        u = np.zeros(nu); v = np.zeros(nv); w = np.zeros(nw); b = np.zeros(nc)
        gp = np.zeros(nu + nv + nw)          # grad of the CURRENT pressure, lagged
        Wc = np.zeros(nc, complex)

        # SSPRK104 stage abscissae (low-storage form of Ketcheson 2008)
        cst = np.array([0, 1, 2, 3, 4, 2, 3, 4, 5, 6]) / 6.0

        def F(uu, vv, ww, bb, ts):
            """du/dt etc. with pressure held at its lagged value."""
            bz = Ic_z @ bb
            du = -gp[:nu]
            dv = -gp[nu:nu+nv]
            dw = -gp[nu+nv:] + bz + Fw * np.sin(om * ts) - Sw * ww / tau
            wc_ = Idf_w @ ww
            db = -Nb**2 * wc_ - Sc * bb / tau
            return du, dv, dw, db

        t0 = time.time()
        for it in range(1, nt + 1):
            tn = (it - 1) * dt
            # ---- SSPRK104, low storage: only q1 carries state through stages
            q1 = (u.copy(), v.copy(), w.copy(), b.copy())
            q2 = (u.copy(), v.copy(), w.copy(), b.copy())
            for i in range(5):
                f = F(*q1, tn + cst[i]*dt)
                q1 = tuple(q1[j] + dt*f[j]/6.0 for j in range(4))
            q2 = tuple(q2[j]/25.0 + 9.0*q1[j]/25.0 for j in range(4))
            q1 = tuple(15.0*q2[j] - 5.0*q1[j] for j in range(4))
            for i in range(5, 9):
                f = F(*q1, tn + cst[i]*dt)
                q1 = tuple(q1[j] + dt*f[j]/6.0 for j in range(4))
            f = F(*q1, tn + cst[9]*dt)
            q1 = tuple(q2[j] + 0.6*q1[j] + 0.1*dt*f[j] for j in range(4))
            us, vs, ws, b = q1
            us *= um; vs *= vm; ws *= wm

            # ---- incremental projection: solve only for the CORRECTION phi
            rhs = (D @ np.concatenate([us, vs, ws])) / dt; rhs[ip] = 0
            phi = solve(rhs)
            gphi = G @ phi
            u = (us - dt*gphi[:nu]) * um
            v = (vs - dt*gphi[nu:nu+nv]) * vm
            w = (ws - dt*gphi[nu+nv:]) * wm
            gp = gp + gphi                    # p^{n+1} = p^n + phi

            tt = it * dt
            if it % max(1, int(round(T/dt))) == 0:
                ke = 0.5*(u@u + v@v + w@w)
                pe = 0.5*(b@b)/Nb**2
                print(f"      period {int(round(tt/T)):3d}  KE {ke:.4e}  PE {pe:.4e}  "
                      f"E {ke+pe:.4e}  |grad p| {np.abs(gp).max():.3e}"); sys.stdout.flush()
            if tt / T > nper - 2:
                Wc += (Idf_w @ w) * np.exp(-1j * om * tt) * dt

        results.append((sig, np.abs(Wc)))
    return results, Cp, a


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
    ap.add_argument('--sig', type=float, nargs='+', default=[0.10, 0.06, 0.04])
    ap.add_argument('--spg', type=float, default=0.06)
    ap.add_argument('--tau', type=float, default=0.5,
                    help='sponge damping timescale = tau/N; smaller = stronger')
    ap.add_argument('--spp', type=int, default=300, help='timesteps per forcing period')
    ap.add_argument('--dtfac', type=float, default=0.02,
                    help='dt = dtfac/N; symplectic-Euler stability limit is 2')
    ap.add_argument('--patch', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--threads', default=None,
                    help='MKL threads; set CONE3D_THREADS instead if importing')
    g = ap.parse_args()
    print(f"\n=== cone3d {g.mode}  a={g.a}  omega/N={g.ratio}  nper={g.nper}  "
          f"sig={g.sig}  spg={g.spg}  {'PATCHED' if g.patch else 'STOCK'} GI13 ===")
    res, C, a = run(g.mole, g.mode, g.a, g.amp, g.ratio, g.nper, g.patch, g.sig, g.spg, g.dtfac, g.spp, g.tau)
    print("\n" + "=" * 62 + "\n  SUMMARY\n" + "=" * 62)
    for sig, A in res:
        print(f"\n>>>>>> source sig = {sig}")
        report(A, C, a, g.ratio)
    if g.out:
        np.savez_compressed(g.out, sigs=[s for s, _ in res],
                            **{f"A{i}": A for i, (_, A) in enumerate(res)},
                            xc=C[0], yc=C[1], zc=C[2], ratio=g.ratio, a=a)
        print(f"\n  saved {g.out}")
    print()
