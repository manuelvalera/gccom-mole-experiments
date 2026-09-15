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


T0 = time.time()
_LAST_BEAT = [time.time()]


def el():
    """elapsed wall time since launch, as mm:ss or h:mm:ss"""
    s_ = time.time() - T0
    h, r = divmod(int(s_), 3600)
    m, sec = divmod(r, 60)
    return f"{h}:{m:02d}:{sec:02d}" if h else f"{m:02d}:{sec:02d}"


def beat(msg, every=None, force=False):
    """print a heartbeat at most every `every` seconds (default 5 min)"""
    now = time.time()
    if every is None:
        every = globals().get('_BEAT_EVERY', 300.0)
    if force or now - _LAST_BEAT[0] >= every:
        _LAST_BEAT[0] = now
        print(f"    [{el()} elapsed] {msg}", flush=True)


def pull(oc, nm):
    oc.eval(f"[ii,jj,vv]=find({nm}); sz=size({nm});")
    return sp.csr_matrix(
        (oc.pull('vv').ravel(),
         (oc.pull('ii').ravel().astype(int) - 1,
          oc.pull('jj').ravel().astype(int) - 1)),
        shape=tuple(oc.pull('sz').ravel().astype(int)))



# Map variants, chosen to separate the two candidate instabilities:
#
#   full : X=p+A sin2piq, Y=q+A sin2pip, Z=r+A/2 sin2pip sin2piq
#          both effects active (the original).
#   horiz: X=p+A sin2piq, Y=q+A sin2pip, Z=r
#          zeta IS vertical, so buoyancy is added along the right direction;
#          but the u/v faces are curved surfaces, so the no-flux masks sit on
#          the wrong places.  Isolates the BOUNDARY-MASK candidate.
#   vert : X=p, Y=q, Z=r+A/2 sin2pip sin2piq
#          u/v faces are true planes, so the masks are correct; but zeta is
#          NOT vertical, so buoyancy is applied along a tilted direction.
#          Isolates the BUOYANCY-DIRECTION candidate.
MAPS = ('full', 'horiz', 'vert')


def build(mole, mode, a, amp, patched, korder=2, alpha=1e-4, gmap='full'):
    oc = Oct2Py(); oc.eval("more off; warning('off','all');")
    oc.addpath(mole)
    if patched:
        d = os.path.join(tempfile.gettempdir(), 'cone3d_patch')
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'GI13.m'), 'w').write(PATCHED_GI13)
        oc.addpath(d)
    which = oc.feval("which", "GI13").strip()

    N = a + 1
    XY = 0.0 if gmap == "vert" else 1.0
    ZZ = 0.0 if gmap == "horiz" else 1.0
    oc.eval(f"k={korder}; a={a}; N={N}; amp={amp}; h=1/a; XY={XY}; ZZ={ZZ};")
    if mode == 'cart':
        oc.eval("G=grad3D(k,a,h,a,h,a,h); D=div3D(k,a,h,a,h,a,h);")
    else:
        # NOTE -- Legacy pair is hard-wired on purpose.
        #
        # On MOLE HEAD the 4-arg calls dispatch differently than they look:
        #   grad3DCurv(k,X,Y,Z) -> grad3DCurvLegacy   (19-pt, uses GI13)
        #   div3DCurv (k,X,Y,Z) -> div3DCurvLegacy, which is NOT the same
        #                          operator as the 2021 div3DCurv.
        # Measured at a=20: grad is byte-identical between the two versions,
        # but div3DCurv goes from 48,000 to 172,564 nnz (max|diff| = 23.6),
        # taking L from 12.8 to 39.5 nnz/row.  Composing HEAD's dense div
        # with grad3DCurvLegacy is not an adjoint pair, so D*G is not a
        # Laplacian: the projection stops projecting and the run diverges
        # (KE ~ 1e45 within one forcing period).
        #
        # GI13 is called only from grad3DCurvLegacy, so the Legacy pair is the
        # correct -- and only meaningful -- configuration for testing the patch.
        # If you want to exercise HEAD's newer 6-arg curvilinear path, that is a
        # different experiment and needs its own validation.
        oc.eval("""
        [p,q,r]=meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
        X=p+amp*sin(2*pi*q)*XY; Y=q+amp*sin(2*pi*p)*XY;
        Z=r+0.5*amp*sin(2*pi*p).*sin(2*pi*q)*ZZ;
        if exist('grad3DCurvLegacy') == 2
          G=grad3DCurvLegacy(k,X,Y,Z); D=div3DCurvLegacy(k,X,Y,Z);
        else
          G=grad3DCurv(k,X,Y,Z); D=div3DCurv(k,X,Y,Z);   % 2021 snapshot
        end
        """)
    oc.eval(f"B=robinBC3D(k,a,h,a,h,a,h,{alpha},1);")
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
    # Computational (logical) coords are returned too: the sponge MUST be built
    # in computational space.  In physical space a curvilinear map pushes points
    # outside the unit box (measured x in [-0.087, 1.087] at amp=0.10), so the
    # distance-to-wall goes negative, the sponge coefficient hits 7.5 instead of
    # capping at 1, dt*S/tau exceeds its stability limit, and the sponge itself
    # drives the blow-up.
    comp = (fl(Pw), fl(Qw), fl(Rw), fl(Pc), fl(Qc), fl(Rc))

    def dmap(P, Q, R):
        """analytic derivatives of the map, (X_xi,X_eta,X_zeta), etc."""
        tp = 2*np.pi; one = np.ones_like(P); zero = np.zeros_like(P)
        if mode == 'cart':
            return (one, zero, zero), (zero, one, zero), (zero, zero, one)
        return ((one,              XY*amp*tp*np.cos(tp*Q), zero),
                (XY*amp*tp*np.cos(tp*P), one,              zero),
                (ZZ*0.5*amp*tp*np.cos(tp*P)*np.sin(tp*Q),
                 ZZ*0.5*amp*tp*np.sin(tp*P)*np.cos(tp*Q), one))

    return ops, (fl(Xw), fl(Yw), fl(Zw)), (fl(Xc), fl(Yc), fl(Zc)), which, comp, dmap


def run(mole, mode, a, amp, ratio, nper, patched, sigs, spg,
        dtfac, spp, taufac, korder, nu_visc, alpha, gmap='full',
        nomask=False, contra=False):
    print(f"  [{el()}] building operators in Octave "
          f"({'curvilinear -- 19-point stencil, slow' if mode == 'curv' else 'Cartesian'}) ...",
          flush=True)
    (G, D, B, Ic, Idf), Wp, Cp, which, comp, dmap = build(mole, mode, a, amp,
                                                    patched, korder, alpha, gmap)
    m = n = o = a
    nu, nv, nw = (m+1)*n*o, m*(n+1)*o, m*n*(o+1)
    nc = (m + 2) ** 3
    print(f"  [{el()}] GI13 : {which}")
    print(f"  grid : {a}^3   faces {nu+nv+nw}   centres {nc}   k={korder}")

    L = (D @ G + B).tocsc()
    asym = abs(L - L.T).max() / abs(L).max()
    # Pure Neumann makes L singular.  Abouali & Castillo (2014) Eq.10 give the
    # Robin form a*f + b*grad(f).n = g; a small alpha regularises the null
    # space.  Measured: alpha=1e-4 takes cond 1.3e14 -> 1.6e2 and the
    # projection residual to machine epsilon.  No pinning needed.
    print(f"  [{el()}] L: n={L.shape[0]} nnz={L.nnz} ({L.nnz/L.shape[0]:.1f}/row) "
          f"nonsym {asym:.2e}  Robin alpha={alpha:g}")
    print(f"  [{el()}] factorising -- this is the slow step "
          f"({'curvilinear: ~15 nnz/row' if mode == 'curv' else 'Cartesian: ~6.5 nnz/row'}) ...",
          flush=True)
    t = time.time(); solve = None
    try:
        import pypardiso
        _ps = pypardiso.PyPardisoSolver(); _Lc = L.tocsr(); _ps.factorize(_Lc)
        print(f"  [{el()}] PARDISO factorise : {time.time()-t:.1f}s")
        solve = lambda r: _ps.solve(_Lc, r)
    except Exception as e:
        print(f"  (pypardiso unavailable: {str(e)[:45]} -- SuperLU)")
    if solve is None:
        lu = spl.splu(L)
        print(f"  [{el()}] splu : {time.time()-t:.1f}s")
        solve = lu.solve

    # Is D*G actually a Laplacian?  The projection residual below cannot see
    # a mismatched (D, G) pair -- it holds by construction for ANY D and G.
    # Apply L to a field with a known Laplacian instead.
    try:
        s_ = np.linspace(0, 1, a + 1)
        cc_ = np.concatenate(([0], s_[:-1] + 0.5 / a, [1]))
        Pq, Qq, Rq = np.meshgrid(cc_, cc_, cc_, indexing='ij')
        if mode == 'curv':      # map to PHYSICAL coords, same map as build()
            Xq = Pq + amp*np.sin(2*np.pi*Qq)
            Yq = Qq + amp*np.sin(2*np.pi*Pq)
            Zq = Rq + 0.5*amp*np.sin(2*np.pi*Pq)*np.sin(2*np.pi*Qq)
        else:
            Xq, Yq, Zq = Pq, Qq, Rq
        fq = (Xq**2 + Yq**2 + Zq**2).transpose(2, 1, 0).ravel()
        vq = ((D @ G) @ fq).reshape(a+2, a+2, a+2)[3:-3, 3:-3, 3:-3] - 6.0
        lap_err = float(np.sqrt(np.mean(vq**2)))
        print(f"  [{el()}] D*G on f=x^2+y^2+z^2 (true 6): rms err {lap_err:.3e}  "
              f"{'OK' if lap_err < 1.0 else '*** D and G ARE NOT AN ADJOINT PAIR ***'}")
    except Exception as e:
        print(f"  (Laplacian check skipped: {str(e)[:40]})")

    fk = np.random.default_rng(1).standard_normal(G.shape[0])
    rr = D @ fk
    dv = np.linalg.norm(D @ (fk - G @ solve(rr))) / np.linalg.norm(rr)
    print(f"  [{el()}] check: divergence after projection {dv:.3e}  "
          f"{'OK' if dv < 1e-8 else '*** PROJECTION IS NOT WORKING ***'}")
    sys.stdout.flush()

    Nb = 1.0; om = ratio * Nb
    xw, yw, zw = Wp; xc, yc, zc = Cp
    pw, qw, rw, pc, qc, rc = comp          # computational coords, always [0,1]
    dW = np.minimum.reduce([pw, 1-pw, qw, 1-qw, rw, 1-rw])
    dC = np.minimum.reduce([pc, 1-pc, qc, 1-qc, rc, 1-rc])
    Sw = np.clip((np.maximum(0, spg - dW) / spg) ** 2, 0.0, 1.0)
    Sc = np.clip((np.maximum(0, spg - dC) / spg) ** 2, 0.0, 1.0)
    lim = dt_probe = min(dtfac/Nb, 2*np.pi/(ratio*Nb)/spp)
    print(f"  [{el()}] sponge: max coeff {max(Sw.max(), Sc.max()):.2f} (must be <=1), "
          f"dt*S/tau = {dt_probe*max(Sw.max(), Sc.max())/(taufac/Nb):.2f} "
          f"{'OK' if dt_probe*max(Sw.max(),Sc.max())/(taufac/Nb) < 2 else '*** UNSTABLE ***'}")
    tau = taufac / Nb

    T = 2*np.pi/om; dt = min(dtfac/Nb, T/spp); nt = int(nper*T/dt)
    print(f"  dt   : {dt:.4g}  (dt*N={dt*Nb:.3g})  {nt} steps, {T/dt:.0f}/period")

    Ic_z  = Ic[nu+nv:, :].tocsr()
    Idf_u = Idf[:, :nu].tocsr()
    Idf_v = Idf[:, nu:nu+nv].tocsr()
    Idf_w = Idf[:, nu+nv:].tocsr()
    if nu_visc > 0:
        LAPc = (D @ G).tocsr()
        Ic_x = Ic[:nu, :].tocsr(); Ic_y = Ic[nu:nu+nv, :].tocsr()

    # ---------------------------------------------------------------------
    # No-flux at the physical boundary.
    #
    # Zeroing the CARTESIAN component at faces of constant xi/eta/zeta is only
    # correct on a Cartesian grid.  On a curved grid those faces are sheets
    # that wander up to `amp` off the true boundary, so the constraint lands on
    # the wrong surfaces -- measured to destabilise the run above amp ~ 0.05.
    #
    # Correct condition: annihilate the CONTRAVARIANT flux n.(u,v,w), with n
    # the metric cofactor normal to that coordinate surface:
    #     xi   faces   n = r_eta x r_zeta = (A,B,C)   ->  u = -(B v + C w)/A
    #     eta  faces   n = r_zeta x r_xi  = (D,E,F)   ->  v = -(D u + F w)/E
    #     zeta faces   n = r_xi x r_eta   = (G,H,I)   ->  w = -(G u + H v)/I
    # The two components not stored on that face are averaged in from their
    # own face sets.  Diagonal cofactors stay O(1) (min 0.74 at amp=0.10), so
    # each solve is well conditioned.
    # ---------------------------------------------------------------------
    sN = np.linspace(0, 1, a + 1); cN = sN[:-1] + 0.5 / a
    iu3 = np.arange(nu).reshape(o, n, m+1)
    iv3 = np.arange(nv).reshape(o, n+1, m)
    iw3 = np.arange(nw).reshape(o+1, n, m)

    def cof(Pg, Qg, Rg):
        (Xx, Xe, Xz), (Yx, Ye, Yz), (Zx, Ze, Zz) = dmap(Pg, Qg, Rg)
        A_ = Ye*Zz - Ze*Yz; B_ = Ze*Xz - Xe*Zz; C_ = Xe*Yz - Ye*Xz
        D_ = Yz*Zx - Zz*Yx; E_ = Zz*Xx - Xz*Zx; F_ = Xz*Yx - Yz*Xx
        G_ = Yx*Ze - Zx*Ye; H_ = Zx*Xe - Xx*Ze; I_ = Xx*Ye - Yx*Xe
        return (A_, B_, C_), (D_, E_, F_), (G_, H_, I_)

    bnd = []
    for side, ip_ in ((0, 0), (1, m)):                       # xi faces
        Qg, Pg = np.meshgrid(cN, cN, indexing='ij')          # (k,j) over centres
        Kg, Jg = np.meshgrid(np.arange(o), np.arange(n), indexing='ij')
        (A_, B_, C_), _, _ = cof(np.full_like(Qg, sN[ip_]), Pg, Qg)
        tgt = iu3[:, :, ip_].ravel()
        ic = 0 if side == 0 else m-1
        vsrc = (iv3[Kg, Jg, ic].ravel(), iv3[Kg, Jg+1, ic].ravel())
        wsrc = (iw3[Kg, Jg, ic].ravel(), iw3[Kg+1, Jg, ic].ravel())
        bnd.append(('u', tgt, vsrc, wsrc, B_.ravel()/A_.ravel(), C_.ravel()/A_.ravel()))
    for side, jp_ in ((0, 0), (1, n)):                       # eta faces
        Qg, Pg = np.meshgrid(cN, cN, indexing='ij')
        Kg, Ig = np.meshgrid(np.arange(o), np.arange(m), indexing='ij')
        _, (D_, E_, F_), _ = cof(Pg, np.full_like(Qg, sN[jp_]), Qg)
        tgt = iv3[:, jp_, :].ravel()
        jc = 0 if side == 0 else n-1
        usrc = (iu3[Kg, jc, Ig].ravel(), iu3[Kg, jc, Ig+1].ravel())
        wsrc = (iw3[Kg, jc, Ig].ravel(), iw3[Kg+1, jc, Ig].ravel())
        bnd.append(('v', tgt, usrc, wsrc, D_.ravel()/E_.ravel(), F_.ravel()/E_.ravel()))
    for side, kp_ in ((0, 0), (1, o)):                       # zeta faces
        Qg, Pg = np.meshgrid(cN, cN, indexing='ij')
        Jg, Ig = np.meshgrid(np.arange(n), np.arange(m), indexing='ij')
        _, _, (G_, H_, I_) = cof(Pg, Qg, np.full_like(Qg, sN[kp_]))
        tgt = iw3[kp_, :, :].ravel()
        kc = 0 if side == 0 else o-1
        usrc = (iu3[kc, Jg, Ig].ravel(), iu3[kc, Jg, Ig+1].ravel())
        vsrc = (iv3[kc, Jg, Ig].ravel(), iv3[kc, Jg+1, Ig].ravel())
        bnd.append(('w', tgt, usrc, vsrc, G_.ravel()/I_.ravel(), H_.ravel()/I_.ravel()))

    def apply_bc(uu, vv, ww):
        """set the stored component on each boundary face so n.(u,v,w) = 0"""
        for kind, tgt, s1, s2, c1, c2 in bnd:
            if kind == 'u':
                a1 = 0.5*(vv[s1[0]] + vv[s1[1]]); a2 = 0.5*(ww[s2[0]] + ww[s2[1]])
                uu[tgt] = -(c1*a1 + c2*a2)
            elif kind == 'v':
                a1 = 0.5*(uu[s1[0]] + uu[s1[1]]); a2 = 0.5*(ww[s2[0]] + ww[s2[1]])
                vv[tgt] = -(c1*a1 + c2*a2)
            else:
                a1 = 0.5*(uu[s1[0]] + uu[s1[1]]); a2 = 0.5*(vv[s2[0]] + vv[s2[1]])
                ww[tgt] = -(c1*a1 + c2*a2)
        return uu, vv, ww

    if nomask:
        um = np.ones(nu); vm = np.ones(nv); wm = np.ones(nw)
        apply_bc = lambda uu, vv, ww: (uu, vv, ww)
        print(f"  [{el()}] boundary: OFF (sponge only)")
    elif contra:
        um = np.ones(nu); vm = np.ones(nv); wm = np.ones(nw)
        print(f"  [{el()}] boundary: CONTRAVARIANT no-flux "
              f"(n.(u,v,w)=0 from metric cofactors)")
    else:
        um = np.ones(nu); vm = np.ones(nv); wm = np.ones(nw)
        um[iu3[:, :, 0].ravel()] = 0; um[iu3[:, :, -1].ravel()] = 0
        vm[iv3[:, 0, :].ravel()] = 0; vm[iv3[:, -1, :].ravel()] = 0
        wm[iw3[0].ravel()] = 0; wm[iw3[-1].ravel()] = 0
        apply_bc = lambda uu, vv, ww: (uu, vv, ww)
        print(f"  [{el()}] boundary: Cartesian masks "
              f"(correct only on a Cartesian grid)")

    cst = np.array([0, 1, 2, 3, 4, 2, 3, 4, 5, 6]) / 6.0    # SSPRK104 abscissae
    results = []
    for sig in sigs:
        print(f"\n  --- sig={sig}: {sig*a:.1f} cells wide, "
              f"{(0.5-spg)/sig:.1f} source-widths of travel ---")
        if nu_visc > 0:
            print(f"  viscosity nu={nu_visc:g}")
        print("  scheme: SSPRK104 + incremental projection, 1 solve/step")
        sys.stdout.flush()
        Fw = 1e-3*np.exp(-(((xw-.5)**2 + (yw-.5)**2 + (zw-.5)**2)/sig**2))
        u = np.zeros(nu); v = np.zeros(nv); w = np.zeros(nw); b = np.zeros(nc)
        gp = np.zeros(nu+nv+nw)
        Uc = np.zeros(nc, complex); Vc = np.zeros(nc, complex)
        Wc = np.zeros(nc, complex)

        def F(uu, vv, ww, bb, ts):
            du  = -gp[:nu]
            dv_ = -gp[nu:nu+nv]
            dw  = -gp[nu+nv:] + Ic_z@bb + Fw*np.sin(om*ts) - Sw*ww/tau
            db  = -Nb**2*(Idf_w@ww) - Sc*bb/tau
            if nu_visc > 0:
                du  = du  + nu_visc*(Ic_x @ (LAPc @ (Idf_u@uu)))
                dv_ = dv_ + nu_visc*(Ic_y @ (LAPc @ (Idf_v@vv)))
                dw  = dw  + nu_visc*(Ic_z @ (LAPc @ (Idf_w@ww)))
                db  = db  + nu_visc*(LAPc@bb)
            return du, dv_, dw, db

        t0 = time.time(); ref_w = [None]; diverged = False; div_tT = None
        beat(f"time loop starting: {nt} steps", force=True)
        for it in range(1, nt+1):
            tn = (it-1)*dt
            q1 = (u.copy(), v.copy(), w.copy(), b.copy())
            q2 = tuple(x.copy() for x in q1)
            for i in range(5):
                f = F(*q1, tn+cst[i]*dt)
                q1 = tuple(q1[j] + dt*f[j]/6.0 for j in range(4))
            q2 = tuple(q2[j]/25.0 + 9.0*q1[j]/25.0 for j in range(4))
            q1 = tuple(15.0*q2[j] - 5.0*q1[j] for j in range(4))
            for i in range(5, 9):
                f = F(*q1, tn+cst[i]*dt)
                q1 = tuple(q1[j] + dt*f[j]/6.0 for j in range(4))
            f = F(*q1, tn+cst[9]*dt)
            q1 = tuple(q2[j] + 0.6*q1[j] + 0.1*dt*f[j] for j in range(4))
            us, vs, ws, b = q1
            us = us*um; vs = vs*vm; ws = ws*wm

            rhs = (D @ np.concatenate([us, vs, ws]))/dt
            gphi = G @ solve(rhs)
            u = (us - dt*gphi[:nu])*um
            v = (vs - dt*gphi[nu:nu+nv])*vm
            w = (ws - dt*gphi[nu+nv:])*wm
            u, v, w = apply_bc(u, v, w)
            gp = gp + gphi

            tt = it*dt
            if it % 20 == 0:
                mw = np.abs(w).max()
                if ref_w[0] is None and tt/T > 0.5:
                    ref_w[0] = max(mw, 1e-30)
                bad = (not np.isfinite(mw)) or \
                      (ref_w[0] is not None and mw > 1e4*ref_w[0])
                if bad:
                    print(f"\n  *** DIVERGED at step {it} (t/T={tt/T:.2f}): "
                          f"max|w|={mw:.3e}, {mw/ref_w[0]:.1e}x its value at "
                          f"t/T=0.5 -- aborting this sig ***\n", flush=True)
                    diverged = True; div_tT = tt/T
                    break
            if it % 20 == 0:
                rate = (time.time() - t0) / it
                beat(f"step {it}/{nt}  ({100*it/nt:.0f}%)  "
                     f"ETA {(nt-it)*rate/60:.1f} min")
            if it % max(1, int(round(T/dt))) == 0:
                ke = 0.5*(u@u + v@v + w@w); pe = 0.5*(b@b)/Nb**2
                print(f"      [{el()}] period {int(round(tt/T)):3d}  "
                      f"KE {ke:.4e}  PE {pe:.4e}  "
                      f"|grad p| {np.abs(gp).max():.3e}", flush=True)
                _LAST_BEAT[0] = time.time()
            if tt/T > nper - 2:
                # Internal-wave velocity is ALONG the beam (c_p . c_g = 0), so
                # |w| alone carries a sin^2(theta) weighting that biases an
                # angular peak-find toward steeper angles.  Use TOTAL speed.
                ph = np.exp(-1j*om*tt)*dt
                Uc += (Idf_u@u)*ph; Vc += (Idf_v@v)*ph; Wc += (Idf_w@w)*ph
        spd = np.sqrt(np.abs(Uc)**2 + np.abs(Vc)**2 + np.abs(Wc)**2)
        if diverged or not np.isfinite(spd).all() or spd.max() == 0:
            print(f"  RESULT: DIVERGED at t/T={div_tT:.2f}" if div_tT
                  else "  RESULT: no usable signal")
            spd = None
        results.append((sig, spd))
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
    print(f"\n  angular profile of |SPEED| (not |w|), deg from horizontal")
    for q, vv in zip(range(0, 91, 5), prof):
        bar = '#'*int(45*vv/mx) if np.isfinite(vv) else ''
        print(f"    {q:3d} {bar}{'  <== theory' if abs(q-tgt) < 2.5 else ''}")

    R3 = np.sqrt((xc-.5)**2 + (yc-.5)**2 + (zc-.5)**2)
    A3 = np.degrees(np.arctan2(np.abs(zc-.5), np.hypot(xc-.5, yc-.5)))
    def ridge(rr, half=0.011):
        s_ = np.abs(R3-rr) < half
        if s_.sum() < 60: return None
        aa, vv = A3[s_], F[s_]
        k_ = (aa > 5) & (aa < 85)
        if k_.sum() < 30: return None
        g = np.arange(5, 85, 0.5)
        sm = np.array([vv[k_][np.abs(aa[k_]-q) < 2.5].mean()
                       if (np.abs(aa[k_]-q) < 2.5).sum() else np.nan for q in g])
        return None if np.all(np.isnan(sm)) else (g[np.nanargmax(sm)],
                                                  np.nanmax(sm))
    print(f"\n  ridge angle vs radius   (near field contaminates small r)")
    print(f"    {'r':>6} {'ridge deg':>10} {'peak val':>11}")
    for rr in np.arange(0.12, 0.42, 0.03):
        got = ridge(rr)
        if got: print(f"    {rr:>6.2f} {got[0]:>10.2f} {got[1]:>11.3e}")

    est = [ridge(rr, 0.012)[0] for rr in np.linspace(0.20, 0.36, 10)
           if ridge(rr, 0.012)]
    meas = float(np.median(est)) if est else float('nan')
    print(f"\n  CONE HALF-ANGLE  measured {meas:.2f}   theory {tgt:.2f}   "
          f"error {meas-tgt:+.2f} deg ({100*(meas-tgt)/tgt:+.1f}%)")
    return meas


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=['cart', 'curv'])
    ap.add_argument('--mole', required=True)
    ap.add_argument('--a', type=int, default=64)
    ap.add_argument('--ratio', type=float, default=0.6)
    ap.add_argument('--amp', type=float, default=0.10)
    ap.add_argument('--nper', type=int, default=6)
    ap.add_argument('--sig', type=float, nargs='+', default=[0.07])
    ap.add_argument('--spg', type=float, default=0.05)
    ap.add_argument('--tau', type=float, default=0.5)
    ap.add_argument('--dtfac', type=float, default=2.0)
    ap.add_argument('--spp', type=int, default=60)
    ap.add_argument('--order', type=int, default=2, choices=[2, 4, 6])
    ap.add_argument('--nu', type=float, default=0.0)
    ap.add_argument('--contra', action='store_true',
                    help='contravariant no-flux BC from metric cofactors')
    ap.add_argument('--nomask', action='store_true',
                    help='drop the hard no-flux masks; rely on the sponge')
    ap.add_argument('--map', dest='gmap', default='full',
                    choices=['full', 'horiz', 'vert'],
                    help='curvilinear map variant (see MAPS note in source)')
    ap.add_argument('--alpha', type=float, default=1e-4,
                    help='Robin coefficient regularising the pressure system')
    ap.add_argument('--patch', action='store_true')
    ap.add_argument('--out', default=None)
    ap.add_argument('--beat', type=float, default=300.0,
                    help='heartbeat interval in seconds (default 300 = 5 min)')
    g = ap.parse_args()
    _BEAT_EVERY = g.beat
    print(f"\n=== cone3d {g.mode}  a={g.a}  omega/N={g.ratio}  nper={g.nper}  "
          f"k={g.order}  nu={g.nu:g}  alpha={g.alpha:g}  "
          f"map={g.gmap}  {'PATCHED' if g.patch else 'STOCK'} GI13 ===")
    res, C, a = run(g.mole, g.mode, g.a, g.amp, g.ratio, g.nper, g.patch,
                    g.sig, g.spg, g.dtfac, g.spp, g.tau, g.order, g.nu, g.alpha, g.gmap, g.nomask, g.contra)
    print("\n" + "="*62 + "\n  SUMMARY\n" + "="*62)
    for sig, A in res:
        print(f"\n>>>>>> source sig = {sig}")
        if A is None:
            print("  no usable field (run diverged)")
            continue
        report(A, C, a, g.ratio)
    if g.out:
        np.savez_compressed(g.out, sigs=[s for s, _ in res],
                            **{f"A{i}": A for i, (_, A) in enumerate(res)
                               if A is not None},
                            xc=C[0], yc=C[1], zc=C[2], ratio=g.ratio, a=a)
        print(f"\n  saved {g.out}")
    print(f"\n  total wall time {el()}\n")
