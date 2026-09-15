#!/usr/bin/env python3
"""Assembled-Laplacian test: does D*G reproduce a known Laplacian on a
curvilinear grid?  Stock GI13 vs the patch.  No time stepping involved."""
import os, sys, tempfile, numpy as np, scipy.sparse as sp
from oct2py import Oct2Py
PATCH = r"""
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
def get(oc, nm):
    oc.eval(f"[ii,jj,vv]=find({nm}); sz=size({nm});")
    return sp.csr_matrix((oc.pull('vv').ravel(),
        (oc.pull('ii').ravel().astype(int)-1, oc.pull('jj').ravel().astype(int)-1)),
        shape=tuple(oc.pull('sz').ravel().astype(int)))

def trial(mole, a, amp, patched, korder=2):
    oc = Oct2Py(); oc.eval("more off; warning('off','all');"); oc.addpath(mole)
    if patched:
        d = os.path.join(tempfile.gettempdir(), 'lapt'); os.makedirs(d, exist_ok=True)
        open(os.path.join(d, 'GI13.m'), 'w').write(PATCH); oc.addpath(d)
    N = a + 1
    oc.eval(f"k={korder}; a={a}; N={N}; amp={amp};")
    oc.eval("""
    [p,q,r]=meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
    X=p+amp*sin(2*pi*q); Y=q+amp*sin(2*pi*p);
    Z=r+0.5*amp*sin(2*pi*p).*sin(2*pi*q);
    if exist('grad3DCurvLegacy')==2
      G=grad3DCurvLegacy(k,X,Y,Z); D=div3DCurvLegacy(k,X,Y,Z);
    else
      G=grad3DCurv(k,X,Y,Z); D=div3DCurv(k,X,Y,Z);
    end""")
    G, D = get(oc, 'G'), get(oc, 'D'); oc.exit()
    s = np.linspace(0, 1, N); c = np.concatenate(([0], s[:-1]+0.5/a, [1]))
    P_, Q_, R_ = np.meshgrid(c, c, c, indexing='ij')
    X = P_+amp*np.sin(2*np.pi*Q_); Y = Q_+amp*np.sin(2*np.pi*P_)
    Z = R_+0.5*amp*np.sin(2*np.pi*P_)*np.sin(2*np.pi*Q_)
    f = (X**2+Y**2+Z**2).transpose(2,1,0).ravel()
    v = ((D@G)@f).reshape(a+2,a+2,a+2)[3:-3,3:-3,3:-3] - 6.0
    return float(np.sqrt(np.mean(v**2))), float(np.abs(v).max())


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--mole', required=True, help='MOLE src/matlab_octave path')
    ap.add_argument('--a', type=int, nargs='+', default=[12, 16, 20, 24, 28])
    ap.add_argument('--amp', type=float, nargs='+', default=[0.0, 0.02, 0.05, 0.10, 0.15, 0.20])
    ap.add_argument('--k', type=int, nargs='+', default=[2, 4, 6])
    g = ap.parse_args()

    print("\nMOLE GI13: assembled 3-D curvilinear Laplacian L = D*G")
    print("applied to f = x^2+y^2+z^2 (exact value 6).  Static test, no solver.\n")

    print("=== A. grid refinement (amp=0.10, k=2) ===")
    print(f"{'a':>4} {'stock rms':>12} {'patched rms':>12} {'ratio':>8}"
          f" {'stock max':>12} {'patched max':>12}")
    rows = []
    for a in g.a:
        s = trial(g.mole, a, 0.10, False); p = trial(g.mole, a, 0.10, True)
        print(f"{a:>4} {s[0]:>12.4e} {p[0]:>12.4e} {s[0]/p[0]:>7.1f}x"
              f" {s[1]:>12.4e} {p[1]:>12.4e}")
        rows.append((a, s[0], p[0]))
    print("\nobserved order of convergence:")
    for i in range(1, len(rows)):
        r = np.log(rows[i][0] / rows[i-1][0])
        print(f"  {rows[i-1][0]}->{rows[i][0]}:  "
              f"stock {np.log(rows[i-1][1]/rows[i][1])/r:+5.2f}   "
              f"patched {np.log(rows[i-1][2]/rows[i][2])/r:+5.2f}")

    print("\n=== B. grid distortion (a=20, k=2) ===")
    print(f"{'amp':>6} {'stock rms':>12} {'patched rms':>12} {'ratio':>8}")
    for amp in g.amp:
        s = trial(g.mole, 20, amp, False); p = trial(g.mole, 20, amp, True)
        r = s[0]/p[0] if p[0] > 0 else float('inf')
        print(f"{amp:>6.2f} {s[0]:>12.4e} {p[0]:>12.4e} {r:>7.1f}x")

    print("\n=== C. operator order (a=20, amp=0.10) ===")
    print(f"{'k':>3} {'stock rms':>12} {'patched rms':>12} {'ratio':>8}")
    for k in g.k:
        try:
            s = trial(g.mole, 20, 0.10, False, k); p = trial(g.mole, 20, 0.10, True, k)
            print(f"{k:>3} {s[0]:>12.4e} {p[0]:>12.4e} {s[0]/p[0]:>7.1f}x")
        except Exception as e:
            print(f"{k:>3}  failed: {str(e)[:55]}")
    print()
