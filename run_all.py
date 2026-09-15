#!/usr/bin/env python3
"""
run_all.py -- runs every MOLE check in one go and prints a single report.

Usage:
    python run_all.py --mole /path/to/mole/src/matlab_octave

Prints a self-contained block you can paste anywhere.  Takes about a minute.
"""
import argparse, os, sys, time, textwrap
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def banner(t):
    print("\n" + "=" * 74)
    print(t)
    print("=" * 74)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mole", required=True,
                    help="path to MOLE's matlab_octave directory")
    ap.add_argument("--patch", default=os.path.join(ROOT, "patch"),
                    help="directory holding the patched GI13.m")
    ap.add_argument("--grids", default=os.path.join(ROOT, "grids"))
    a = ap.parse_args()

    from oct2py import Oct2Py
    oc = Oct2Py()
    oc.eval("more off;")
    oc.addpath(a.mole)
    print("MOLE     :", a.mole)
    print("GI13 seen:", oc.feval("which", "GI13"))
    try:
        v = oc.eval("version", verbose=False)
    except Exception:
        v = "?"
    print("Octave   :", v)

    # ---------------------------------------------------------------- 1
    banner("1.  GI13 index map  (no numerics -- pure index bookkeeping)")
    oc.eval("m=4; n=3; o=3;")
    oc.eval("[I,J,K]=ndgrid(1:m,1:n+1,1:o); v=I(:)+100*J(:)+10000*K(:);")
    oc.eval("W=reshape(GI13(speye(m*(n+1)*o),m,n,o,'Gn')*v, m+1,n,o);")
    W = oc.pull("W")
    bad = 0
    for k in range(3):
        for j in range(3):
            for i in range(5):
                val = W[i, j, k]
                if val == 0:
                    bad += 1
                    continue
                kk = int(val // 10000)
                jj = int((val - 10000 * kk) // 100)
                ii = int(val - 10000 * kk - 100 * jj)
                if kk != k + 1:
                    if bad < 5:
                        print(f"    out({i+1},{j+1},{k+1}) <- src({ii},{jj},{kk})"
                              f"   *** wrong zeta-plane")
                    bad += 1
    print(f"\n    {bad} of 45 output entries read from the wrong zeta-plane")
    print("    (expected 15 with stock GI13, 0 with the patch)")

    # ---------------------------------------------------------------- 2
    banner("2.  grad3DCurv convergence  (curvilinear grid, k=2)")
    setup = """
    function [X,Y,Z]=cn(N,amp)
      [a,b,c]=meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
      X=a+amp*sin(2*pi*b); Y=b+amp*sin(2*pi*a); Z=c+0.5*amp*sin(2*pi*a).*sin(2*pi*b);
    end
    function [x,y,z]=cc(N,amp)
      s=linspace(0,1,N); s=[0, s(1:end-1)+0.5/(N-1), 1];
      [a,b,c]=meshgrid(s,s,s);
      x=a+amp*sin(2*pi*b); y=b+amp*sin(2*pi*a); z=c+0.5*amp*sin(2*pi*a).*sin(2*pi*b);
    end
    """
    oc.eval(setup)

    def sweep():
        rows = []
        for N in (13, 21, 33, 49):
            oc.eval(f"N={N}; amp=0.10; k=2;")
            oc.eval("[X,Y,Z]=cn(N,amp); [xc,yc,zc]=cc(N,amp);")
            oc.eval("f=reshape(permute(xc.^2+yc.^2+zc.^2,[2 1 3]),[],1);")
            oc.eval("G=grad3DCurv(k,X,Y,Z); T=G*f;")
            oc.eval("s=linspace(0,1,N); cx=s(1:end-1)+0.5/(N-1);")
            oc.eval("[aa,bb,~]=meshgrid(s,cx,cx); xf=aa+amp*sin(2*pi*bb);")
            oc.eval("gx=T(1:N*(N-1)*(N-1)); "
                    "eg=sqrt(mean((gx-reshape(permute(2*xf,[2 1 3]),[],1)).^2));")
            oc.eval("L=div3DCurv(k,X,Y,Z)*G; v=reshape(L*f,N+1,N+1,N+1);")
            oc.eval("d=v(3:end-2,3:end-2,3:end-2)(:)-6; el=sqrt(mean(d.^2));")
            rows.append((N, float(oc.pull("eg")), float(oc.pull("el"))))
        return rows

    def show(rows, tag):
        print(f"\n  {tag}")
        print(f"  {'n':>5} {'grad rms':>13} {'order':>7} {'L=D*G rms':>13} {'order':>7}")
        pg = pl = pn = None
        for N, eg, el in rows:
            if pg is None:
                og = ol = "   -"
            else:
                r = np.log((N - 1) / (pn - 1))
                og = f"{np.log(pg/eg)/r:5.2f}"
                ol = f"{np.log(pl/el)/r:5.2f}"
            print(f"  {N:>5} {eg:>13.4e} {og:>7} {el:>13.4e} {ol:>7}")
            pg, pl, pn = eg, el, N

    show(sweep(), "STOCK GI13")
    oc.addpath(a.patch)          # addpath prepends -> patch now shadows MOLE
    print("\n  now using:", oc.feval("which", "GI13"))
    show(sweep(), "PATCHED GI13")

    print("\n  regression, Cartesian grid (must stay at roundoff):")
    for N in (11, 17, 25):
        oc.eval(f"N={N}; k=2;")
        oc.eval("[X,Y,Z]=meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));")
        oc.eval("s=linspace(0,1,N); s=[0, s(1:end-1)+0.5/(N-1), 1]; [x,y,z]=meshgrid(s,s,s);")
        oc.eval("f=reshape(permute(x.^2+y.^2+z.^2,[2 1 3]),[],1);")
        oc.eval("L=div3DCurv(k,X,Y,Z)*grad3DCurv(k,X,Y,Z); v=reshape(L*f,N+1,N+1,N+1);")
        oc.eval("d=v(3:end-2,3:end-2,3:end-2)(:)-6;")
        print(f"    n={N:>2}  rms {float(oc.pull('d').std()):.3e}")

    # ---------------------------------------------------------------- 3
    banner("3.  gridGen vs jacobian2D -- array-orientation mismatch")
    oc.addpath(os.path.join(a.grids, "seamount"))
    cwd = os.getcwd()
    os.chdir(a.grids)
    try:
        oc.eval("global BT BB; BT=2.0; BB=5.5;")
        oc.eval(f"cd('{a.grids}'); [Xg,Zg]=gridGen('TFI','seamount',99,38,false);")
        print("    gridGen returns  :", tuple(int(v) for v in oc.pull("Xg").shape),
              " (xi-nodes as ROWS)")
        print("    grad2DCurv wants : rows = eta, cols = xi   [n,m]=size(X)\n")
        for expr, tag in [("Xg,Zg", "as returned"), ("Xg.',Zg.'", "transposed ")]:
            oc.eval(f"J=jacobian2D(2,{expr});")
            J = oc.pull("J").ravel()
            s = "ALL +" if (J > 0).all() else ("ALL -" if (J < 0).all() else "MIXED")
            print(f"    {tag} : jacobian2D -> {s:<8}"
                  f"|J| {np.abs(J).min():.4g} .. {np.abs(J).max():.4g}")
    finally:
        os.chdir(cwd)

    # ---------------------------------------------------------------- 4
    banner("4.  TTM grid folding on the thesis Eq. 2.24 seamount")
    os.chdir(a.grids)
    try:
        print(f"    {'iters':>6} {'Jacobian':>18} {'min J':>13}")
        for it in (5, 10, 25, 50, 100, 300):
            oc.eval(f"cd('{a.grids}'); [X,Z]=gridGen('TTM','seamount',99,38,false,{it});")
            X = oc.pull("X"); Z = oc.pull("Z")
            xx = np.gradient(X, axis=0); xh = np.gradient(X, axis=1)
            zx = np.gradient(Z, axis=0); zh = np.gradient(Z, axis=1)
            J = xx * zh - xh * zx
            s = "ALL +" if (J > 0).all() else f"MIXED ({int((J<0).sum())} neg)"
            print(f"    {it:>6} {s:>18} {J.min():>13.4g}")
    finally:
        os.chdir(cwd)

    banner("done")


if __name__ == "__main__":
    main()
