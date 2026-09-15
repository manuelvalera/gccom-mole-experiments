#!/usr/bin/env python3
"""
pardiso_tune.py -- find the best MKL thread count for the pressure solve.

MKL reads its thread settings at import time, so each configuration has to run
in a fresh process.  This script re-launches itself as a subprocess per setting.

USAGE
    python pardiso_tune.py --mole C:/path/to/mole/src/matlab_octave --a 48

On a hybrid Intel CPU (12th gen and later: P-cores + E-cores) the best value is
usually the number of PHYSICAL P-CORES, not the total thread count.  MKL splits
work statically, so E-cores become stragglers the whole team waits on.
i7-14700F: 8 P-cores (16 threads) + 12 E-cores = 20 cores / 28 threads -> try 8.
"""
import argparse, os, subprocess, sys, time


def child(mole, a, nsolve):
    import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
    from oct2py import Oct2Py
    oc = Oct2Py(); oc.eval("warning('off','all');"); oc.addpath(mole)

    def get(nm):
        oc.eval(f"[ii,jj,vv]=find({nm}); sz=size({nm});")
        return sp.csr_matrix(
            (oc.pull('vv').ravel(),
             (oc.pull('ii').ravel().astype(int) - 1,
              oc.pull('jj').ravel().astype(int) - 1)),
            shape=tuple(oc.pull('sz').ravel().astype(int)))

    oc.eval(f"k=2; a={a}; h=1/a; G=grad3D(k,a,h,a,h,a,h); "
            f"D=div3D(k,a,h,a,h,a,h); B=robinBC3D(k,a,h,a,h,a,h,0,1);")
    G, D, B = get('G'), get('D'), get('B'); oc.exit()
    nc = (a + 2) ** 3
    L = (D @ G + B).tocsc()
    bb = np.random.default_rng(7).standard_normal(nc)
    u = bb - L @ spl.lsqr(L, bb, atol=1e-10, btol=1e-10, iter_lim=3000)[0]
    ip = int(np.argmax(np.abs(u)))
    L = L.tolil(); L[ip, :] = 0; L[ip, ip] = 1; L = L.tocsr()

    import pypardiso
    ps = pypardiso.PyPardisoSolver()
    t = time.time(); ps.factorize(L); tf = time.time() - t
    rhs = np.random.default_rng(1).standard_normal(nc); rhs[ip] = 0
    ps.solve(L, rhs)                       # warm up
    t = time.time()
    for _ in range(nsolve):
        ps.solve(L, rhs)
    ts = (time.time() - t) / nsolve * 1000
    try:
        import psutil
        rss = psutil.Process().memory_info().rss / 1e9
    except Exception:
        rss = float('nan')
    print(f"RESULT {os.environ.get('MKL_NUM_THREADS','?')} {tf:.2f} {ts:.2f} {rss:.2f}",
          flush=True)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--mole', required=True)
    ap.add_argument('--a', type=int, default=48)
    ap.add_argument('--nsolve', type=int, default=15)
    ap.add_argument('--threads', type=int, nargs='+',
                    default=[1, 2, 4, 6, 8, 12, 16, 20, 28])
    ap.add_argument('--_child', action='store_true')
    g = ap.parse_args()

    if g._child:
        child(g.mole, g.a, g.nsolve)
        sys.exit(0)

    print(f"\nPARDISO thread sweep, {g.a}^3   ({(g.a+2)**3} unknowns)")
    print(f"{'threads':>8} {'factorise s':>12} {'solve ms':>10} {'peak RSS GB':>12} {'speedup':>9}")
    base = None
    for t in g.threads:
        env = dict(os.environ)
        env['MKL_NUM_THREADS'] = str(t)
        env['OMP_NUM_THREADS'] = str(t)
        env['MKL_DYNAMIC'] = 'FALSE'
        cmd = [sys.executable, os.path.abspath(__file__), '--_child',
               '--mole', g.mole, '--a', str(g.a), '--nsolve', str(g.nsolve)]
        try:
            out = subprocess.run(cmd, env=env, capture_output=True, text=True,
                                 timeout=3600).stdout
            line = [l for l in out.splitlines() if l.startswith('RESULT')]
            if not line:
                print(f"{t:>8}   (no result -- likely out of memory)")
                continue
            _, th, tf, ts, rss = line[0].split()
            tf, ts, rss = float(tf), float(ts), float(rss)
            if base is None:
                base = ts
            print(f"{t:>8} {tf:>12.2f} {ts:>10.2f} {rss:>12.2f} {base/ts:>8.2f}x")
        except subprocess.TimeoutExpired:
            print(f"{t:>8}   (timed out)")
    print("\nUse the thread count with the lowest solve time -- that is 96% of "
          "each timestep.\nSet it before launching:  $env:MKL_NUM_THREADS=\"<n>\"; "
          "$env:OMP_NUM_THREADS=\"<n>\"\n")
