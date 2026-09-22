#!/usr/bin/env python3
"""solver_bench.py -- which pressure solver is fastest on THIS matrix, and is
it still exact?

  python iwbcurv.py ... --save-matrix L896.npz --nper 1     (dump it once)
  python solver_bench.py L896.npz
  python solver_bench.py L896.npz --backends superlu pardiso cupy64 cupy32r

The model factors the pressure matrix L once and then solves with it every
time step. This times exactly that -- factor once, solve many -- for each
backend, on the real matrix, and checks every answer against a known solution
(b = L @ x_true), because a solver that is fast but loses the projection's
~1e-15 accuracy is not usable here.

Backends (each is skipped with a reason if its library is missing):

  superlu   SciPy SuperLU, single-threaded CPU        (the model's default)
  pardiso   MKL PARDISO via pypardiso, threads = MKL_NUM_THREADS
  cupy64    CuPy splu: factor on the CPU, triangular solves on the GPU, FP64
  cupy32r   the same in FP32 on the GPU, corrected to FP64 accuracy by
            iterative refinement with an FP64 residual. FP32 is the GeForce
            card's strong suit; whether refinement converges depends on how
            well conditioned L is, and the script reports that rather than
            assuming it.
  cudss64   NVIDIA cuDSS via nvmath-python (full GPU factor + solve), FP64.
            Written against the documented nvmath API but NOT tested here;
            if it fails, the error is printed so it can be fixed.

GPU SAFETY. This is a light, memory-bandwidth-bound load, nothing like a game
or a training run. Even so: the script reads temperature and power from
nvidia-smi before and after each GPU backend, refuses to start one if the card
is already above 80 C, and each timing loop is a few seconds long.
"""
import argparse, os, subprocess, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def gpu_status():
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=temperature.gpu,power.draw,power.limit,"
             "utilization.gpu,memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10).stdout.strip().splitlines()[0]
        t, p, pl, u, mu, mt = [v.strip() for v in out.split(",")]
        return float(t), f"{t} C  {float(p):.0f}/{float(pl):.0f} W  util {u}%  mem {mu}/{mt} MiB"
    except Exception as e:                                  # no nvidia-smi, no GPU
        return None, f"(nvidia-smi unavailable: {e.__class__.__name__})"


def report(name, tf, ts, x, x_true, A, b, note=""):
    res = np.linalg.norm(A @ x - b) / np.linalg.norm(b)
    err = np.linalg.norm(x - x_true) / np.linalg.norm(x_true)
    # both must hold: a small residual alone is not enough on an ill-conditioned
    # matrix (FP32 refinement reached 1.4e-12 residual with a 7e-4 error)
    ok = "OK " if (res < 1e-11 and err < 1e-6) else "BAD"
    print(f"  {name:9s} factor {tf:7.2f} s   solve {1e3*ts:8.2f} ms   "
          f"residual {res:.1e} {ok}  error {err:.1e}  {note}")
    return ts


# ---------------------------------------------------------------- backends
def bench_superlu(A, b, rep):
    Ac = A.tocsc()
    t = time.perf_counter(); lu = spl.splu(Ac); tf = time.perf_counter() - t
    lu.solve(b)
    t = time.perf_counter()
    for _ in range(rep):
        x = lu.solve(b)
    return tf, (time.perf_counter() - t) / rep, x, ""


def bench_pardiso(A, b, rep):
    from pypardiso import PyPardisoSolver
    Ar = sp.csr_matrix(A); Ar.sort_indices()
    ps = PyPardisoSolver(mtype=11)
    t = time.perf_counter(); ps.factorize(Ar); tf = time.perf_counter() - t
    def solve(rhs):
        ps.set_phase(33); return ps._call_pardiso(Ar, rhs)
    solve(b)
    t = time.perf_counter()
    for _ in range(rep):
        x = solve(b)
    return tf, (time.perf_counter() - t) / rep, x, f"threads={os.environ.get('MKL_NUM_THREADS', 'default')}"


def _cupy():
    import cupy as cp
    import cupyx.scipy.sparse as csp
    import cupyx.scipy.sparse.linalg as cspl
    return cp, csp, cspl


def bench_cupy64(A, b, rep):
    cp, csp, cspl = _cupy()
    Ag = csp.csr_matrix(sp.csr_matrix(A, dtype=np.float64))
    t = time.perf_counter(); lu = cspl.splu(Ag); cp.cuda.Device().synchronize()
    tf = time.perf_counter() - t
    x = cp.asnumpy(lu.solve(cp.asarray(b)))                  # warm-up / analysis
    t = time.perf_counter()
    for _ in range(rep):                                     # includes both PCIe transfers
        x = cp.asnumpy(lu.solve(cp.asarray(b)))
    cp.cuda.Device().synchronize()
    return tf, (time.perf_counter() - t) / rep, x, \
        f"GPU mem {cp.get_default_memory_pool().used_bytes()/2**20:.0f} MiB"


def bench_cupy32r(A, b, rep, tol=1e-13, maxit=8):
    cp, csp, cspl = _cupy()
    A64 = csp.csr_matrix(sp.csr_matrix(A, dtype=np.float64))
    A32 = csp.csr_matrix(sp.csr_matrix(A, dtype=np.float32))
    t = time.perf_counter(); lu = cspl.splu(A32); cp.cuda.Device().synchronize()
    tf = time.perf_counter() - t
    bg = cp.asarray(b); nb = float(cp.linalg.norm(bg))

    def solve(bh):
        bg_ = cp.asarray(bh)
        x = lu.solve(bg_.astype(cp.float32)).astype(cp.float64)
        its = 0
        for its in range(1, maxit + 1):
            r = bg_ - A64 @ x                                # FP64 residual on the GPU
            if float(cp.linalg.norm(r)) / nb < tol:
                break
            x += lu.solve(r.astype(cp.float32)).astype(cp.float64)
        return cp.asnumpy(x), its

    x, its = solve(b)
    t = time.perf_counter()
    for _ in range(rep):
        x, its = solve(b)
    cp.cuda.Device().synchronize()
    return tf, (time.perf_counter() - t) / rep, x, f"refinement steps {its}"


def bench_cudss64(A, b, rep):
    import nvmath
    from nvmath.sparse.advanced import DirectSolver
    Ar = sp.csr_matrix(A, dtype=np.float64); Ar.sort_indices()
    b2 = b.reshape(-1, 1).copy()
    t = time.perf_counter()
    solver = DirectSolver(Ar, b2)
    solver.plan(); solver.factorize()
    tf = time.perf_counter() - t
    x = np.asarray(solver.solve()).ravel()
    t = time.perf_counter()
    for _ in range(rep):
        solver.reset_operands(b=b2)
        x = np.asarray(solver.solve()).ravel()
    tsol = (time.perf_counter() - t) / rep
    solver.free()
    return tf, tsol, x, ""


BACKENDS = {"superlu": bench_superlu, "pardiso": bench_pardiso,
            "cupy64": bench_cupy64, "cupy32r": bench_cupy32r, "cudss64": bench_cudss64}
GPU = {"cupy64", "cupy32r", "cudss64"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--backends", nargs="+", default=list(BACKENDS))
    ap.add_argument("--repeat", type=int, default=20)
    g = ap.parse_args()

    A = sp.load_npz(g.matrix).tocsr()
    rng = np.random.default_rng(0)
    x_true = rng.standard_normal(A.shape[0])
    b = A @ x_true
    print(f"\n{g.matrix}: n = {A.shape[0]}, nnz = {A.nnz}, {g.repeat} solves per backend\n")

    times = {}
    for name in g.backends:
        if name in GPU:
            temp, txt = gpu_status()
            print(f"  [gpu before {name}: {txt}]")
            if temp is not None and temp > 80:
                print(f"  {name:9s} SKIPPED: GPU already at {temp:.0f} C"); continue
        try:
            tf, ts, x, note = BACKENDS[name](A, b, g.repeat)
            times[name] = report(name, tf, ts, x, x_true, A, b, note)
        except ImportError as e:
            print(f"  {name:9s} not available ({e})")
        except Exception as e:
            print(f"  {name:9s} FAILED: {e.__class__.__name__}: {str(e)[:200]}")
        if name in GPU:
            print(f"  [gpu after  {name}: {gpu_status()[1]}]")

    if "superlu" in times:
        print("\n  speed-up of the per-step solve relative to SuperLU:")
        for k, v in times.items():
            print(f"     {k:9s} {times['superlu']/v:5.2f}x")
    print("\n  Only backends marked OK keep the projection exact enough to use.")


if __name__ == "__main__":
    main()
