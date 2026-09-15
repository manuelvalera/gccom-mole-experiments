"""Verify the row-replacement projection identity without MOLE or Octave.

Builds a synthetic problem with the same STRUCTURE as iwbcurv:
  - D : nc x (nu+nw), zero on the ghost ring (MOLE's convention)
  - G : (nu+nw) x nc
  - B : boundary rows only (the Robin stand-in)
  - C : nbed x (nu+nw), one row per bed face, w - rat*avg(u) = 0
and checks that after the splice the projected field satisfies BOTH
  D (f - G phi) = 0   and   C (f - G phi) = 0
to machine precision, with no tuned coefficient.

The point of the test is the plumbing I wrote -- index pairing, row zeroing,
the Pb splice, and the claim that the equilibration cancels.  It says nothing
about MOLE's operators, which is the part that needs the real run.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from iwbcurv import splice_bed_rows

rng = np.random.default_rng(0)
m, n = 9, 7                                   # xi, eta cells
nu_, nw_ = (m + 1) * n, m * (n + 1)
nc = (m + 2) * (n + 2)

iu2 = np.arange(nu_).reshape(n, m + 1)        # (eta, xi), xi fastest
iw2 = np.arange(nw_).reshape(n + 1, m)
ic2 = np.arange(nc).reshape(n + 2, m + 2)

# --- synthetic D with an empty ghost ring, and a dense-ish G ----------------
interior = ic2[1:-1, 1:-1].ravel()
rows, cols, vals = [], [], []
for k, c in enumerate(interior):
    j = rng.choice(nu_ + nw_, 6, replace=False)
    rows += [c] * 6; cols += list(j); vals += list(rng.standard_normal(6) / 1e-2)
D = sp.csr_matrix((vals, (rows, cols)), shape=(nc, nu_ + nw_))

G = sp.random(nu_ + nw_, nc, density=0.05, random_state=1, format='csr')
G.data *= 10.0

ghost = np.setdiff1d(np.arange(nc), interior)
B = sp.csr_matrix((rng.standard_normal(ghost.size) + 5.0,
                   (ghost, ghost)), shape=(nc, nc))

# --- C, built from bnd2-shaped index arrays ---------------------------------
bnd2 = []
for jp_, jc in ((0, 0), (n, n - 1)):
    rat = rng.standard_normal(m)
    bnd2.append((iw2[jp_, :], iu2[jc, :-1], iu2[jc, 1:], rat))

_r, _rows, _cols, _vals = 0, [], [], []
for tgt, s1, s2, rat in bnd2:
    rr_ = np.arange(_r, _r + tgt.size)
    _rows += [rr_, rr_, rr_]
    _cols += [nu_ + tgt, s1, s2]
    _vals += [np.ones(tgt.size), -0.5 * rat, -0.5 * rat]
    _r += tgt.size
nbed = _r
C = sp.csr_matrix((np.concatenate(_vals),
                   (np.concatenate(_rows), np.concatenate(_cols))),
                  shape=(nbed, nu_ + nw_))

# C must reproduce the bc_residual loop exactly (the guard in iwbcurv)
uu, ww = rng.standard_normal(nu_), rng.standard_normal(nw_)
ref = np.concatenate([ww[t] - r * 0.5 * (uu[a] + uu[b]) for t, a, b, r in bnd2])
err = np.abs(ref - C @ np.concatenate([uu, ww])).max()
print(f"C vs bnd2 loop                    {err:.2e}")
assert err < 1e-13

bedidx = np.concatenate([ic2[0, 1:-1], ic2[n + 1, 1:-1]])
assert bedidx.size == nbed
print(f"D*G on bed ghost rows             {np.abs((D @ G)[bedidx]).max() if (D @ G)[bedidx].nnz else 0.0:.2e}")

L0 = (D @ G + B).tocsr()
f = rng.standard_normal(nu_ + nw_)

# --- baseline: no bed rows.  div is enforced, bed is not. -------------------
lu0 = spl.splu(L0.tocsc())
p0 = f - G @ lu0.solve(D @ f)
print(f"\nno constraint : div {np.linalg.norm(D @ p0)/np.linalg.norm(D @ f):.2e}"
      f"   bed {np.abs(C @ p0).max()/np.abs(f).max():.2e}")

# --- spliced, both equilibrations -------------------------------------------
sols = {}
for eq in ('auto', 'unit'):
    L, bedsc = splice_bed_rows(L0, C @ G, bedidx, eq)
    rhs = D @ f
    rhs[bedidx] = bedsc * (C @ f)
    phi = spl.splu(L.tocsc()).solve(rhs)
    pf = f - G @ phi
    dv = np.linalg.norm(D @ pf) / np.linalg.norm(D @ f)
    bd = np.abs(C @ pf).max() / np.abs(f).max()
    sols[eq] = phi
    print(f"constraint {eq:5s}: div {dv:.2e}   bed {bd:.2e}"
          f"   scale range x{bedsc.min():.1e}..{bedsc.max():.1e}")
    assert dv < 1e-8 and bd < 1e-8, f"{eq} failed"

d = np.abs(sols['auto'] - sols['unit']).max() / np.abs(sols['unit']).max()
print(f"\nauto vs unit solution difference  {d:.2e}   "
      f"{'INVARIANT (not a tuned coefficient)' if d < 1e-8 else '*** NOT INVARIANT ***'}")
assert d < 1e-8

# --- pairing guards fire ----------------------------------------------------
for bad, why in ((bedidx[:-1], 'wrong count'),
                 (np.r_[bedidx[:-1], bedidx[0]], 'duplicate index')):
    try:
        splice_bed_rows(L0, C @ G, bad)
        print(f"guard MISSED: {why}")
    except ValueError as e:
        print(f"guard fires on {why:16s}: {e}")

print("\nall checks passed")
