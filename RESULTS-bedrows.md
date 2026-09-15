# Bed BC redesign — results

Replaces HANDOFF §7.3. The blocker is cleared.

Environment: Ubuntu 24.04, Octave 8.4.0 (same as the Linux box in §2),
MOLE cloned from `csrc-sdsu/mole` at session time, SuperLU (`splu`) — not
PARDISO, so wall times here are pessimistic.

---

## 1. What was wrong

Three different operators were in play and none of them matched:

| where | operator |
|---|---|
| LHS bed rows of `L` | `robinBC2D`: `alpha*phi + dphi/dn`, computational eta, spacing `dz` |
| RHS bed rows | physical contravariant flux `n.u*` divided by `h_eta` |
| `bc_residual()` | `bnd2`: `w - (z_xi/x_xi)*avg(u)` |

`--bsc` was one scalar trying to reconcile operators that differ point by
point on a curved bed. It cannot, which is the whole of §7.3: a residual floor
no `bsc` gets under, an accuracy lower bound and a stability upper bound, and
those bounds closing at Lx = 6 km.

The evidence was already in the old logs and did not need a new run: `div` at
1e-16 and `bed` at 6.8e-3 in the *same converged solve* means the equation
being solved is not the condition being measured.

## 2. What replaced it — `--bedmode constraint` (now the default)

`C` is built directly from the `bnd2` arrays that `bc_residual()` uses, so the
equation solved and the residual reported are one operator. Then

```
L[bedidx, :] = C*G          rhs[bedidx] = C*u* / dt
=>  C u = C(u* - dt*G*phi) = C u* - dt*(C u*/dt) = 0
```

exactly, on any grid, with nothing to tune. Left/right ghost rows keep their
Robin rows so `alpha` still removes the constant nullvector.

**Verified, not assumed:** `D*G` is exactly empty on the entire ghost ring —
0 nnz on all 162 boundary rows at 48x33, printed every run. Row replacement
discards nothing. Note `div2DCurv(k,X,Y)` with three arguments dispatches to
`div2DCurvLegacy`, so this is the legacy operator; the 5-arg version zeroes
boundary rows explicitly and the legacy one inherits it from
`divNonPeriodic`.

## 3. Acceptance test (§7.3) — PASSES

40 periods, Lx = 6 km, nx = 256, nz = 101, no `bsc`:

| omega/N | div | bed | 40 periods | angle | theory | err | fit rms |
|---|---|---|---|---|---|---|---|
| 0.2 | 1.7e-16 | 1.5e-15 | survives | 12.78 | 11.54 | +1.24 (+10.8%) | 101 m |
| 0.4 | 1.3e-16 | 1.6e-15 | survives | 27.16 | 23.58 | +3.58 (+15.2%) | 124 m |
| 0.6 | 1.6e-16 | 3.3e-15 | survives | 38.45 | 36.87 | +1.58 (+4.3%) | 9.6 m |
| 0.8 | 1.5e-16 | 8.5e-15 | survives | 60.81 | 53.13 | +7.68 (+14.4%) | 100 m |

Both residuals are **flat** over 40 periods — no doubling, no drift.
`max|u|` steady (0.6: 1.165e-02 at t/T=4 to 1.154e-02 at t/T=40).

A/B on the identical grid, only the bed BC differing:

```
--bedmode projection --bsc 12   DIVERGED t/T = 9.00, max|u| = 5.2e-01
--bedmode constraint            40 periods, bed 3.3e-15
```

t/T = 9.00 reproduces the handoff's 9.0 for omega/N = 0.6 exactly.

Controls still bracket:

```
--bedmode post    div 6.5e-04   bed 0.00e+00
--bedmode none    div 1.6e-16   bed 6.1e-01  (growing)
```

## 4. The equilibration is not a disguised `--bsc`

Constraint rows are scaled to match the median row-inf-norm of `L`, purely so
`splu` does not see 1/dx vs 1/dx^2 imbalance. Scaling a row and its rhs by the
same number cannot change the solution, and `--bedeq unit` proves it on the
real operators:

```
--bedeq auto  (x8.0e-02 .. 8.7e-02)   38.45 deg   max|u| 1.154e-02
--bedeq unit  (x1.0)                  38.45 deg   max|u| 1.154e-02
```

Identical to every digit printed. Synthetic test agrees to 2e-14.

## 5. Guards added, all of which fail loudly

* `C` is checked against the `bc_residual` loop on a random field at startup;
  `SystemExit` if it disagrees. This is the guard against §4 — `C` is built
  from flat indices with no reshape anywhere, and this proves it.
* `D*G` on the bed ghost rows is printed with a loud warning if non-empty.
* Bed-row/constraint pairing: count and distinctness, both raise.
* **Static projection test before any timestep**: projects a random field and
  reports both residuals, `SystemExit` if either fails. This is the §11 point
  — the whole design is validated in under a second, and no time stepping is
  needed to find out it is wrong.

## 6. What this did NOT fix

The angle errors are still +4% to +15% and this is a separate problem
(HANDOFF §8.2). Do not read the table above as improved physics:

* Three of the four have fit residuals near 100 m against the 4.9 m of §2.6.
  A residual that large means the tracked maxima are not on a line and the
  angle is not a measurement. Only 0.6 (9.6 m) is a real one.
* `max/mean` is 2.7-3.7 at 0.4/0.6/0.8, which by the criterion in the script's
  own output means reflections dominate. At 40 periods with 4.5 bounces across
  6 km the accumulation window is full of them. The auto window is set from
  the first bounce, but the field it is tracking is not.
* omega/N = 0.2 is essentially unchanged: 12.78 (+10.8%) over 40 periods here
  against 12.80 (+10.9%) over 20 in §2.6. Expected — 0.2 was already stable,
  so the constraint had nothing to fix there. That the number did not move is
  a consistency check, not a disappointment.

So §8.2 did not resolve itself once §8.1 was fixed. It is now the top item,
and it is cleanly separated from the stability question for the first time.

## 7. Suggested next steps

1. Look at the field (`view.py` on `iwb_r06_constraint.npz`) before anything
   else. §6 of the handoff, every time.
2. The reflection contamination is a measurement/domain question, not a
   solver one: either lengthen the domain, shorten the accumulation window, or
   both. Cheap to test now that runs do not diverge.
3. `cone3d` above amp~0.15 (§8.4) was hypothesised to share a root cause with
   the bed BC. That is now directly testable — port the constraint.
4. Monterey (§9) is no longer blocked on stability.

## 8. Files

* `iwbcurv.py` — patched. `--bedmode {constraint,projection,post,none}`,
  default `constraint`; `--bedeq {auto,unit}`. `oct2py` now imported inside
  `main()` so the module is importable without Octave.
* `bedrows.patch` — diff against the uploaded version.
* `test_bedrows.py` — the splice algebra without MOLE or Octave.
* `iwb_r06_constraint.npz` — omega/N = 0.6, 40 periods, for `view.py`.
