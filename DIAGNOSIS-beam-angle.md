# Beam-angle bias — resolved. HPGF is not needed.

Supersedes the earlier version of this file, whose §3-§5 were **wrong** and
are retracted in §3 below.

Bottom line: the residual beam-angle error is ordinary discretization error.
It converges at roughly second order under joint refinement, `r_x` is not
involved, and there is no missing physics. HPGF should not be built for this
solver.

---

## 1. Step 1 — the contamination was lateral, and the sponge fixes it

`Lsl = LX/10` and `TAUS = 100` were hardcoded. Exposed as `--lsl` / `--taus`.
At the worst-affected frequency (omega/N = 0.8, 20 periods):

| lsl | taus | fit rms | error |
|---|---|---|---|
| 0.10 | 100 | 92.2 m | +7.95 |
| 0.20 | 100 | 26.5 m | +2.60 |
| 0.20 | 50 | 5.7 m | +1.61 |
| 0.30 | 25 | 5.2 m | +1.02 |

Energy was returning from the lateral boundaries because the sponge was too
narrow and too weak. Nothing to do with the solver.

`omega/N = 0.2` additionally needed a bigger box: its bounce is 4900 m, so at
Lx = 6 km the fit window ended 300 m from the wall, inside the sponge itself.
At Lx = 12 km the error dropped from +0.88 to +0.59.

**Protocol adopted**: `lsl = 0.2`, `taus = T/30`, `spp = 120` (keeps
`dt/taus = 0.25`), `Lx` sized so the window clears the sponge.

### A silent stability limit, now guarded

The sponge relaxation is EXPLICIT, so `dt/taus` is a stability number:

```
dt/taus = 1.0  runs
dt/taus = 2.5  max|u| = 1.0e+00     blow-up
dt/taus = 5.0  max|u| = 5.1e+08     blow-up
```

Same mechanism as HANDOFF §4 item 4, and there was no guard — it looked like a
physics divergence. Now printed every run and refused above 1.0, with the
required `--taus` or `--spp` in the message.

## 2. Result of step 1 — all four frequencies became trustworthy

| omega/N | fit rms before | after | max/mean before | after |
|---|---|---|---|---|
| 0.2 | 3.7 m | 4.7 m | 10.3 | 18.7 |
| 0.4 | 90.8 m | 4.5 m | — | 13.6 |
| 0.6 | 5.4 m | 3.5 m | 3.9 | 12.3 |
| 0.8 | 92.2 m | 6.4 m | — | 10.9 |

Plotting the field was not necessary to arbitrate anything after this.

## 3. RETRACTION — `r_x` is not the driver

The earlier version claimed the bias correlated with `r_x` across two knobs.
Both knobs (eta stretching, and `nx`) also changed resolution, so the
correlation was confounded. The control was never run. It has now been run:
vary `r_x` by a factor of ten via ridge height, at **fixed `dz` = 10 m**:

| ab | r_x | error |
|---|---|---|
| 20 m | 0.709 | +1.84 |
| 10 m | 0.353 | +1.99 |
| 5 m | 0.176 | +2.16 |
| 2 m | 0.070 | +1.92 |

Flat. `r_x` moves by 10x and the answer does not move. The earlier §3, and the
inference in §4-§5 that followed from it, are withdrawn.

## 4. What it actually is — discretization error, converging

Refining both directions together (ab = 10 so `r_x` never binds, omega/N = 0.6):

| nx | nz | dx | dz | error |
|---|---|---|---|---|
| 256 | 101 | 23.5 m | 10.0 m | +1.99 |
| 384 | 151 | 15.7 m | 6.7 m | +1.06 |
| 512 | 201 | 11.7 m | 5.0 m | +0.37 |

Roughly second order. The single-direction sweeps looked first-order only
because the unrefined direction's error dominated.

Confirmed at all four frequencies on the benchmark ridge (ab = 20), halving
`dz` alone:

| omega/N | theory | nz=101 | nz=201 | fit rms | max/mean |
|---|---|---|---|---|---|
| 0.2 | 11.54 | +0.59 | **+0.32** | 2.2 m | 33.7 |
| 0.4 | 23.58 | +1.22 | **+0.50** | 1.8 m | 20.0 |
| 0.6 | 36.87 | +1.84 | **+0.71** | 1.3 m | 14.8 |
| 0.8 | 53.13 | +1.43 | **+0.57** | 2.7 m | 11.5 |

Roughly halves everywhere. This also retires the "effective N is 5% low"
reading: as `sin(phi)` excess it went 5.0/4.9/4.2/1.8% at nz=101 to
2.8/2.0/1.6/0.7% at nz=201. It scales with resolution, so it was never a
fixed offset.

`div` and `bed` stay at 1e-16 / 1e-15 throughout, including at nz=201 and
`r_x = 1.42`.

## 5. Why HPGF does not apply here

`r_x` is Shchepetkin & McWilliams' criterion for the **hydrostatic**
pressure-gradient force in sigma-coordinate models, where the horizontal PGF
is computed as a difference of two large near-cancelling terms along sigma
surfaces. HPGF (the density-Jacobian formulation) fixes that specific
cancellation.

`iwbcurv` is **non-hydrostatic**. Pressure comes from the full Poisson
projection; buoyancy enters as a body force `dt*(Ic_z @ b)` on the w-faces.
There is no hydrostatic PGF in the algorithm, so there is nothing for HPGF to
reformulate.

Measured directly with the new `--forcegrid`, benchmark ridge, omega/N = 0.6:

| nz | r_x | fit rms | max/mean | error |
|---|---|---|---|---|
| 101 | 0.709 | 3.5 m | 12.3 | +1.84 |
| 151 | 1.065 **VIOLATED** | 1.9 m | 14.0 | +1.19 |
| 201 | 1.422 **VIOLATED** | 1.3 m | 14.8 | +0.71 |

Violating the criterion by 42% makes the answer strictly **better** on every
measure. Nothing degrades.

**Consequence for Monterey (§9).** The handoff expects `r_x > 1` on a canyon
flank to be a blocker requiring HPGF. On this evidence it is not a blocker for
this solver, and the grid space §9 assumed was closed is open. The `r_x`
abort should become a warning. The folded-Jacobian check should stay an abort
— that one is real.

## 6. Open, in priority order

1. **Re-validate against Garcia et al. Fig 8/9 amplitudes.** §8.2's other
   half — beams ~3x weaker than theirs — has not been retested since the
   sponge changed, and the sponge directly sets amplitude. Do this first.
2. **Why `ab = 20` does not beat `ab = 10` at fixed `dz`** (+1.84 vs +1.99).
   A taller ridge should be better resolved relative to the beam. Small, but
   it is the one number that does not fit the convergence picture.
3. **Turn the `r_x` abort into a warning**, keep the Jacobian abort.
4. **Monterey.** No longer blocked on stability or on `r_x`.
5. **`cone3d` above amp~0.15** (§8.4) — port `--bedmode constraint`.

## 7. Do not re-test

Measurement (tracker verified to ±0.08 deg against a synthetic beam on the
real grid), timestep (angle converged by spp=60), operator order (k=2/4/6
within 0.25 deg), curvilinear metrics (bias anti-correlated with distortion),
`r_x` (§3). All eliminated with runs on record for **this** solver, not
inherited from `cone3d`.

## 8. New flags

```
--lsl        sponge width as a fraction of Lx      (was hardcoded Lx/10)
--taus       sponge relaxation time in s           (was hardcoded 100)
--order      MOLE operator order k, 2/4/6          (was hardcoded 2)
--ab --lb    ridge height / half-width in m        (were module constants)
--forcegrid  run despite a failed r_x check
--bedmode    constraint (default) / projection / post / none
--bedeq      auto / unit
```

Defaults reproduce the pre-session behaviour exactly.
