# Roadmap: from a verified linear solver to nearshore internal bores

Target: reproduce Figure 10 of Walter, Woodson, Arthur, Fringer & Monismith
(2012), *J. Geophys. Res.* 117 C07017 — canonical (ξ ≈ 0.2) against
non-canonical (ξ ≈ 2) bore behaviour on the southern Monterey Bay shelf, with
the same transect and stratification, using mimetic curvilinear operators.

## Where we are

**Verified, 2-D, linear, curvilinear** (this repository):

- second order against the exact seiche dispersion relation, with rotation
  (matched at every latitude from 0 to 85°, same error as without rotation)
- exact dynamical similarity between domains of different depth
- beam angle converging to linear theory at both frequencies tested
- the Monterey transect running with measured bathymetry and N(z), M2 forcing,
  transport-conserving barotropic forcing and rotation
- 5.2 ms/step at 315k cells on the GPU (cuDSS), fields matching the CPU to 1e-13

**Missing for the target: nonlinear advection.** A bore is a steepened front —
the paper reports Fr ≈ 0.7 and O(1 m) overturns. No linear solver produces one.

## What the 2021 work already established

From the dissertation (Valera 2021, CGU/SDSU) and the `reference2021/` archive:

- The **Fortran GCCOM cannot practically run this case**. §2.3.3 documents that
  the HPGF algorithm requires a tiny Δt or it produces artifacts, that this was
  traced to HPGF specifically, and that repair attempts failed; §2.4.3.1 records
  that Monterey and La Jolla canyon runs were attempted and abandoned for that
  reason. HPGF also forces sigma grids and pins each water column to one
  processor.
- The **mimetic GCCOM has a working nonlinear formulation**: gradient-based
  momentum (thesis eq. 4.3–4.7), mimetic upwind built from `interpol3D` and the
  one-sided `d2D`/`d3D` operators, buoyancy as a body force (no HPGF, hence no
  sigma-grid restriction). Lock exchange validated in 2-D and 3-D (Fr 0.705 at
  401×6×101, 0.29% from theory); seamount qualitative; **internal-wave beam
  failed** — beams form at the right angle then dissipate, in both 2-D and 3-D.

So the remaining obstacle is narrow and documented, and the reference
implementation survives.

## What is wrong with the 2021 code, specifically

Four defects found by reading `reference2021/`, each independently fixable:

1. **The time scheme is forward Euler.** `SSPRK101D.m`/`SSPRK102D.m` take `RHS`
   as a fixed array and never re-evaluate it. Put through symbolic algebra, the
   ten stages collapse to exactly `p + dt·RHS` — not approximately, identically.
   Every result from this code is first order in time, and the "order-of-magnitude
   wider time steps" came from dropping HPGF, not from SSP. Fix: pass a function
   handle and evaluate per stage.
2. **The one-sided operators are first order.** `d2D`/`d3D` build `[-1, 1]/dx`
   between adjacent nodes. First-order upwind diffusion erodes a narrow beam over
   twenty tidal periods, which is the thesis's own diagnosis of the failure.
   Fix: second- and fourth-order one-sided stencils. MOLE today has only
   `sidedNodal.m` (1-D, first order), so this is also the missing library piece.
3. **The upwind coefficients mix velocity components.** In `beam2Dv3.m` both
   momentum blocks set `a_pls = max(Ui,0)` and `a_neg = min(Uj,0)`, then use
   that pair for the ξ- AND the η-direction terms. Correct upwinding needs
   `Ui` for the ξ term and `Uj` for the η term. As written, vertical advection is
   weighted by the horizontal velocity. Two lines.
4. **The buoyancy coupling is one-way at the boundary.** `Main.m` builds
   `batv(i,j,k) = (buoy(i,j,k) + buoy(i,j,k-1))/2` over `k = 2:o+1`, so the bed
   face never receives buoyancy, while the first interior cell is advected by a
   `w` that includes it. Centre→face by plain averaging, face→centre by upwind
   gradient: not adjoint in any inner product. This is the same class of defect
   found independently in this repository's linear solver, where it drove a
   numerical instability over sloping topography — invisible on a flat bottom,
   active on any slope. Fix: the area-weighted adjoint, already verified here.

Separately: the October 2021 curvilinear seamount used `grad3DCurv`, which calls
`GI13` — broken in MOLE until csrc-sdsu/mole#466. Those results were computed
with a gradient that does not converge on curvilinear grids.

## Stage 1 — nonlinear advection in 2-D — DONE

Implemented in `iwbcurv.py` (`--advect scalar|full`, `--advscheme
upwind1|upwind2|minmod`, `--advbg split|total`) and verified:

| check | result |
|---|---|
| flux form reproduces `div(u)` for a constant field | exact, every scheme |
| scalar order on a smooth bump | 0.5 / 1.3-1.6 / 1.2-1.5 |
| monotonicity at a step | bounded for upwind1 and minmod; upwind2 overshoots 0.21 |
| momentum: uniform flow, `u = x`, solid-body rotation | 1e-15 |
| momentum on a smoothly stretched grid, refined | order 1.83 -> 1.97 |
| lock release, 101 -> 801 cells (minmod) | Fr 0.654, 0.663, 0.674, 0.686 |
| lock release at 401x101, upwind2 | **Fr 0.699** vs theory 0.7071, published 0.705 |
| nonlinear beam, 25 periods, 512x201 | stable, 53.55 deg vs 53.13, max/mean 13.8 |
| regression, `--advect none` | bit-identical to every earlier result |

Findings worth carrying forward:

- **First-order advection alone does not explain the 2021 beam failure.** With
  `upwind1` -- the accuracy of `d2D`/`d3D` -- beams survive 25 periods at
  53.21 deg. The structural difference is what the scheme is applied TO: this
  solver advects only the perturbation and adds the background analytically as
  `N^2 w`, while the 2021 code advects the full temperature field, so its
  numerical diffusion acts on the stratification itself. `--advbg total`
  diverges at t/T = 1.8 where `split` is stable, which points the same way --
  though that implementation does not yet treat the background flux through a
  sloping bed consistently, so it is a lead rather than a result.
- **minmod reverts to first order at a front.** In the lock release it landed
  exactly on `upwind1` (0.6742 against 0.6741) while `upwind2` reached 0.699.
  Monotonicity is bought at the price of accuracy precisely where a bore lives,
  so the bore case may need a higher-order limiter.
- **More diffusion can flatter the diagnostics.** `upwind1` gave the straightest
  beam and the smallest angle bias, because it smooths the far field where the
  sag lives. Not accuracy.
- **`alpha` is scale-dependent, not a constant.** It regularises a matrix whose
  entries scale as 1/dx^2: 1e-6 is right at 6 km and breaks the projection at
  0.8 m, where ~1 gives residuals of 1e-14. A dimensionless form, normalised by
  the mean Laplacian diagonal, would remove the trap.

Scripts: `advect_test.py` (operator tests), `lock_test.ps1` (benchmark),
`beam_nonlinear.ps1` (scheme comparison on the beam).

## Stage 1 as originally planned

Port the `Main.m` structure into `iwbcurv.py`, carrying over what is verified
here (projection, bed constraint, energy-consistent buoyancy, sponges, GPU path)
and fixing the four defects above. Verification ladder, in order:

1. **Linear limit.** At small amplitude the nonlinear code must reproduce the
   present linear solver. Every result in this repository becomes a test.
2. **Passive scalar in a uniform flow.** Exact; isolates scalar advection and
   its monotonicity.
3. **Lock release** — Garcia §3.1–3.2, with published mimetic-GCCOM numbers to
   match: Fr 0.705 (3-D no-slip), 0.6882/0.7002 (2-D), theory 0.7071.
4. **Internal-wave beam, nonlinearly.** The case the 2021 model failed and this
   one already passes linearly. Beams must survive twenty periods at the correct
   angle.
5. **Energy budget.** Free runs and growth probes, as used for the buoyancy fix.

Design choices to settle: flux form versus vector-invariant momentum; limiter
versus centred scalar advection (a bore is exactly where this shows, so make it
a flag); explicit viscosity and diffusivity (`--nu`, `--kappa`) both to match
Walter's 1e-4 m²/s and as a fallback.

## Stage 2 — mode-1 forcing

Impose the mode-1 internal wave at the offshore boundary using the vertical
structure from `mode1.py`. For this stratification: **c₁ = 0.209 m/s, λ = 9.4 km**,
and ξ ≈ 2 on the real slope (s ≈ 0.04) needs **a ≈ 3.7 m**. ξ ≈ 0.2 would need
370 m, so the canonical case requires a gentler synthetic slope, as in the paper.
Verify by propagating the mode across a flat domain at the computed speed.

## Stage 3 — the comparison

Extend the transect offshore to 20 km, run both slopes for six tidal periods,
and reproduce: the cross-shelf temperature snapshots, virtual thermistors at
2/4/6 m above bed at the 15 m isobath, bore speed against the observed 0.12 m/s
and Fr ≈ 0.7, and the qualitative split (bolus at low ξ, sloshing at high ξ).

Their time step was 1 s, 2.7e5 steps. Ours is limited by advective CFL
(dt < 33 s), buoyancy stiffness (N_max·dt < 1 → dt < 65 s) and front resolution.
At dt ≈ 10–25 s that is 10k–27k steps: minutes on the GPU. If it holds, running
a published nonhydrostatic case orders of magnitude faster is a result in itself.

## Stage 4 — three dimensions

A slab first — a few cells across, 3-D operators, periodic in the third
direction, reproducing the 2-D answer — which exercises `GI13`, `jacobian3D` and
the 3-D projection on a problem whose answer is known. Both of those paths were
broken in stock MOLE until #466 and #468, so a slab built earlier would have been
quietly wrong. Then real bathymetry, at which point the direct solver stops
fitting and PETSc with multigrid — the framework built in Valera et al. (2019) —
becomes the right tool.

## Still open, for Jared

The archive answers the scheme questions; two things remain:

- does the current GCCOM tree still build and run, and does the Bitbucket
  repository the thesis references still exist
- do the Monterey transect and forcing files survive anywhere (nothing
  Monterey-related is in `reference2021/`)

Worth telling him either way: the buoyancy pairing defect found here is present
in the 2021 mimetic code too, in a different form, and any code that pairs
interpolators this way over topography inherits it.
