# HANDOFF v3 — iwbcurv / GCCOM-MOLE revival

Written to be self-contained: a fresh assistant with no prior context should be
able to work from this alone. Supersedes HANDOFF v2 and the original HANDOFF.

---

## 0. What this is

`iwbcurv.py` is a 2D, **linear**, non-hydrostatic, Boussinesq internal-wave
solver on a fully curvilinear grid, built on the MOLE mimetic operator library
(Castillo–Grone). It is a testbed for reviving the mimetic GCCOM coastal ocean
model. The validation target is Garcia et al. (2019), *J. Comput. Sci.*
30:143–156, §3.3 — internal wave beams radiated by an oscillating tidal flow
over a Gaussian ridge.

**The solver has no advection term.** `us = u − dt·gp − sponge`. Everything is
linear in the forcing amplitude; a 4× change in `--u0` gives bit-identical
output. Do not propose nonlinear explanations for anything.

Algorithm: incremental projection. Predict `u*` with lagged pressure and
buoyancy, solve `L φ = D u*/dt`, correct `u = u* − dt G φ`, then step buoyancy
`b −= dt N² (Idf_w w)`. Time scheme is **second order** (verified on the
seiche, below). MOLE supplies `grad2DCurv`, `div2DCurv` (note: the 3-argument
call dispatches to `div2DCurvLegacy`), `robinBC2D`, `interpol2D`,
`interpolD2D`, `gridGen`.

## 1. Environment

MOLE commit `b0b9ff7` (2026-09-03), Octave 8.4.0, Python 3.11–3.12, numpy 2.4,
scipy 1.17, SuperLU via `scipy.sparse.linalg.splu`. Reproduced independently on
Linux and Windows/anaconda with matching numbers, including `DIVERGED at
t/T = 9.00` for the pre-fix bed penalty.

`--mole` must point at MOLE's `src/matlab_octave` (the directory containing
`gridGen.m`). Octave's `addpath` does **not** error on a missing directory; the
solver now checks for `gridGen.m` up front.

## 2. What is fixed and validated

### The bed boundary condition — DONE

The old `--bedmode projection` was a penalty with three mutually inconsistent
operators: `robinBC2D`'s computational-eta normal derivative on the LHS, the
physical contravariant flux over `h_eta` on the RHS, and `bnd2`'s
`w = (z_xi/x_xi)·avg(u)` in `bc_residual()`. The scalar `--bsc` was trying to
reconcile operators that differ point by point on a curved bed. It cannot.

`--bedmode constraint` (now default) builds `C` from the same `bnd2` arrays the
residual uses. `D*G` is empty on the ghost ring (0 nnz on all 162 boundary rows,
verified), so the bottom/top ghost rows carry `C*G` with rhs `C u*/dt`, giving
`C u = C(u* − dt G φ) = 0` exactly, with no coefficient.

Result: `div ~1e-16`, `bed ~1e-15`, flat over 40 periods, all four frequencies.
The old penalty at `bsc = 12` diverges at `t/T = 9.00` on the identical grid.
Row equilibration on the constraint rows is bookkeeping, not tuning:
`--bedeq auto` and `--bedeq unit` agree to every printed digit.

The constraint operator itself is **second order** — measured directly against
an analytic flow exactly tangent to the bed (`bedorder.py`): observed orders
2.09, 2.03, 2.15.

### Beam angle vs Garcia et al. Fig 9

Best available numbers, at 512×201 (1024×201 at ω/N = 0.2), Lx = 6 km (12 km at
0.2), sponge verified independent, fit window reaching ~0.75 of a bounce:

| ω/N | ours | GCCOM Fig 9 | nonhydro−hydro gap |
|---|---|---|---|
| 0.2 | +0.43 | +1.81 | 0.23 |
| 0.4 | +0.17 | +2.00 | 1.78 |
| 0.6 | +0.09 | +1.07 | 5.91 |
| 0.8 | −0.69 | +1.14 | 14.47 |

GCCOM's points are **digitized** from Fig 9 (no table exists), good to ~±0.1°;
treat differences under 0.3° as ties. The ω/N = 0.2 point cannot discriminate
hydrostatic from nonhydrostatic at all — the gap there is smaller than either
model's bias.

**These numbers are not converged.** See §3.

## 3. THE OPEN PROBLEM — read this before doing anything else

### 3.1 The beam angle converges at first order to a nonzero limit

Four-grid series, Lx = 6 km, `lsl 0.20`, window `[133,1000]`, ω/N = 0.6:

| nx×nz | 256×101 | 384×151 | 512×201 | 640×251 |
|---|---|---|---|---|
| bias | +1.86 | +0.70 | +0.09 | −0.22 |

Fitting `e = A·h^p + e0` in `h = dz`: `p = 1` gives `e0 = −1.63` and reproduces
all four points to **0.023°**; `p = 2` is five times worse. Same shape at
ω/N = 0.8 (`+1.64, +0.01, −0.69, −1.13`). So the +0.09 in the table above is the
curve **crossing zero on its way down**, not a converged answer.

### 3.2 The seiche test localises it to the core discretization

`--seiche I J` runs a free standing internal wave in a flat rectangular box —
no forcing, no sponge, no topography, no fit window, no tracker — against the
exact `ω = N k/√(k²+p²)`, `k = Iπ/Lx`, `p = Jπ/D0`. `--bulge 0` makes the
*domain* rectangular while leaving the *grid* curvilinear (xi-lines still tilt
148 m, because `betaOf` stretches top and bottom differently).

Mode (4,2), `spp = 400`, relative frequency error:

| grid | curvilinear (wander 148 m) | near-Cartesian (bt=bb=0.1, wander 0.5 m) |
|---|---|---|
| 128×51 | +2.91e-3 | +5.04e-3 |
| 192×76 | +3.95e-3 | +6.02e-3 |
| 256×101 | +4.32e-3 | +6.36e-3 |

**The error grows under refinement, on both grids.** The scheme converges to
something, but not to the analytic dispersion relation. This is almost
certainly the root cause of §3.1 and of every unexplained offset in this
project, and it explains why seven targeted hypotheses died — each tested a
component that was not the cause.

Timestep is **not** involved: at fixed 128×51, `spp` 100/200/400/800 gives
3.066e-3, 2.941e-3, 2.910e-3, 2.902e-3 — second-order convergence to ~2.90e-3.
Use `spp = 400`; temporal error is then ~1e-5.

### 3.3 What to try next, in order

All are cheap on the seiche (seconds to a minute per run).

1. **Mode-number scan** — `--seiche 2 1`, `4 2`, `8 4`, `4 8` at fixed grid.
   If the error scales with `k·dz` it is an interpolation/symbol problem; if it
   is roughly constant it is something else.
2. **`--order 4`** on the seiche. The beam-angle version of this test showed
   nothing, but the beam was confounded; this one is not.
3. **`--alpha`** sensitivity. With no sponge, the lateral walls are held *only*
   by the Robin pressure rows, and `alpha = 1e-4` is the sole stand-in for a
   solid wall. If the frequency error moves with `alpha`, the lateral boundary
   treatment is implicated.
4. **The buoyancy interpolation.** `Idf_w · Ic_z` has symbol exactly
   `cos²(k·dz/2)` (measured), so `N_eff = N·cos(k·dz/2)`. That is second order,
   *negative*, and vanishing — it explains the coarse-grid overshoot in §3.1 but
   **cannot** produce a growing positive error, so it is a contributor, not the
   answer.
5. If the seiche error is traced and fixed, **re-run §3.1's four-grid series**.
   Every beam number in §2 must be requoted afterwards.

Do not publish Fig 9 until this resolves, or publish it with the resolution
stated and error bars spanning the extrapolation.

## 4. Dead ends — do not re-test

Each eliminated with runs on record for **this** solver (several were
previously "known" only from the 3D `cone3d` code, which is not the same
discretization).

| hypothesis | evidence against |
|---|---|
| Measurement / tracker bias | `synth_tracker.py`: analytic beam of known angle on the real grid through the identical fit block → ±0.08° at all four frequencies, flat from 20 cells down to 0.6 cells across the beam, and under 40% cross-beam asymmetry |
| Sub-cell peak quantization | parabolic peak refinement shifts the angle by 0.02–0.06° |
| Timestep | angle identical at `spp` 60/120/240; seiche second order in time |
| Operator order (beam) | `k` = 2/4/6 within 0.25° |
| Curvilinear metrics | bias is *anti*-correlated with grid distortion — worst on the near-Cartesian grid; seiche error also worse there |
| `r_x` / HPGF | `r_x` varied 10× (0.709→0.070) at fixed `dz` changes the answer <0.3°; running at `r_x = 1.42` via `--forcegrid` makes every measure *better*. `r_x` is Shchepetkin & McWilliams' criterion for the **hydrostatic** PGF in sigma models; this solver is non-hydrostatic and takes pressure from the Poisson projection, so there is no hydrostatic PGF to reformulate. **Do not build HPGF for this solver.** |
| Near-field curvature | inner window edge from 0.10→0.40 of a bounce, outer fixed: 0.07° at ω/N = 0.4, 0.28° at 0.6 |
| Grid-locked oscillation on the maxima | residual about the fitted line is one arch, not an oscillation; amplitude scales as `dx^1.4` |
| Bed constraint operator | second order, measured (`bedorder.py`) |
| Nonlinearity / finite amplitude | solver has no advection; 4× `u0` change is bit-identical |

## 5. Things that control the measurement (learned the hard way)

1. **Domain size.** At Lx = 3 km the sponge needed to suppress lateral returns
   (`lsl 0.30–0.40`) reaches within 600 m of a `[200,500]` window; ω/N = 0.8
   then moves 1.14° between two sponge settings that both pass the fit test.
   At 6 km with `lsl 0.20` the spread is 0.13–0.26°. **Warning:** at 3 km,
   ω/N = 0.4 had a sponge spread of only 0.18° and was still wrong by 1.83°
   against the 6 km answer at identical `dx` and `dz`. *Local parameter
   insensitivity is not correctness.*
2. **Fit window — the outer edge.** The angle depends strongly on how far out
   the window reaches (2.12° spread at ω/N = 0.6) and weakly on where it starts
   (0.28°). The paper's fixed `[200,500]` m is 0.04–0.10 of a bounce at
   ω/N = 0.2 and 0.27–0.67 at 0.8 — never the same measurement twice. Define
   windows as bounce fractions.
3. **Sponge stability.** The relaxation is explicit, so `dt/taus` is a
   stability number: 1.0 runs, 2.5 → `max|u| = 1.0`, 5.0 → `5.1e+08`. Guarded
   now. This previously looked like a physics blow-up.
4. **Topography sampling.** `gridGen` samples the Gaussian at node positions,
   so the *represented* crest is 8.6 / 16.3 / 19.0 / 19.5 / 19.7 m at
   nx = 64 / 128 / 256 / 384 / 512 against a requested `ab = 20`. Below
   nx ≈ 256 refinement changes the topography, not just the truncation error.

## 6. Tooling

| file | purpose |
|---|---|
| `iwbcurv.py` | the solver. Flags added this session: `--bedmode --bedeq --lsl --taus --order --ab --lb --u0 --bulge --seiche --forcegrid` |
| `paper_repro.ps1` / `.sh` | sections 0–13; logs to `paper-logs/s{sec}_{nnn}.log` |
| `compare_paper.py` | tabulates logs vs theory and vs digitized Fig 9; convergence series grouped by (Lx, lsl, window); power-law extrapolation with uncertainty; sponge-sensitivity report |
| `view.py` | fields on the real grid; geometry read from the npz; `--xlim`, `--total` saturation warning |
| `fig9.py` | reproduces Fig 9 from the logs, with discriminating-power reporting |
| `gridplot.py` | the curvilinear grid, ridge close-up, wall bulge |
| `bedorder.py` | truncation order of the bed constraint against an exactly-tangent flow |
| `wiggle.py` | residual of the tracked maxima about the fitted line |
| `synth_tracker.py` | tracker vs a known-angle beam |
| `test_bedrows.py` | splice algebra, no Octave needed |

`grids/iwbridge/left.m` and `right.m` now read a global `BULGE` (default 0.15)
so the driver can set it. Hand-editing those files once left them at
`bulge = 0`, which silently turned several runs into a different experiment
whose numbers were reported before the corruption was noticed.

Harness lessons worth keeping: never pipe a run's stderr into a filter (a
failing run then prints nothing and the sweep looks instant); prefix log names
with the section number so re-running one section cannot clobber the rest;
refuse to launch with a blank argument.

## 7. Monterey — planned, not started

Full analysis in `MONTEREY-SETUP.md`. Headline: the regime is not a rescaled
toy problem.

* **Rotation is not optional.** At 36.8°N, `f/ω = 0.62` for M2. Beam slope is
  `√((ω²−f²)/(N²−ω²))`, **22% below** the non-rotating value at every N. The
  solver has no `f`. A 2D transect with rotation needs the third velocity
  component (`du/dt += f·v`, `dv/dt = −f·u`); it does not touch the projection.
* **Supercritical.** Canyon flanks give γ = 3–27 against the toy ridge's 0.42.
  Energy reflects rather than radiating; **the beam-angle diagnostic does not
  transfer** and Monterey needs its own.
* **Scale.** Bounce is 27–136 km, so the domain is 100–150 km and the problem
  is 200k–800k unknowns. PARDISO matters here in a way it did not for the toy.
  `pyamg` is recorded as a dead end.
* **`bnd2` divides by `x_xi`.** Fine on a 0.32 slope; degenerate on a
  near-vertical canyon wall. Rewrite `C` as `n_x·u + n_z·w = 0` with the
  cofactor normal — same condition, bounded at any slope, and only changes how
  three coefficients are computed. `bedorder.py` is the harness to verify it.
* **N(z).** `b −= dt·N²·wc` uses a scalar; Monterey N varies ~5× with depth,
  which is what refracts the beam onto the canyon wall.
* Suggested stage 1 is a **synthetic canyon** with realistic slope and depth:
  it exercises all of the above and keeps an analytic expectation in reach,
  whereas real bathymetry does not.

Data received but not yet examined: `mry_offshore_rgb_root.zip` and
`MRY_profile.txt`. The latter is a 7-Gaussian fit over roughly x ∈ [−110, 0]
with `a1 ≈ 12.7`, which looks like **temperature (°C) vs depth (m)**, not
bathymetry. Several coefficients (`a2`, `a3`, `a5`) have confidence bounds
straddling zero, so the individual Gaussians are not identified even if the sum
fits. If it is a temperature profile it gives N(z) directly. Confirm what it is
before using it.

## 8. Process notes

Seven hypotheses died this session — `r_x`, near-field window position,
measurement convention, a grid mode, the constraint's half-cell offset,
quantization, nonlinearity. Every one looked convincing on partial data, and
the failure mode was identical each time: **a correlation measured across a
knob that also moved something else.** The two that survived scrutiny (the
sponge at 3 km, and now the seiche) were found by running a control that held
everything else fixed.

Concretely: before believing any trend here, ask what else changed when the
knob moved, and run the control first. And prefer a diagnostic with an exact
answer (the seiche) over a derived one with four confounds (the beam angle) —
switching to the seiche found in twenty minutes what the beam angle had hidden
for the whole session.
