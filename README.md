# gccom-mole-experiments

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22884121.svg)](https://doi.org/10.5281/zenodo.22884121)

Verification and validation of a two-dimensional, non-hydrostatic, Boussinesq
internal-wave solver on **fully curvilinear grids**, built on the
[MOLE](https://github.com/csrc-sdsu/mole) mimetic operator library. It is a
testbed for reviving the mimetic formulation of the General Curvilinear Coastal
Ocean Model (GCCOM), and its validation target is the internal-wave-beam
benchmark of Garcia et al. (2019), *J. Comput. Sci.* 30:143–156, §3.3.

The solver (`iwbcurv.py`) was built linear and verified that way first; every
section below marked *linear* used `--advect none`, which remains the default
and is bit-identical to those results. Nonlinear advection, rotation, viscosity
and mode-1 boundary forcing were added afterwards, each verified in isolation,
and are used in the last two sections to reproduce the nearshore internal bores
of Walter et al. (2012) on the measured Monterey Bay transect.

## Status (v1.3.0, in progress)

| | status |
|---|---|
| Nonlinear advection | **Verified.** Machine-precision exact cases; second order on stretched grids; lock release Fr = 0.699 against 0.705 published (§5) |
| Rotation | **Verified.** Rotating seiche matches the exact dispersion relation at every latitude, same error as without rotation |
| Nonlinear internal-wave beam | **Runs.** 25 periods, steady, 53.55° against 53.13° — the case the 2021 mimetic GCCOM could not sustain |
| Nearshore bores, Monterey Bay | **Magnitude and vertical structure reproduced** at the 15 m isobath at Walter et al.'s own resolution (dx = 5 m, dz = 1 m); the event asymmetry is **canonical**, not the non-canonical shape they observed (§6) |
| Core discretization | **Verified.** Second order against the exact seiche dispersion relation |
| Dynamical similarity | **Verified.** Two domains at different depths give identical dimensionless results |
| Lateral boundary parameter `alpha` | **Fixed.** The old default dominated the error; now `1e-6` |
| Start-up transient | **Fixed.** A smooth forcing ramp (`--ramp 3`); settled by ~15 periods |
| Energy leak at a sloping bed | **Fixed.** `--buoy energy`; confirmed against a no-stratification control |
| Beam angle, near field | **Converges to theory**: +0.12° at ω/N = 0.8, +0.37° at 0.6 |
| Disagreement with the published Fig 9 | **Explained**: a fit window fixed in metres samples different parts of the beam |
| Speed | 28× faster per run with cuDSS on the GPU |

## 1. The core discretization is verified

A free standing internal wave (`--seiche I J`) in a flat rectangular box — on a
grid that is still curvilinear — has the exact frequency
ω = N k / √(k² + p²). The solver converges to it at **second order**:

<img src="docs/figures/seiche_convergence.png" width="470">

| grid | dz (m) | relative frequency error | local order |
|---|---|---|---|
| 128×51 | 20.00 | −1.8100e-03 | |
| 192×76 | 13.33 | −7.7169e-04 | 2.10 |
| 256×101 | 10.00 | −4.0843e-04 | 2.21 |
| 384×151 | 6.67 | −1.4903e-04 | 2.49 |

**Dynamical similarity.** Two domains that differ in depth by 1.5× but match in
every dimensionless group (`ab/D0`, `Lb/D0`, `Lx/D0`, cells per bounce) give
results identical to two decimals in every distance band, and the same
beam-angle bias:

<img src="docs/figures/similarity.png" width="480">

## 2. Three defects, found and fixed

**The lateral Robin coefficient.** `alpha` only regularises an otherwise
singular pure-Neumann problem at the side walls. At the old default `1e-4` it
dominated the seiche error, which *grew* under refinement. Below about `1e-6`
the error plateaus at its true discretization value.

<img src="docs/figures/alpha_scan.png" width="470">

**The start-up transient.** The tide used to switch on abruptly, kicking the
domain with a broadband transient that rings down slowly — and more slowly on
finer grids, which is the kind of error that can fake a convergence trend.
`--ramp 3` brings the forcing on over three periods; the field is then settled
by ~15 periods.

<img src="docs/figures/spinup.png" width="480">

**An energy leak at a sloping bed.** On fine grids over a steeper ridge, a
spurious mode grew at the crest and radiated beams at a frequency nobody was
forcing (~0.37 N, a shallow ~22° angle), eventually diverging:

<img src="docs/figures/instability_fields.png" width="800">

Ablation located it: the growth was independent of the time step, of `alpha`, of
the operator order and of the grid stretching, but **vanished with `--N 0`** —
so it lived in the buoyancy coupling. The cause is a one-way exchange at the
boundary. MOLE's interpolator pair is not adjoint there: the bed face takes its
buoyancy from the boundary point alone, while the first interior cell takes half
of *its* forcing from the bed face. On a flat bottom `w_bed = 0` and the term
vanishes — which is why the seiche never saw it — but over a slope
`w_bed = slope · u`, and energy enters at the ridge.

`--buoy energy` rebuilds the w-forcing as the exact adjoint of `Idf_w` in the
physical-area inner product (row sums exactly 1, adjoint to ~1e-15). Free runs,
no tide, so any growth is the instability:

<img src="docs/figures/buoyancy_fix.png" width="560">

The fixed solver decays at the same rate as one with no buoyancy coupling at
all, which is as complete as this test can show. A 20 m ridge that used to
diverge at period 79 now runs 150 periods cleanly.

## Validation, animated

The three cases the solver is verified against, run with nonlinear advection
(`--advect full`, second-order upwinding):

| | |
|---|---|
| **Lock release.** The interface collapses into a gravity current and Kelvin-Helmholtz billows roll up along it. Front speed Fr = 0.69 against the inviscid 0.7071 and the 2021 mimetic GCCOM's 0.705. | <img src="docs/animations/lock_release.gif" width="380"> |
| **Seiche.** A standing internal wave in a flat box, the case with an exact frequency. At small amplitude the nonlinear code reproduces the linear result to 0.04%. | <img src="docs/animations/seiche.gif" width="380"> |
| **Internal-wave beam.** Radiation from a ridge under tidal forcing — the experiment the 2021 mimetic GCCOM could not sustain. Steady for 25 periods at 53.55 deg against 53.13 deg theory. | <img src="docs/animations/beam.gif" width="380"> |

Rebuild them with `animate3.ps1` (runs and frames) then `animate_web.ps1`
(sizes for the web). The solver writes frames with `--frames DIR`; `animate.py`
renders any of buoyancy, vertical velocity, horizontal velocity or speed.

## 3. Beam angle

Measured with the perpendicular-slice energy centroid over 0.10–0.40 of a
bounce, with all three fixes in place, six grids per frequency:

<img src="docs/figures/beam_convergence.png" width="520">

| grid | 256×101 | 384×151 | 512×201 | 640×251 | 768×301 | 896×351 | limit |
|---|---|---|---|---|---|---|---|
| ω/N = 0.8 | +2.83 | +1.52 | +0.99 | +0.74 | +0.57 | +0.47 | **+0.12°** (p = 1.63) |
| ω/N = 0.6 | +2.61 | +1.17 | +0.68 | +0.55 | +0.52 | +0.45 | **+0.37°** (p = 2.65) |

**Why the published comparison disagrees.** The measured angle depends on where
along the beam it is fitted, because beyond about half a bounce the beam
approaches the surface and overlaps its own reflection:

<img src="docs/figures/window_dependence.png" width="520">

Garcia et al. fit over a window fixed in metres, `[200,500]`. That is 0.15–0.38
of a bounce at ω/N = 0.6 — the clean near field — but 0.27–0.67 at ω/N = 0.8,
reaching into the contaminated region. In their own window this solver gives
+0.37° at 0.6 and −1.39° at 0.8; over a consistent fraction of a bounce it
gives +0.37° and +0.12°. The disagreement is the convention, not the model.

**Measurement.** `centroid_track.py` tracks the beam by the energy-weighted
centroid on slices perpendicular to the theory direction, which is far more
stable than the column-wise maximum (linear-fit residuals 2.4 m against 3.3 m,
and the argmax series cannot be extrapolated at all). Windows are stated as
fractions of a bounce, D₀/tan θ.

## 4. Speed

<img src="docs/figures/speed.png" width="470">

| | solve | everything else | per step | 20-period run at 896×351 |
|---|---|---|---|---|
| SuperLU, CPU | 124 ms | 23.9 ms | 147.9 ms | 6.1 min |
| cuDSS, CPU loop | 13.1 | 21.4 | 34.6 ms | 1.4 min |
| **cuDSS, GPU loop** | **4.4** | **0.8** | **5.2 ms** | **0.3 min** |

`--device gpu` keeps every field and operator on the card, so nothing crosses
the PCIe bus per step. Saved fields match the CPU run to ~1e-13, and the seiche
reproduces its CPU value exactly. Benchmarked against the alternatives on the
real pressure matrix (`solver_bench.py`): cuDSS is 25× faster per solve than
SuperLU at the same 3e-16 residual, while CuPy's CPU-factor/GPU-solve route is
30× *slower* and its single-precision variant loses too much accuracy to use.
A grid/operator cache removes the ~3 min of Octave grid generation per run.

## 5. Nonlinear advection

Flux form for buoyancy through the mimetic divergence; advective form for
momentum using **contravariant** velocities built from the grid metrics, so the
upwind direction is the one the flow takes through the cell. `--advscheme`
selects `upwind1` (first order, the accuracy of the 2021 one-sided operators),
`upwind2`, `minmod` or `vanleer`. `advect_test.py` checks the operators alone:
exact to 1e-15 on uniform flow, `u = x` and solid-body rotation, and second
order (1.83 → 1.97) on a smoothly stretched grid.

The lock release (`lock_test.ps1`) reproduces the configuration of the 2021
mimetic GCCOM — 0.8 × 0.0656 m, g' = 0.0224 m/s² — and resolves the
Kelvin–Helmholtz billows along the interface:

| scheme, 401×101 | front Fr | vs theory 0.7071 |
|---|---|---|
| upwind1 | 0.6741 | −4.7% |
| minmod | 0.6742 | −4.7% |
| upwind2 | **0.6991** | −1.1% |
| mimetic GCCOM, 2021 | 0.705 | −0.3% |

First-order advection costs 3.7% here; minmod lands on it because at a genuine
discontinuity the limiter reverts to first order. With `--advect full` the beam
of §3 runs steadily for 25 periods at 53.55° — the experiment §4.2.3 of the
2021 dissertation reports as failing. First-order advection alone does not
reproduce that failure here, which points at the other structural difference:
this solver advects only the perturbation and adds the background
stratification analytically, while the 2021 code advected the full temperature
field.

## 6. Nearshore internal bores in southern Monterey Bay

Target: Figure 10 of Walter, Woodson, Arthur, Fringer & Monismith (2012),
*J. Geophys. Res.* 117, C07017 — virtual thermistors at 2, 4 and 6 m above bed
at the 15 m isobath, where a shoaling internal tide arrives as cold bottom
surges. They separate a *canonical* shape (abrupt cold front, gradual warming)
from the *non-canonical* one they observed (sharp drop, continued slow cooling,
then an abrupt warm front), and attribute the difference to the internal
Iribarren number.

**Setup.** The measured transect extended offshore to 20 km (`mry_setup.py
--extend 20 --trim 10`), measured N(z), M2 with rotation at 36.8°N, a mode-1
wave imposed at the offshore boundary (`--mode1 3.7`; c₁ = 0.236 m/s,
λ = 10.5 km, ξ ≈ 2), ν = 1e-4 m²/s as in their SUNTANS runs, and a virtual
mooring sampled every step (`--moor 15 2 4 6`). The grid is **sigma**
(`--bt 0.01 --bb 0.01`): see "The grid" below for why that matters.

**The wave converts from a linear tide to bottom-trapped surges as it shoals.**
The ratio of the 2 mab to 6 mab temperature range crosses 1 between the 80 and
75 m isobaths and climbs to 1.68 at 60 m. The thermistor record at 60 m is
identical at 768×151 and 1536×301.

| isobath | 80 m | 75 m | 70 m | 65 m | 60 m |
|---|---|---|---|---|---|
| range 2 mab / 6 mab | 0.87 | 1.04 | 1.11 | 1.28 | 1.68 |

**The shape at 15 m depends on how mixing is represented.** With a constant
scalar diffusivity κ = 1e-4 the events are canonical — sharp onset, smooth
recovery. Every attempt to lower κ failed at the same moment: the surge
overturns the column at the run-up (N² reaches −1e-3 s⁻² over hundreds of
cells) and the overturn grows at √|N²| with nothing to stop it. Walter et al.
measured O(1 m) overturns at exactly this point. Convective adjustment
(`--kconv 0.01`: vertical diffusion only where N² < 0) lets the background go
to **zero, as in their setup**, and the record changes character: sharp drop,
continued cooling with small internal waves riding on it, then an abrupt warm
front — their non-canonical description almost word for word. The constant κ
had been smearing the warm front.

On the sigma grid with κ = 0 and implicit convective adjustment
(`--kconv 0.1`), bottom-trapped cold events arrive once per tidal period.
**Magnitude** is set by the imposed amplitude; at a = 5 m (ξ ≈ 1.8, chosen, as
they did, to match the observed temperature change) it converges with
resolution and lands near their value:

| 15 m isobath, 2 mab | a = 3.7 m | a = 5.0 m |
|---|---|---|
| 768×151 | 0.202 °C | 0.409 °C |
| 3981×89 (their dx, dz) | 0.215 °C | **0.422 °C** |
| Walter et al. (observed bores) | | ~0.5 °C |

**The event shape does not match.** `events.py` measures, for each event, the
time from onset to minimum and back, and the fastest cooling and warming rates.
Validated on synthetic records of known shape, it reports every configuration
here as canonical or near-symmetric — cooling at least as fast as warming:

| run | cooling | warming | fastest warm / fastest cool |
|---|---|---|---|
| 3981×89, a = 5, κ = 0 + convection | 3.4 h | 3.0 h | 0.25 |
| 768×151, a = 5, κ = 0 + convection | 3.5 h | 2.7 h | 0.33 |
| 3981×89, a = 3.7, κ = 0 + convection | 2.4 h | 3.9 h | 0.09 |
| 768×151, a = 3.7, constant κ, skewed grid | 3.5 h | 3.8 h | 0.74 |

Walter et al.'s non-canonical events warm by ≥ 1 °C in about five minutes,
~0.2 °C/min; the fastest warming here is ~0.006 °C/min. Earlier readings of the
plots as showing an abrupt warm front were wrong, and the metric is what caught
it. Events also last 3–5 h against their 6–20 h.

The leading candidate is geometric: their domain runs to the shoreline and the
non-canonical warm front is the drainback after the surge runs up the slope
past the mooring and stops. Ours ends at the 10 m contour, 177 m inshore of the
15 m mooring, so there is almost no run-up beyond it to drain back from. The
event duration is plausibly structural too: a monochromatic M2 wave gives one
event per 12.4 h, while their events show no fixed tidal phasing and are
attributed to upwelling and bay-scale seiching on 5–10 day scales.

**The grid.** The TFI generator's position-dependent stretching, tuned for the
ridge benchmark, puts different node distributions on the top and bottom
boundaries so interior lines lean. On a 1000 m-deep ridge domain that is ~150 m
of drift and harmless. On an 88 m-deep, 20 km shelf the same machinery drifts
~490 m, and the cells were a median **64° from orthogonal offshore and 81–83° on
the slope** — the solver's own `grid angle` line reported 5.8°–170.5° throughout.
Several operators neglect cross-derivatives on the assumption of near-orthogonal
cells, and every run-up failure sat where the skew was worst. With uniform
stretching (`--bt 0.01 --bb 0.01`) the lines are vertical and the skew drops to a
median 4–7° offshore and 17° on the slope; what remains beyond 9 km is the
measured seabed itself steepening toward the 10 m contour, which terrain-following
coordinates cannot avoid. On the skewed grid every run at 3981×89 died at the
run-up; on the sigma grid the same configuration runs all eight periods. Results
computed on the skewed mry20 grid earlier in this work should be read with that
in mind.

Convective adjustment is implicit — one tridiagonal solve per column per step —
because strengths that keep up with an overturn (0.1–1 m²/s, the range ocean
models use) would need hundreds to thousands of explicit sub-steps in the
thinnest cells.

Things learned the hard way, all now diagnosed by the solver itself: a grid
generator validated on one geometry can be quietly wrong on another;
`alpha` scales as 1/dx² and must change with the domain size; the shoreward
sponge must not cover the slope or the mooring; frames 47 minutes apart cannot
show a five-minute front, hence `--moor`; and explicit diffusion is limited by
the thinnest cells, not the typical ones, hence the metric-based stability
number and `--nusub`.

## MOLE fixes arising from this work

| issue | pull request | subject |
|---|---|---|
| csrc-sdsu/mole#453 | #466, merged | `GI13` index map — `grad3DCurv` did not converge on curvilinear grids |
| csrc-sdsu/mole#454 | #468, merged | warn when the grid handed to the curvilinear operators is left-handed |
| csrc-sdsu/mole#455 | #469, merged | `ttm` initial guess — elliptic grids came out folded |
| csrc-sdsu/mole#456 | #470 | orientation guard for 3-D, vanishing Jacobians, and the Legacy jacobians |

The 2-D results here use **stock MOLE**: none of those code paths is exercised
by `iwbcurv.py`. The buoyancy-coupling defect in §2 is not in a MOLE operator —
`interpol2D` and `interpolD2D` are each correct — but in how the pair is used
together at a sloping boundary.

## Reproducing

MOLE at upstream `main` with #466, #468 and #469 (the `mole` submodule), GNU Octave 8.4, Python
3.11–3.12, numpy, scipy, oct2py; optionally pypardiso, and cupy + nvmath-python
for the GPU path.

```powershell
git clone --recurse-submodules https://github.com/manuelvalera/gccom-mole-experiments.git
cd gccom-mole-experiments
$env:MOLE_SRC = "$PWD\mole\src\octave"        # src\octave, not src\matlab_octave
.\preflight_bulge.ps1                          # grid files must honour --bulge
.\validate.ps1                                 # sections A-C: ~15 min on a GPU
python centroid_track.py val-npz
```

| flag | purpose |
|---|---|
| `--ramp N` | bring the tide on smoothly over N periods |
| `--buoy energy` | energy-consistent buoyancy coupling (§2) |
| `--solver cudss --device gpu` | GPU factorization and an all-GPU time loop |
| `--solver pardiso` | multithreaded CPU solver (`MKL_NUM_THREADS=4` was fastest here) |
| `--u0 0 --noise A --probe` | free run and growth-rate report — any growth is an instability |
| `--cache DIR`, `--no-cache` | grid/operator cache |
| `--advect none\|scalar\|full`, `--advscheme` | nonlinear advection (§5) |
| `--lat`, `--fcor` | rotation |
| `--nu`, `--kappa`, `--nusub` | viscosity and diffusivity, sub-cycled against the thinnest cells |
| `--kconv` | convective adjustment where N² < 0, implicit per column (§6) |
| `--bt`, `--bb` | grid stretching; `0.01 0.01` gives a sigma grid, which the shelf needs |
| `--mode1 A`, `--lock G` | mode-1 boundary forcing; lock release |
| `--moor ISOBATH MAB...` | virtual mooring sampled every step |
| `--savestate T`, `--restart FILE` | checkpoint and resume (resumes by time, so `--spp` may change) |
| `--trace N`, `--checkevery N` | step-level diagnostics, including CFL and minimum N² |
| `--frames DIR` | snapshots for `animate.py` |

| study | scripts |
|---|---|
| everything, in order | `validate.ps1` |
| seiche verification | `seiche_study.ps1`, `seiche_alpha.ps1`, `seiche_report.py` |
| run length and the ramp | `time_test.ps1`, `ramp_test.ps1`, `time_compare.py` |
| beam measurement | `centroid_track.py`, `refit.py`, `curvature.py`, `crest_check.py` |
| the instability | `reproducer.ps1`, `ablate.ps1`, `instab_test.ps1` |
| similarity | `reflect_test.ps1` |
| performance | `solver_bench.py`, `speed_check.ps1`, `gpu_check.ps1`, `gpu_bench.ps1` |
| nonlinear operators | `advect_test.py`, `lock_test.ps1`, `beam_nonlinear.ps1` |
| mode-1 forcing | `mode1.py`, `phase_speed.py` |
| Monterey bores | `mry_setup.py`, `mry_bore.ps1`, `mry_fig10.ps1`, `thermistors.py`, `snapzoom.py` |
| animations | `animate.py`, `animate3.ps1`, `animate_web.ps1` |

Saved fields (`*.npz`) and the cache are not tracked; the scripts regenerate
them. Run logs are included. `STYLE.md` and `style_check.py` describe the
one-statement-per-line convention required for contributions to MOLE itself;
it is not enforced on the scripts in this repository.

## Citation

Cite the concept DOI, which always resolves to the latest version:
[10.5281/zenodo.22884121](https://doi.org/10.5281/zenodo.22884121). Individual
versions: v1.0.1 [10.5281/zenodo.22884885](https://doi.org/10.5281/zenodo.22884885).
See also `CITATION.cff`.

Please also cite Garcia et al. (2019) for the benchmark and the MOLE JOSS paper
(Corbino, Dumett & Castillo, 2024) for the operator library.

## License

GPL-3.0-or-later, matching MOLE.
