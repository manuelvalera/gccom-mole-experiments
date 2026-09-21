# gccom-mole-experiments

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22884121.svg)](https://doi.org/10.5281/zenodo.22884121)

Verification and validation of a two-dimensional, non-hydrostatic, Boussinesq
internal-wave solver on **fully curvilinear grids**, built on the
[MOLE](https://github.com/csrc-sdsu/mole) mimetic operator library. It is a
testbed for reviving the mimetic formulation of the General Curvilinear Coastal
Ocean Model (GCCOM), and its validation target is the internal-wave-beam
benchmark of Garcia et al. (2019), *J. Comput. Sci.* 30:143–156, §3.3.

The solver (`iwbcurv.py`) is **linear**: there is no advection term, so every
result scales with the forcing amplitude.

## Status (v1.0.1)

**Verified against an exact solution.** A free standing internal wave
(`--seiche I J`) in a flat rectangular box on a curvilinear grid converges at
**second order to zero** against the exact frequency
ω = N k / √(k² + p²):

| grid | dz (m) | relative frequency error | local order |
|---|---|---|---|
| 128×51 | 20.00 | −1.81e-03 | |
| 192×76 | 13.33 | −7.72e-04 | 2.10 |
| 256×101 | 10.00 | −4.08e-04 | 2.21 |
| 384×151 | 6.67 | −1.49e-04 | 2.49 |

Mode (4,2), `alpha = 1e-6`, time step converged (second order in time). Fitting
e = A·dzᵖ + e₀: p = 2 leaves a residual of 7.8e-08 with e₀ = +5.9e-05; p = 1 is
about 1500× worse.

**The lateral Robin coefficient matters.** `alpha` only regularises an
otherwise singular pure-Neumann Poisson system. At the former default `1e-4` it
dominated the seiche error (+4.3e-03, *growing* under refinement); below ~1e-6
the error plateaus at its true discretisation value. The default is now `1e-6`.
Beam angles move by at most 0.1° between the two values, because the sponge also
holds the lateral boundaries in the forced problem.

**Beam angle — converges, with one open question.** At ω/N = 0.8, measured by
the energy centroid on slices perpendicular to the beam over 0.10–0.40 of a
bounce, the bias against theory is

| grid | 256×101 | 384×151 | 512×201 | 640×251 | 768×301 | 896×351 |
|---|---|---|---|---|---|---|
| bias (°) | +2.76 | +1.45 | +0.82 | +0.46 | +0.16 | +0.02 |

This converges at **first order** (free-exponent fit p = 1.05) to a limit near
**−1.0°**. The first-order behaviour was confirmed by prediction: from the first
four grids, p = 1 predicted +0.19 and +0.00 at the two finest; p = 2 predicted
+0.43 and +0.35; the runs gave +0.16 and +0.02.

The offset is not in the core discretisation (the seiche shows that), not in the
lateral boundary parameter, and not in the represented ridge crest (fixing the
crest at exactly 20.00 m on every grid changes the limit from −1.10° to −1.08°).
Its origin is **open** in this version. Candidates under test: contamination of
the rms-speed diagnostic by the barotropic tide, the terrain-following metrics
near the ridge, and flank sampling.

**Measurement convention.** The measured angle depends on where along the beam
it is fitted; the extrapolated limit spreads by 1.2° (ω/N = 0.6) and 2.2°
(ω/N = 0.8) across fit windows. Windows here are stated as fractions of a
bounce, D₀ / tan θ. Garcia et al. use a fixed [200, 500] m, which is
0.15–0.38 of a bounce at ω/N = 0.6 but 0.27–0.67 at 0.8, so it samples
different parts of the beam at different frequencies.

## MOLE fixes arising from this work

| issue | pull request | subject |
|---|---|---|
| csrc-sdsu/mole#453 | #466 | `GI13` index map — `grad3DCurv` did not converge on curvilinear grids |
| csrc-sdsu/mole#454 | #468 | warn when the grid handed to the curvilinear operators is left-handed |
| csrc-sdsu/mole#455 | #469 | `ttm` initial guess — elliptic grids came out folded |
| csrc-sdsu/mole#456 | — | orientation guard for 3-D and vanishing Jacobians |

The 2-D results in this repository use **stock MOLE**: none of those code paths
(`GI13`, `ttm`, the 3-D Jacobian) is exercised by `iwbcurv.py`.

## Reproducing

Environment: MOLE at commit `1d009d14` (2026-09-14; the `mole` submodule),
GNU Octave 8.4, Python 3.11–3.12, numpy 2.4, scipy 1.17, oct2py.

```powershell
git clone --recurse-submodules https://github.com/manuelvalera/gccom-mole-experiments.git
cd gccom-mole-experiments
$env:MOLE_SRC = "$PWD\mole\src\octave"        # note: src\octave, not src\matlab_octave
.\preflight_bulge.ps1                          # grid files must honour --bulge
```

| step | command | what it establishes |
|---|---|---|
| seiche sweep | `.\seiche_study.ps1` then `python seiche_report.py seiche-logs` | mode, refinement, alpha, order, dx/dz, bed controls |
| alpha | `.\seiche_alpha.ps1` then `python seiche_report.py alpha-logs` | alpha floor, true order, beam sensitivity, long runs |
| beam series | `.\refit_fields.ps1`, `.\finer.ps1` | saved fields at alpha = 1e-6, six grids |
| window / curvature | `python refit.py refit-npz`, `python curvature.py refit-npz` | window dependence, sag |
| tracker | `python centroid_track.py refit-npz` | argmax vs perpendicular energy centroid, extrapolated limits |
| crest | `python crest_check.py refit-npz`, `.\fixedcrest.ps1` | fixed-crest control |
| ridge height | `.\ab_test.ps1`, `python ab_compare.py 10=ab10-npz 20=refit-npz 40=ab40-npz` | tide contamination of the diagnostic |

Saved fields (`*.npz`) are not tracked in git; the scripts regenerate them. Run
logs are included.

## Citation

Cite the concept DOI, which always resolves to the latest version: [10.5281/zenodo.22884121](https://doi.org/10.5281/zenodo.22884121). See also `CITATION.cff`. Please cite Garcia et al.
(2019) for the benchmark and the MOLE JOSS paper (Corbino, Dumett & Castillo,
2024) for the operator library.

## License

GPL-3.0-or-later, matching MOLE.
