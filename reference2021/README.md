# reference2021 — the mimetic GCCOM as it stood in October 2021

Read-only reference, not part of the build. These are the files from
`gccom-mole-2021.zip` that the four surviving drivers actually depend on, plus
the one-sided operators and upwind helpers they were built from. The full
archive has 587 files across many generations (v8–v14, `fullymimetic`,
`gridflip`, `petrelis`); this is the subset that matters.

## The drivers, newest first

| file | date | what it is |
|---|---|---|
| `seamount_test_stratified_curvilinear.m` | 2021-10-18 | last file touched: stratified seamount on a 3-D curvilinear grid. Uses `grad3DCurv`, so it ran with the broken `GI13` (fixed in csrc-sdsu/mole#466) |
| `LockRelease3Dv7fullymimeticmodular.m` | 2021-10-14 | the validated case; output named `LE3D_401x6x101`, matching Table 4.1 of the dissertation (Fr 0.705, 0.29% from theory) |
| `Main.m` | 2021-10-13 | the shared time loop: mimetic upwind advection, `InterpVtoV` between Arakawa positions, predictor, Poisson, correction, flux-form upwind temperature, EOS, buoyancy body force |
| `beam2Dv3.m` | 2021-08-13 | 2-D internal-wave beam: Gaussian ridge, sponge, tidal forcing, Coriolis including the non-traditional term. The experiment §4.2.3 reports as failing |

## Supporting pieces

`d2D.m`, `d2Dc.m`, `d3D.m` — the one-sided gradient operators from the thesis,
designed by the author and implemented by J. Corbino. **Not in MOLE today**;
MOLE has only `sidedNodal.m` (1-D, first order).

`MimeticUpwind*.m` — the upwind assembly. `InterpVtoV*.m` — interpolation
between Arakawa C positions. `applyboundaries*.m` — boundary conditions per
experiment. `SSPRK10*.m` — the time integrator. `CalcMetrics*.m`, `Neumann*.m`,
`expand2D.m` — helpers.

## Known defects (see ROADMAP.md for detail and fixes)

1. `SSPRK101D.m` / `SSPRK102D.m` reduce **exactly** to forward Euler: `RHS` is
   passed as a fixed array and never re-evaluated. Verified symbolically.
2. `d2D`/`d3D` are first order (`[-1, 1]/dx`) — the thesis's own diagnosis of
   why the beams dissipate.
3. `beam2Dv3.m` upwind coefficients mix components: `a_pls = max(Ui,0)` with
   `a_neg = min(Uj,0)`, applied to both directions.
4. Buoyancy reaches the w-faces by plain two-point averaging with the bed face
   omitted — not adjoint to the advection, the same class of defect found
   independently in this repository's linear solver.

None of these is fundamental. All four are small.
