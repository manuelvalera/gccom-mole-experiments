# GCCOM–MOLE handoff (v4, 7 October 2026)

Upload this file at the start of the next chat. It replaces the earlier handoff
notes; the README and ROADMAP in the repository carry the full story.

## Working conventions

- Commands to run go at the start of every response.
- Code: one statement per line, Python and Octave alike.
- Windows PowerShell, conda env `mole`, repo at
  `C:\Users\dap21\gccom-mole-experiments-2026`
  (github.com/manuelvalera/gccom-mole-experiments, concept DOI 10.5281/zenodo.22884121).
- `$env:MOLE_SRC` points at the MOLE submodule's `src/octave`.
- PowerShell traps met so far: `"$case:"` in a string is read as a drive
  variable (write `${case}`); argparse needs `--n2=-1e-5` / `--dz=-15` for
  negative values; `Tee-Object -Append` writes UTF-16 (capture to a variable,
  then `Out-File -Encoding utf8`); right-click in the console pastes and runs the
  clipboard; a `Select-String` filter can hide the one error line, so check
  unfiltered output when a command returns in seconds.

## Where things stand

**Headline result (README section 6).** On Walter et al. (2012)'s own two
bathymetries and stratification, recovered from their archived SUNTANS setup
(`rectddatafromsuntansmodel.zip`), at their resolution, `iwbcurv.py` reproduces
both bore regimes at the 15 m mooring:

| | steep step, xi ~ 2 | gentle slope, xi ~ 0.2 |
|---|---|---|
| cooling / warming | 4.6 h / 1.2 h | 0.17 h / 6.3 h |
| ratio -> regime | 3.7 -> non-canonical | 0.03 -> canonical |
| 2 / 4 / 6 mab range | 1.06 / 0.70 / 0.37 degC | 1.39 / 1.40 / 1.40 degC |

Remaining gap: the non-canonical warm return takes ~1 h (fastest warming
0.03-0.04 degC/min) where theirs takes minutes (~0.2 degC/min). Leading
suspect: amplitude, a = 10 m run against ~16 m implied by their xi values.

**Run configuration that works** (both Walter cases, 6 periods, ~1.2 h each):

```
--grids grids --gridname walter_xi2 | walter_xi02 --Lx 20000 --D0 81.0
--nprofile mry_N.txt --omega 1.405e-4 --lat 36.8 --forcegrid
--solver cudss --device gpu --buoy energy --advect full --advscheme upwind2
--nu 1e-4 --kappa 0 --kconv 0.1 --nuh 0.05 --kappah 0.05 --u0 0
--mode1 10 --mode1force uw
--nx 3981 --nz 81 --spp 24000 --nper 6 --ramp 1 --taus 1500 --lsl 0.20 --lslr 0.005
--ulim 5000 --bt 0.01 --bb 0.01 --bulge 0
--savestate 2.6 --moor 15 2 4 6 --frames <dir> --framerate 2400 --failframes 60
```

Grids are built with `walter_bathy.py --case noncanonical|canonical` then
`mry_setup.py --bathy <csv> --out grids\<name> --smooth 0`. Checkpoints from
these runs: `mry-fig10\walter_xi2_full_state.npz`, `walter_xi02_full_state.npz`
(t/T = 2.6; restart refuses a different grid via a fingerprint).

**Outputs in the repo:** `docs/figures/walter_xi2_full_15m.png`,
`walter_xi02_full_15m.png`, `walter_bores.gif` (side-by-side animation from
`compare_anim.py`), `grid_walter_bulge015.png`, `grid_walter_bulge0.png`.

## What we learned (do not relearn)

- **Grid.** Two ridge-benchmark defaults wreck shelf grids: position-dependent
  stretching (fix `--bt 0.01 --bb 0.01`) and `--bulge 0.15` (fix `--bulge 0`).
  The "wander" diagnostic misses a symmetric bow; `grid_view.py --skew` on the
  grid a run actually used is the check. Clean sigma grid: worst cell 4.7 deg.
  `--sigma` rebuilds a true sigma grid from any TFI output.
- **Orientation.** The solver forces in the LEFT sponge; `mry_setup.py` now
  always puts the deep end on the left.
- **Bathymetry decides the shape.** On our measured transect the shape never
  went non-canonical (duration ratios 1.0-1.2 with their stratification) through
  changes of resolution, rotation, mixing, stratification and run-up room: it
  climbs 50 -> 15 m over ~1 km where theirs does it in ~400 m.
- **Stratification.** Their model profile is a fit to the MBARI C1 cast of
  1 April 2010; `MRY_profile.txt` / `mry_N.txt` match it within 0.1 degC. The
  C1 18 May 2010 cast (`ctd_profile.py`) is a realism experiment: right bore
  size (0.5 degC), soliton fission, still canonical on our transect.
- **Resolution.** 26 m cells cannot carry the steepening fronts on either
  profile; 5 m can. Halving dt does not help; the limit is spatial.
- **Mixing.** Convective adjustment (implicit, per column) for overturns;
  horizontal mixing 0.05 m2/s for a two-grid-interval mode in x. Both are
  disclosed departures from SUNTANS.
- **Measurement.** Classify events by where the time goes (cooling vs warming
  duration), not by the fastest rates; merge overlapping events
  (`events.py`). Several plot readings were wrong once measured.
- **The mimetic core never failed**: divergence and bed constraint held at
  1e-13 in every run. All failures came from the grid, the wavemaker, or the
  hand-built non-mimetic advection and mixing.

## Next, in order

1. **a = 16 m at full resolution, with an energy-budget diagnostic.** Add to
   `iwbcurv.py`: d/dt(KE + APE) minus forcing work minus explicit dissipation
   (viscous, diffusive, convective, sponges); the residual is numerical
   dissipation. Report it per period and near the warm front. Then run both
   Walter cases at a = 16 (consider `--spp 48000` if CFL climbs).
2. **Measured-transect near-coast reruns** on a `--bulge 0` grid (all earlier
   shelf runs used the bowed grid).
3. **Tag v1.3.0** on Zenodo once 1-2 are in the README.
4. **MOLE submodule** to upstream main (local `GI13.m` edit and untracked test
   are an old draft of #466). Do this AFTER the full-resolution runs: changing
   MOLE source invalidates every cached grid.
5. **Write to the paper's authors** with the result and the warm-front gap.
6. **Tensor Laplacian** (`D_curv M(K, g) G_curv` with metric cross terms) in
   this repo first, replacing the diagonal-only viscosity; precedents:
   Shashkov & Steinberg (1996), Hyman, Shashkov & Steinberg (1997), Boada et
   al. (arXiv 2407.01443). Then energy-consistent (skew-symmetric) advection
   with separate SPD upwind dissipation, as part of Stage 4 (3-D). Decision:
   keep this work in-house for now; no outside contact planned yet.
7. Older items still open: PR #470 (3-D orientation guard + Legacy jacobians)
   awaiting merge; CITATION.cff / .zenodo.json ORCID 0000-0001-8766-2450;
   Zenodo v1.0.0 removal ticket.

## Tools added in this stretch

`walter_bathy.py`, `ctd_profile.py`, `shift_profile.py`, `events.py`,
`thermistors.py` (per-step moorings, cuts diverging tails), `instab_anim.py`
(N2 < 0 and CFL markers), `grid_view.py`, `compare_anim.py`, `snapzoom.py`.
Solver flags: `--moor`, `--savestate/--restart` (by time, grid fingerprint),
`--failframes`, `--trace` (CFL, min N2), `--kconv` (implicit), `--nuh/--kappah`,
`--mode1force uw`, `--sigma`, `--nusub`.

`mry_setup.py` copies `betaOf.m` and `stretch.m` from `grids/iwbridge`, so grid
directories only build inside the repository.
