# HANDOFF v2 — iwbcurv, September 2026

Supersedes the previous HANDOFF for everything it covers. §5 of the old
document (dead ends) is still valid and is extended here; §7.3 (the bed BC
redesign) is DONE; §8.1 (stability) is CLOSED; §8.2 (beam angle) is
essentially closed with one check outstanding.

Environment used throughout: MOLE commit `b0b9ff7` (2026-09-03), Octave
8.4.0, Python 3.11–3.12, numpy 2.4, scipy 1.17, SuperLU. Reproduced
independently on two machines (Linux container and Windows/anaconda) with
matching numbers, including `DIVERGED t/T = 9.00` for the old penalty.

---

## 1. What was fixed

### The bed boundary condition (was §7.3)

The old `--bedmode projection` was a penalty in which three different
operators were in play: `robinBC2D`'s computational-eta normal derivative on
the LHS, the physical contravariant flux over `h_eta` on the RHS, and
`bnd2`'s `w = (z_xi/x_xi)·avg(u)` in `bc_residual()`. `--bsc` was one scalar
trying to reconcile operators that differ point by point on a curved bed.

Replaced by `--bedmode constraint` (now default). `C` is built from the same
`bnd2` arrays the residual uses; `D*G` is empty on the ghost ring (verified,
0 nnz on all 162 boundary rows), so the bottom/top ghost rows carry `C*G`
with rhs `C*u*/dt`. Then `C u = C(u* − dt·G·φ) = 0` exactly, no coefficient.

Results: `div ~1e-16`, `bed ~1e-15`, flat over 40 periods, all four
frequencies, Lx = 6 km. The old penalty at `bsc = 12` diverges at
`t/T = 9.00` on the identical grid.

Row equilibration on the constraint rows is bookkeeping, not a tuned
coefficient: `--bedeq auto` and `--bedeq unit` agree to every printed digit.

### Guards added, all of which fail loudly

* `C` checked against the `bc_residual` loop on a random field at startup.
* `D*G` on the bed ghost rows printed; warns if non-empty.
* Static projection test before any timestep — both residuals, exits if either
  fails. The whole design validates in under a second.
* Sponge stability: the relaxation is EXPLICIT, so `dt/taus` is a stability
  number. Measured: 1.0 runs, 2.5 → `max|u| = 1.0`, 5.0 → `5.1e+08`. Now
  refused above 1.0 with the required `--taus`/`--spp` printed. This
  previously looked like a physics blow-up.
* `--mole` path checked for `gridGen.m` before Octave is touched (Octave's
  `addpath` does not error on a missing directory).

## 2. Current results — Garcia et al. (2019) §3.3

Fully curvilinear grid (the paper used sigma), Lx = 6 km (12 km at
ω/N = 0.2), finest grid, sponge verified independent, fit window reaching
~0.75 of a bounce:

| ω/N | ours | GCCOM Fig 9 (digitized) | theory |
|---|---|---|---|
| 0.2 | +0.43 | +1.81 | 11.54 |
| 0.4 | +0.17 | +2.00 | 23.58 |
| 0.6 | +0.09 | +1.07 | 36.87 |
| 0.8 | −0.69 | +1.14 | 53.13 |

mean |bias| 0.35° vs 1.51°.

GCCOM's points are DIGITIZED from Fig 9 (no table exists), good to ~±0.1°;
recovered x positions were 0.202/0.399/0.601/0.799 against nominal, which is
the accuracy check. Differences under ~0.3° are a tie.

## 3. The three things that actually control the measurement

Learned the hard way; all three had to be right simultaneously.

1. **Domain size.** At Lx = 3 km this grid needs `Lsl = 0.30–0.40 Lx` to
   suppress lateral returns, which puts the sponge within 600 m of a
   `[200,500]` window. At 3 km, ω/N = 0.8 moves 1.14° between two sponge
   settings that both pass the fit test. At 6 km with `lsl = 0.20` the spread
   is 0.13–0.26°. **Note carefully:** at 3 km, ω/N = 0.4 had a sponge spread
   of only 0.18° and was still wrong by 1.83° relative to the 6 km answer at
   identical dx and dz. Local parameter insensitivity is not correctness.
2. **Fit window, outer edge.** The measured angle depends strongly on how far
   OUT the window reaches (2.12° of spread at ω/N = 0.6) and only weakly on
   where it starts (0.28° with the outer edge fixed). The paper's fixed
   `[200,500]` m is 0.04–0.10 of a bounce at ω/N = 0.2 and 0.27–0.67 at 0.8 —
   never the same measurement twice. Use windows defined as bounce fractions.
3. **Resolution.** Ordinary convergence, but slow, and the coarse-grid answer
   can sit on the far side of theory from the converged one.

## 4. Dead ends — do not re-test

Extending the old §5. All of these were eliminated with runs on record for
THIS solver, not inherited from `cone3d`.

* **Measurement / tracker.** `synth_tracker.py` puts an analytic beam of known
  angle on the real curvilinear grid through the identical fit block:
  ±0.08° at all four frequencies, flat across beam widths from 20 cells down
  to 0.6, and under 40% cross-beam asymmetry.
* **Timestep.** Angle identical at spp 60/120/240 while `max|u|` still moves.
* **Operator order.** k = 2/4/6 within 0.25°.
* **Curvilinear metrics.** Bias is ANTI-correlated with grid distortion —
  worst on the near-Cartesian grid.
* **`r_x` / HPGF.** `r_x` varied 10× (0.709 → 0.070) at fixed dz changes the
  answer by <0.3°. Running with `r_x = 1.42` via `--forcegrid` makes every
  measure BETTER. `r_x` is Shchepetkin & McWilliams' criterion for the
  hydrostatic PGF in sigma models; this solver is non-hydrostatic and takes
  its pressure from the Poisson projection, so there is no hydrostatic PGF to
  reformulate. **HPGF should not be built for this solver.**
* **Near-field curvature.** Pushing the inner window edge out from 0.10 to
  0.40 of a bounce, outer edge fixed, changes the answer by 0.07° at ω/N =
  0.4 and 0.28° at 0.6.
* **Grid-locked oscillation on the maxima.** The residual about the fitted
  line is one arch, not an oscillation, and its amplitude scales as dx^1.4.

Four hypotheses died in one session (`r_x`, near-field window position,
measurement convention, grid mode). Each looked convincing on partial data.
The pattern in every case: a correlation measured across a knob that also
moved something else. Run the control first.

## 5. Open

1. **Section 13 — is the long-window agreement a limit or a crossing?**
   The grid series for the exact windows in §2 still drifts down:
   ω/N=0.6 `[133,1000]` gives +1.86, +0.70, +0.09 (extrap −0.59); ω/N=0.8
   `[75,562]` gives +1.64, +0.01, −0.69 (extrap −1.22). 640×251 tests it.
   **Everything in §2 is provisional until this returns.**
2. **Barotropic phase.** At the `u_bc = 0` snapshot phase the interior
   depth-mean is −0.27 u0, near-uniform in x, at both frequencies. Garcia et
   al. plot TOTAL `u/u0` at that phase against ±0.15 and their beams are
   visible, so their interior depth-mean is ~0 there. Ours is a quarter of
   the forcing amplitude out of phase. Does not affect the beam angle
   (rms-based, insensitive to a uniform offset), but our Fig 8 only works in
   baroclinic form. Most likely in how the sponge relaxes the interior toward
   `u_bc`, since the sponge sets the barotropic response in a rigid-lid
   closed domain.
3. **`cone3d` above amp ~0.15** — port `--bedmode constraint`.
4. **Monterey.** No longer blocked on stability or on `r_x`. The `r_x` abort
   should become a warning; keep the folded-Jacobian abort, that one is real.
   Before committing to a grid, run `--gridonly --forcegrid` on a real
   transect.

## 6. Tooling

| file | what |
|---|---|
| `iwbcurv.py` | solver. New: `--bedmode --bedeq --lsl --taus --order --ab --lb --forcegrid` |
| `paper_repro.ps1` / `.sh` | sections 0–13; logs to `paper-logs/s{sec}_{nnn}.log` |
| `compare_paper.py` | tabulates logs vs theory and vs digitized Fig 9; convergence series grouped by (Lx, lsl, window); sponge sensitivity |
| `view.py` | fields on the real grid; geometry read from the npz; `--xlim`, saturation warning |
| `fig9.py` | reproduces Fig 9 from the logs |
| `wiggle.py` | residual of the tracked maxima about the fitted line |
| `test_bedrows.py` | splice algebra, no Octave needed |
| `synth_tracker.py` | tracker vs a known-angle beam |

Harness lessons worth keeping: never pipe a run's stderr into a filter (a
failing run then prints nothing and the sweep looks instant); prefix log names
with the section so re-running one section cannot clobber the rest; refuse to
launch with a blank argument.
