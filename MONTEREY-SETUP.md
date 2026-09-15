# Monterey case — parameter regime and what the solver needs

Written before any Monterey code exists, from the M2 numbers at 36.8 N.
Everything here is arithmetic plus a read of `iwbcurv.py`; nothing has been
run. Treat the code-gap list as a proposal, not a finding.

---

## 1. The regime is nothing like the toy problem

| | toy (Garcia §3.3) | Monterey M2 |
|---|---|---|
| forcing | ω/N = 0.2–0.8 | M2, ω = 1.405e-4 |
| rotation | f = 0 | f = 8.74e-5, **f/ω = 0.62** |
| beam slope | 0.75 (at 0.6) | 0.011–0.055 |
| bounce | 1.3 km | **27–136 km** |
| topography | Gaussian, 20 m on 1000 m | canyon, ~1000 m relief |
| criticality γ | 0.42 (subcritical) | **3–27 (supercritical)** |

Beam slope with rotation is `tan φ = sqrt((ω²−f²)/(N²−ω²))`:

| N (1/s) | φ, f=0 | φ with f | change | bounce (D0=1500 m) |
|---|---|---|---|---|
| 0.0020 | 4.03° | 3.16° | −22% | 27 km |
| 0.0035 | 2.30° | 1.80° | −22% | 48 km |
| 0.0050 | 1.61° | 1.26° | −22% | 68 km |
| 0.0070 | 1.15° | 0.90° | −22% | 95 km |

Criticality γ = topographic slope / beam slope:

| topo slope | N=0.002 | N=0.0035 | N=0.005 | N=0.007 |
|---|---|---|---|---|
| 0.05 | 0.9 | 1.6 | 2.3 | 3.2 |
| 0.10 | 1.8 | 3.2 | 4.5 | 6.4 |
| 0.20 | 3.6 | 6.4 | 9.1 | 12.7 |
| 0.30 | 5.4 | 9.5 | 13.6 | 19.1 |

The toy ridge was subcritical, so the beam radiated cleanly away and the
whole validation rested on measuring its angle. Canyon flanks are strongly
supercritical: energy reflects back rather than radiating, and the beam angle
stops being the quantity of interest. **The validated measurement does not
transfer.** Whatever the Monterey target is, it needs its own diagnostic.

## 2. Code gaps, in the order they bite

### (a) No Coriolis — must be added
`f/ω = 0.62`, and rotation cuts the beam slope by 22% at every N. The solver
has no `f`. A 2D transect with rotation needs the third velocity component:
`du/dt −= ... + f·v`, `dv/dt = −f·u`, with `v` co-located with `u`. No
pressure gradient in y for a transect, so this does not touch the projection —
it is an extra field and two terms. Small, but it changes every answer.

### (b) Constant N only — must become N(z)
`b = b − dt·N²·wc` uses a scalar. Monterey N varies ~5× between thermocline
and deep water, and that is exactly what refracts the beam and sets where it
hits the canyon wall. Needs `N²` as an array on cell centres. Small change,
large physical consequence.

### (c) The bed constraint divides by `x_xi`
`bnd2` stores `rat = z_xi/x_xi` and enforces `w = rat·avg(u)`. On the toy
ridge the max slope is 0.32, so `rat` is well behaved. On a near-vertical
canyon wall `x_xi → 0` and `rat` blows up. The constraint should be
rewritten as `n_x·u + n_z·w = 0` using the cofactor normal, which is bounded
for any slope and is the same condition. `C` is already assembled as a
matrix, so this is a change to how its three coefficients are computed, not
to the machinery around it. **This is the one I would fix first**, because it
is cheap and it removes a failure mode rather than adding physics.

### (d) Barotropic forcing over variable depth
`ubc = u0·sin(ωt)` is uniform in x. Over real bathymetry, depth varies by an
order of magnitude and it is the transport that should be uniform, not the
velocity — `u_bc(x) = U(t)/D(x)`. With uniform `u`, the forcing injects
spurious convergence wherever the depth changes.

### (e) Cost
At N = 0.0035 the bounce is 48 km, so a domain that contains one bounce is
100–150 km. Holding cells wide enough that a beam of slope 0.03 crosses
several of them vertically, the workable combinations are roughly
`Lx = 100–150 km, nx = 1024–2048, nz = 201–401`, i.e. 200k–800k unknowns.
SuperLU direct will factor that but it is memory-heavy; PARDISO matters here
in a way it did not for the toy. The handoff records pyamg as a dead end, so
if direct becomes infeasible that is an open problem, not a fallback.

## 3. What I need to start

1. **Transect.** Along-canyon axis or cross-canyon, and the endpoints. The
   solver is 2D; a canyon is not, so which 2D slice is being claimed and what
   it is meant to represent needs deciding before anything is built.
2. **Bathymetry.** Which dataset, and can you export a depth-vs-distance
   profile along that transect? `gridGen` resolves boundary curves by name
   from the working directory, so a real profile means writing a `bottom.m`
   that interpolates a stored array — straightforward once the array exists.
3. **Stratification.** A representative N(z), or the CTD/climatology to build
   one from. Constant N will give the wrong beam paths.
4. **What counts as success.** Reproducing an observed feature (the internal
   bores in the Walter et al. papers the Garcia paper cites), matching a
   published model run, or demonstrating the solver runs stably on real
   bathymetry? These need different domains, diagnostics and validation.
5. **Whether stage 1 is a synthetic canyon.** An idealised Gaussian canyon
   with realistic slope and depth would exercise (a)–(d) and the
   supercritical regime without waiting on data, and would give a case where
   an analytic expectation still exists.

## 4. Suggested order

1. Fix (c), the cofactor-normal constraint — cheap, and re-run the toy
   validation to confirm nothing moved.
2. Add (a) Coriolis, and validate against the rotating dispersion relation on
   the existing Gaussian ridge with f set to Monterey's value. The beam-angle
   test still works there because that ridge is subcritical.
3. Add (b) N(z), validate against a two-layer analytic refraction case.
4. Then (d) and real bathymetry.

Steps 1–3 all have analytic checks available on the toy geometry. Step 4 does
not, which is the argument for doing them in that order.
