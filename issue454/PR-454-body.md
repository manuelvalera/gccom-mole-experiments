## What type of PR is this? (check all applicable)

- [ ] Refactor
- [ ] Feature
- [x] Bug Fix
- [ ] Optimization
- [ ] Example
- [ ] Documentation

## Description

`gridGen` returns `X` with xi-nodes as **rows**. `grad2DCurv`, `div2DCurv` and
`jacobian2D` all do `[n, m] = size(X)` and therefore want rows = eta.
Transposing two axes is an odd permutation, so handedness inverts:

```
gridGen returned size [39 99]     (m=39 horizontal nodes is the ROW count)

as returned  : jacobian2D -> ALL -
transposed   : jacobian2D -> ALL +
magnitudes identical after sorting, max diff 0.000e+00
```

Same grid, same magnitudes to the last bit, opposite sign.

A uniform sign flip is self-consistent for `grad` and `div` alone — `J` flips,
the cofactors flip with it, and gradients still come out exact. It stops being
harmless as soon as anything downstream assumes a fixed physical orientation:
a buoyancy term, a volume weighting, a `sqrt(J)`. And nothing anywhere warned.

### What this PR does, and does not, change

This adds the missing warning. **It does not change `gridGen`'s output**, and
it does not change any numerical result — `J` is byte-identical before and
after.

Changing `gridGen` to return `[n, m]` would be the other half of the fix, but
it is a breaking change for anyone who already compensates by transposing, and
whether that belongs in a patch release or in MOLE 2.0 is a maintainer call.
Happy to open it separately once there's a preference. This PR is deliberately
the uncontroversial half.

### Implementation

New `src/octave/checkGridOrientation.m`, called from `jacobian2D` on both
entry paths (3 arguments dispatches to `jacobian2DLegacy`, 5 arguments computes
inline). `grad2DCurv`, `div2DCurv`, `grad2DCurvLegacy` and `div2DCurvLegacy`
all route through `jacobian2D`, so users who never call it directly still get
told.

Two cases are reported separately:

| condition | meaning | warning id |
|---|---|---|
| `J < 0` everywhere | left-handed grid, domain mirrored | `jacobian2D:leftHandedGrid` |
| `J` changes sign | mapping folds over itself, not invertible there | `jacobian2D:tangledGrid` |

The left-handed message names the likely cause and the fix:

```
The Jacobian is negative everywhere, so the grid is left-handed and the domain
is mirrored.
The curvilinear operators read [n, m] = size(X): rows are eta, columns are xi.
If X and Y came from gridGen, they have xi-nodes as ROWS -- pass X.' and Y.'
instead.
grad and div are self-consistent either way, so this will not show up as an
error; it matters for anything that assumes a physical orientation (buoyancy,
volume weights, sqrt(J)).
Silence with: warning('off', 'jacobian2D:leftHandedGrid')
```

Both are warnings, not errors, and both carry an identifier so they can be
switched off.

The tangled-grid case is slightly beyond what the issue asks for. It was three
extra lines on top of the same `nonzeros(J)` sign count and catches a strictly
worse failure, but say the word and I'll drop it.

## ⚠️ This warning fires on a shipped example

`examples/octave/elliptic2DCurvPeriodic.m` builds a left-handed grid. This is
not a false positive — it is analytic:

```matlab
[T, R] = meshgrid(ts, rs);       % rows = r (eta), cols = theta (xi)
X = R .* cos(T);
Y = R .* sin(T);

J = x_xi y_eta - x_eta y_xi
  = (-r sin T)(sin T) - (cos T)(r cos T)
  = -r                              <- negative everywhere
```

Measured: 0 positive, 441 negative. With `meshgrid(rs, ts)` instead it is
`+r`, 441 positive, 0 negative.

The example still produces the right answer, which is exactly the point the
issue makes — a uniform sign flip is invisible to `grad` and `div`. But it
means this PR makes that example emit a warning.

I have deliberately **not** touched it. Swapping the axes there also swaps
which boundary is periodic (`dc = [0;0;1;1]` becomes `[1;1;0;0]`), which
boundary vectors are which, the `reshape`, and the edge-join before plotting.
That is an invasive rewrite of a working example with no error norm to
validate against, and it is a different change from adding a warning.

Three ways forward, maintainer's choice:
1. Fix the example properly in a follow-up (I'm happy to do it).
2. Add `warning('off', 'jacobian2D:leftHandedGrid')` at the top of that example
   with a comment pointing here.
3. Leave it — examples are not run in CI (only `matlab-tests.yml`, over
   `tests/octave`), so nothing automated breaks.

## Related Issues & Documents

- [x] I have created an Issue that is paired with this PR
- Closes #454

## QA Instructions, Screenshots, Recordings

Verified against `1d009d14` (current `main`) with Octave 8.4.0 on Ubuntu 24.04.

Reproduce the original report:

```matlab
addpath(genpath('src/octave'));
cd src/octave/grids
[X, Y] = gridGen('TFI', 'swan', 39, 99, false);
size(X)                                  % [39 99] -- xi-nodes as rows
d  = nonzeros(jacobian2D(2, X,   Y  ));  % all negative, now warns
dt = nonzeros(jacobian2D(2, X.', Y.'));  % all positive, silent
max(abs(sort(abs(d)) - sort(abs(dt))))   % 0 -- identical magnitudes
```

Confirm no numerical change: `J` from a right-handed grid is unchanged by this
patch (441 nonzeros, same values), and the only edit to `jacobian2D.m` is two
added call lines.

## Keep-open request

- [ ] I am requesting maintainer review for `keep-open`.

Reason:

## Added/updated tests?

- [x] Yes

Adds `tests/octave/testGridOrientation.m`, seven checks:

| test | what it pins |
|---|---|
| `testRightHandedIsSilent` | no warning on a correctly oriented `meshgrid` grid |
| `testTransposedWarns` | `leftHandedGrid` fires when the axes are swapped |
| `testGridGenOutputWarns` | the reported case: `gridGen` output warns, its transpose does not |
| `testMagnitudesAreIdentical` | same `|J|` either way — documents why this is hard to spot |
| `testTangledGridWarns` | `tangledGrid` fires on a folded mesh |
| `testWarningReachesOperatorUsers` | `grad2DCurv` and `div2DCurv` users get it too |
| `testFiveArgPath` | both the 3-arg and 5-arg entry paths are covered |

The silent cases matter as much as the firing ones: a check that warns on
everything is as useless as one that warns on nothing.

## Read Contributing Guide and Code of Conduct

- [x] 📖 I have read the MOLE Contributing Guide
- [x] 📖 I have read the MOLE Code of Conduct

## [optional] Are there any post deployment tasks we need to perform?

The `gridGen` orientation itself, and `examples/octave/elliptic2DCurvPeriodic.m`
— both described above, both awaiting a maintainer preference rather than a
technical answer.
