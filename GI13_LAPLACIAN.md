# GI13: the assembled 3-D curvilinear Laplacian

Follow-up to the four findings already reported. This strengthens finding #1
by testing the **assembled operator** rather than `grad3DCurv` in isolation.

Tested against `csrc-sdsu/mole` HEAD, Octave. The test is static — it applies
a matrix to a vector. No time stepping, no solver, no boundary conditions,
nothing that could confound it.

## The test

Build the 3-D curvilinear Laplacian the way any solver would:

```matlab
G = grad3DCurvLegacy(k, X, Y, Z);
D = div3DCurvLegacy (k, X, Y, Z);
L = D*G;
```

on a smooth curvilinear grid

```
X = p + amp*sin(2*pi*q)
Y = q + amp*sin(2*pi*p)
Z = r + amp/2*sin(2*pi*p)*sin(2*pi*q)
```

and apply it to `f = x^2 + y^2 + z^2`, whose Laplacian is exactly 6 everywhere.
Report the rms error over the interior (3 layers trimmed from each face).

Note the Legacy pair is used explicitly. On HEAD the four-argument
`grad3DCurv` / `div3DCurv` calls dispatch to `*Legacy`, and `GI13` is called
only from `grad3DCurvLegacy`, so this is the pairing the bug lives in.

## A. Grid refinement (amp = 0.10, k = 2)

| a | stock rms | patched rms | ratio | stock max | patched max |
|---|---|---|---|---|---|
| 12 | 4.5124e+00 | 3.0267e-01 | 14.9x | 1.6044e+01 | 7.3049e-01 |
| 16 | 4.5794e+00 | 1.7038e-01 | 26.9x | 2.2015e+01 | 4.2995e-01 |
| 20 | 4.6659e+00 | 1.0909e-01 | 42.8x | 2.7559e+01 | 3.0298e-01 |
| 24 | 4.7678e+00 | 7.6118e-02 | 62.6x | 3.3023e+01 | 2.1427e-01 |
| 28 | 4.8741e+00 | 5.6356e-02 | 86.5x | 3.7990e+01 | 1.6839e-01 |

Observed order of convergence:

```
  12->16:  stock -0.05   patched +2.00
  16->20:  stock -0.08   patched +2.00
  20->24:  stock -0.12   patched +1.97
  24->28:  stock -0.14   patched +1.95
```

Stock does not converge — the order is slightly **negative**, so the error grows
under refinement, and the max error grows monotonically from 16.0 to 38.0. The
patch gives clean second order, and because stock is flat while the patch
converges, the advantage widens with resolution: 15x at a=12, 87x at a=28.

## B. Grid distortion (a = 20, k = 2)

| amp | stock rms | patched rms | ratio |
|---|---|---|---|
| 0.00 | 2.8508e-13 | 2.8508e-13 | 1.0x |
| 0.02 | 5.1958e-01 | 1.3264e-02 | 39.2x |
| 0.05 | 1.4727e+00 | 3.5984e-02 | 40.9x |
| 0.10 | 4.6659e+00 | 1.0909e-01 | 42.8x |
| 0.15 | 5.9874e+01 | 1.7291e+00 | 34.6x |
| 0.20 | 8.6555e+03 | 1.8389e+02 | 47.1x |

At `amp = 0` the two are **bit-identical at 2.85e-13** — a Cartesian grid has no
off-diagonal metric couplings, and those are exactly the terms routed through
`GI13`. That is why the bug is invisible in any Cartesian test, and it is a
useful regression check: the patch changes nothing where nothing should change.

The moment the grid curves at all — `amp = 0.02` — the stock error jumps to
5.2e-01 while the patch stays at 1.3e-02.

## C. Operator order (a = 20, amp = 0.10)

| k | stock rms | patched rms | ratio |
|---|---|---|---|
| 2 | 4.6659e+00 | 1.0909e-01 | 42.8x |
| 4 | 4.8847e+00 | 3.6526e-02 | 133.7x |
| 6 | 4.8888e+00 | 3.6610e-02 | 133.5x |

Raising the operator order does nothing for stock (4.67 -> 4.89, slightly
worse) because the error is an index-mapping fault, not a truncation error. It
helps the patch by 3x from k=2 to k=4, then saturates — consistent with Abouali
& Castillo (2014), who recommend 4th order as the best accuracy/memory
trade-off and find little further gain at 6th.

## Why this matters more than the isolated operator test

The earlier finding showed `grad3DCurv` alone failing to converge. This shows
the failure survives composition into `L = D*G`, which is the object an
incompressible solver actually inverts every timestep. A pressure Poisson
solve built on stock `GI13` is solving with an operator that is not a discrete
Laplacian on any curved grid, and no amount of grid refinement or operator
order fixes it.

## Files

- `lap_test.py` — the test; `trial(mole, a, amp, patched, korder)` returns
  `(rms, max)` error. The patched `GI13` is embedded, written to a temp dir and
  prepended to the Octave path, so nothing in the MOLE tree is modified.
- `gi13_laplacian_results.png` — the three panels above.
