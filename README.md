# MOLE curvilinear operator tests

Four independent checks against MOLE, all runnable in about a minute on a laptop.
No MATLAB required — everything runs in Octave, driven from Python via `oct2py`.

## 1. Environment

```bash
conda env create -f environment.yml
conda activate mole
```

If you'd rather not use the yml:

```bash
conda create -n mole -c conda-forge python=3.11 octave oct2py numpy scipy matplotlib pyamg rasterio
conda activate mole
```

On macOS with Homebrew, `brew install octave` then `pip install oct2py` also works.
Verify Octave is visible to Python:

```bash
python -c "from oct2py import Oct2Py; print(Oct2Py().eval('version'))"
```

## 2. Get MOLE

```bash
git clone https://github.com/csrc-sdsu/mole.git
# operators live in  mole/src/matlab_octave
```

## 3. Run everything

```bash
python py/run_all.py --mole /full/path/to/mole/src/matlab_octave
```

That's the whole thing. It prints one block covering all four checks — paste it
straight back.

Optional flags (both default to this repo's directories, so you normally don't
need them):

```
--patch  ./patch      directory holding the patched GI13.m
--grids  ./grids      directory holding grids/seamount/
```

## What each check does

**1 — GI13 index map.** Feeds `GI13` a field whose values encode their own
indices and decodes where each output entry came from. No floating point
involved. Stock MOLE reads 15 of 45 entries from the wrong zeta-plane; the
patch should give 0.

**2 — grad3DCurv convergence.** Same curvilinear grid, `f = x²+y²+z²`, before
and after the patch. Stock gives order ~0.03 (no convergence). Patched should
give ~2.0. Also checks that a Cartesian grid stays exact, so the patch doesn't
break anything that worked.

**3 — gridGen vs jacobian2D.** `gridGen` returns arrays with xi-nodes as rows;
`grad2DCurv` and `jacobian2D` do `[n,m]=size(X)` and want rows = eta. Feed the
output straight through and the Jacobian is negative everywhere — a
left-handed grid, with no warning from anything. Transposing fixes it. Same
|J| both ways, only the sign differs.

**4 — TTM folding.** Thompson–Thames–Mastin elliptic generation on the thesis
Eq. 2.24 seamount folds the grid within 5 iterations, producing genuinely
negative Jacobians in patches while the rest stays positive. Nothing downstream
notices.

## Applying the patch

`patch/GI13.m` is a drop-in replacement with the same signature. Either drop it
into `mole/src/matlab_octave/`, or put its directory earlier on the Octave path
(`addpath` prepends, so it must be added *after* the MOLE directory to take
precedence — that ordering trips people up).

## Files

```
environment.yml
py/run_all.py                          runs all four checks, one report
patch/GI13.m                           the fix
tests/test_GI13_indexmap.m             check 1, standalone Octave
tests/test_grad3DCurv_convergence.m    check 2, standalone Octave
grids/seamount/                        thesis Eq. 2.24 bathymetry, beta modulated
  bottom.m top.m left.m right.m
  stretch.m                            thesis Eq. 2.25-2.26 clustering
  betaOf.m                             position-dependent beta (thesis 5.1)
```

The two `tests/*.m` files run directly in Octave if you'd rather skip Python —
edit the `MOLE_SRC` path at the top of each, then `octave tests/<name>.m`.

## Known-good output

Checks 1, 3 and 4 are deterministic. Check 2's numbers depend slightly on
Octave version but the *orders* should match: ~0.03 stock, ~2.0 patched.

```
1.  15 of 45 output entries read from the wrong zeta-plane

2.  STOCK GI13     n=49  grad rms 5.4146e-01  order 0.03
    PATCHED GI13   n=49  grad rms 1.4859e-03  order 2.09
    Cartesian regression: rms 6.5e-14 / 0 / 4.9e-13

3.  as returned : jacobian2D -> ALL -   |J| 387.3 .. 2747
    transposed  : jacobian2D -> ALL +   |J| 387.3 .. 2747

4.  5 iters -> 9 negative cells;  50 -> 866;  100 -> 1005
```
