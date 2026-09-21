#!/usr/bin/env python3
"""apply_gi2_row.py -- switch GI13's centre->node boundary row to the GI2
convention, [1, .5, -.5], as requested in the Copilot review of PR #466.

  python apply_gi2_row.py            (run from the root of mole-fork)

Measured on grad3DCurv (gi13_boundary_compare.m), the two rows are
indistinguishable in 3-D: interior faces identical, xi-boundary faces within
1% at n=49 (ratio 1.268 -> 1.073 -> 1.022 -> 1.010), both second order. So the
choice is convention, and matching GI2 keeps the 2-D and 3-D operators
consistent at boundaries.

The replacement is checked, not assumed: if the expected lines are not found
exactly once, nothing is written. Python's text mode handles the CRLF working
copy on Windows transparently.
"""
import sys

PATH = 'src/octave/GI13.m'

OLD_ROWS = ("    A(1, 1) = 1.5;    A(1, 2) = -0.5;\n"
            "    A(N+1, N) = 1.5;  A(N+1, N-1) = -0.5;\n")
NEW_ROWS = ("    A(1, 1) = 1;      A(1, 2) = 0.5;      A(1, 3) = -0.5;\n"
            "    A(N+1, N) = 1;    A(N+1, N-1) = 0.5;  A(N+1, N-2) = -0.5;\n")

OLD_DOC = ("% centre -> node: N values on centres -> N+1 values on nodes.  Midpoint average\n"
           "% in the interior, linear extrapolation at the two ends.\n")
NEW_DOC = ("% centre -> node: N values on centres -> N+1 values on nodes.  Midpoint average\n"
           "% in the interior; at the two ends the three-point row [1, .5, -.5], which is\n"
           "% the boundary row GI2 uses once Q's 0.5 weights are factored out.  Both this\n"
           "% and plain linear extrapolation [1.5, -.5] are exact for linear fields and\n"
           "% second order; this one keeps the 2-D and 3-D operators consistent.\n")

try:
    s = open(PATH).read()
except FileNotFoundError:
    sys.exit(f"{PATH} not found -- run from the root of mole-fork")

for old, name in ((OLD_ROWS, 'boundary rows'), (OLD_DOC, 'P docstring')):
    n = s.count(old)
    if n != 1:
        sys.exit(f"expected the {name} exactly once, found {n}; nothing written. "
                 f"Are you on fix/gi13-index-map, and is it already switched?")

s = s.replace(OLD_ROWS, NEW_ROWS).replace(OLD_DOC, NEW_DOC)
open(PATH, 'w').write(s)
print(f"{PATH}: boundary row switched to [1, .5, -.5] (GI2 convention)")
