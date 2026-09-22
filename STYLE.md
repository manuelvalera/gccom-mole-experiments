# Code style

## One statement per line

Adopted from the MOLE review of csrc-sdsu/mole#466. Do not put several
statements on one line, in Octave/MATLAB or in Python.

```matlab
% yes
Ax = P(m);
Ay = Q(n);
Az = speye(o);

% no
Ax = P(m); Ay = Q(n); Az = speye(o);
```

This also rules out the one-line block forms:

```matlab
% yes
if val == 0
    continue;
end

% no
if val == 0, continue; end
```

The same applies to Python:

```python
# yes
oc = Oct2Py()
oc.addpath(mole)

# no
oc = Oct2Py(); oc.addpath(mole)
```

Semicolons inside matrix literals (`dc = [1;1;1;1];`), inside strings, and in
comments are not statement separators and are fine.

## Checking

```
python style_check.py src/octave tests/octave      # or any path
python style_check.py . --ext .m                   # Octave only
```

Exit status is 1 if anything is flagged, so it can gate a commit. The checker
ignores comments, string literals (single and double quoted) and bracketed
literals.

Anything written before this rule was adopted is not retrofitted wholesale;
files are brought into line when they are next touched. All four MOLE pull
requests (#466, #468, #469, and the one for #456) were reformatted when the
rule was adopted.
