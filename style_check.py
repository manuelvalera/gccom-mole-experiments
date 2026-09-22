#!/usr/bin/env python3
"""style_check.py -- one statement per line.

  python style_check.py src/octave tests/octave        # check paths
  python style_check.py --fix-hint FILE                # show what to split

House rule, from the MOLE review of PR #466: do not put several statements on
one line. Write

    Ax = P(m);
    Ay = Q(n);

not

    Ax = P(m); Ay = Q(n);

It applies to Octave/MATLAB (`;` and the one-line `if x, y; end` form) and to
Python (`;`). Comments, string literals and matrix literals such as
`dc = [1;1;1;1]` are not statement separators and are ignored.

Exit status is 1 if anything is flagged, so it can gate a commit.
"""
import argparse
import os
import re
import sys


def strip_noise_matlab(line):
    """Blank out comments, strings and [...]/{...} literals."""
    out, i, n = [], 0, len(line)
    in_str = None
    depth = 0
    while i < n:
        ch = line[i]
        if in_str:
            out.append(' ')
            if ch == in_str:
                in_str = None
            i += 1
            continue
        if ch == '%':
            break                                   # comment to end of line
        if ch == '"':
            in_str = '"'                            # Octave double-quoted string
            out.append(' ')
            i += 1
            continue
        if ch == "'" and (not out or not re.match(r"[\w\)\]\.']", out[-1])):
            in_str = "'"                            # transpose vs string start
            out.append(' ')
            i += 1
            continue
        if ch in '[{':
            depth += 1
        elif ch in ']}':
            depth = max(0, depth - 1)
        out.append(' ' if (depth and ch == ';') else ch)
        i += 1
    return ''.join(out)


def strip_noise_python(line):
    out, i, n = [], 0, len(line)
    quote = None
    depth = 0
    while i < n:
        ch = line[i]
        if quote:
            out.append(' ')
            if ch == quote and line[i-1] != '\\':
                quote = None
            i += 1
            continue
        if ch == '#':
            break
        if ch in '"\'':
            quote = ch
            out.append(' ')
            i += 1
            continue
        if ch in '[{(':
            depth += 1
        elif ch in ']})':
            depth = max(0, depth - 1)
        out.append(' ' if (depth and ch == ';') else ch)
        i += 1
    return ''.join(out)


def offenders(path):
    matlab = path.endswith('.m')
    strip = strip_noise_matlab if matlab else strip_noise_python
    bad = []
    for num, line in enumerate(open(path, encoding='utf-8', errors='ignore'), 1):
        code = strip(line).rstrip()
        if not code.strip():
            continue
        # a trailing ';' is just a statement terminator; anything after it is a
        # second statement
        if re.search(r';\s*\S', code):
            bad.append((num, line.rstrip(), 'several statements'))
        elif matlab and re.search(r'\b(if|for|while)\b.*,\s*\S+.*;?\s*end\s*$', code):
            bad.append((num, line.rstrip(), 'one-line block'))
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('paths', nargs='+')
    ap.add_argument('--ext', nargs='+', default=['.m', '.py'])
    g = ap.parse_args()

    files = []
    for p in g.paths:
        if os.path.isdir(p):
            for dp, dn, fn in os.walk(p):
                dn[:] = [d for d in dn if d not in ('.git', '__pycache__', 'mole')]
                files += [os.path.join(dp, f) for f in sorted(fn)
                          if os.path.splitext(f)[1] in g.ext]
        elif os.path.splitext(p)[1] in g.ext:
            files.append(p)

    total = 0
    for f in sorted(files):
        bad = offenders(f)
        if bad:
            total += len(bad)
            print(f"\n{f}")
            for num, line, why in bad:
                print(f"  {num:5d}  {why:18s} {line.strip()[:88]}")
    print(f"\n{len(files)} files checked, {total} line(s) with more than one statement")
    return 1 if total else 0


if __name__ == '__main__':
    sys.exit(main())
