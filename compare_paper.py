#!/usr/bin/env python3
"""compare_paper.py -- tabulate our runs against Garcia et al. (2019) Fig 9.

  python compare_paper.py paper-logs

Scans *.log written by paper_repro, pulls (omega/N, measured angle, fit rms,
max/mean, nz) out of each, and prints them beside theory and beside GCCOM's
own published points.

THE GCCOM NUMBERS BELOW ARE DIGITIZED, NOT TABULATED.  Garcia et al. give the
beam angles only as markers in Fig 9; there is no table.  These were measured
from a 600 dpi rasterization of p.151 by locating the plot box from the axis
rules, eroding the line work away so only the filled markers survived, and
mapping marker centroids through the axes.  Recovered x positions were 0.202,
0.399, 0.601, 0.799 against nominal 0.2/0.4/0.6/0.8, which is the accuracy
check on the mapping -- so the angles are good to roughly +/-0.1 deg, no
better.  Treat a difference below ~0.3 deg from GCCOM as a tie.
"""
import os, re, sys, glob
import numpy as np

# omega/N -> (published beam angle, digitized), from Fig 9
GCCOM = {0.2: 13.35, 0.4: 25.58, 0.6: 37.94, 0.8: 54.27}


def theory(r):
    return np.degrees(np.arcsin(r))          # phi = asin(omega/N), Eq. 30


def scan(path):
    txt = open(path, errors='ignore').read().replace('\r', '')
    g = lambda p: (re.search(p, txt).group(1) if re.search(p, txt) else None)
    if 'BEAM ANGLE' not in txt:
        return None
    return dict(
        ratio=float(g(r'omega/N=([0-9.]+)')),
        nx=int(g(r'omega/N=[0-9.]+\s+(\d+)x\d+')),
        nz=int(g(r'omega/N=[0-9.]+\s+\d+x(\d+)')),
        Lx=float(g(r'Lx=([0-9.]+) km')),
        meas=float(g(r'BEAM ANGLE\s+measured ([0-9.]+)')),
        rms=float(g(r'rms residual ([0-9.]+) m')),
        mm=float(g(r'max/mean ([0-9.]+)')),
        lsl=g(r'Lsl=\d+ m \(([0-9.]+) Lx\)'),
        taus=g(r'taus=([0-9.]+) s'),
        win=g(r'fit over \|x\| in \[([0-9,]+)\]'),
        log=os.path.basename(path))


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'paper-logs'
    rows = [r for r in (scan(f) for f in sorted(glob.glob(os.path.join(d, '*.log')))) if r]
    if not rows:
        raise SystemExit(f"no completed runs found in {d}/")
    # sections 1 and 2 run the same calibrated config, so identical runs appear
    # twice.  Collapse them; keeping both padded the table and would skew any
    # average taken over rows.
    seen, uniq = set(), []
    for r in rows:
        k = (r['ratio'], r['nx'], r['nz'], r['lsl'], r['taus'], r['win'], r['meas'])
        if k in seen:
            continue
        seen.add(k); uniq.append(r)
    ndup = len(rows) - len(uniq)
    rows = uniq

    print(f"\n{len(rows)} distinct completed runs in {d}/"
          + (f"  ({ndup} duplicate(s) collapsed)" if ndup else "") + "\n")
    hdr = ("w/N", "Lx km", "grid", "sponge", "win", "ours", "theory",
           "ours-th", "GCCOM-th", "vs GCCOM", "fit rms", "max/mean")
    print("%-5s %-6s %-9s %-11s %-9s %-7s %-7s %-8s %-9s %-9s %-8s %-8s" % hdr)
    print("-" * 108)
    for r in sorted(rows, key=lambda z: (z['ratio'], z['nz'], z['log'])):
        th = theory(r['ratio'])
        gb = GCCOM.get(round(r['ratio'], 1))
        gbias = (gb - th) if gb else float('nan')
        obias = r['meas'] - th
        verdict = ("tie" if abs(abs(obias) - abs(gbias)) < 0.3
                   else ("BETTER" if abs(obias) < abs(gbias) else "worse"))
        flag = "" if (r['rms'] < 10 and r['mm'] > 5) else "  <- fit not usable"
        print("%-5.1f %-6.1f %-9s %-11s %-9s %-7.2f %-7.2f %+-8.2f %+-9.2f %-9s %-8.1f %-8.1f%s"
              % (r['ratio'], r['Lx'], f"{r['nx']}x{r['nz']}",
                 f"{r['lsl']}/{r['taus']}", r['win'], r['meas'], th,
                 obias, gbias, verdict, r['rms'], r['mm'], flag))

    print("\nAcceptance for a usable fit: rms < 10 m and max/mean > 5.")
    print("Rows failing that are measurements of reflections, not of a beam --")
    print("do not compare their angles to anything.")
    print("\nGCCOM points digitized from Fig 9, good to about +/-0.1 deg;")
    print("differences under ~0.3 deg are a tie, not a win.")

    ok = [r for r in rows if r['rms'] < 10 and r['mm'] > 5]
    if not ok:
        return

    def mean_bias(sel):
        return np.mean([abs(r['meas'] - theory(r['ratio'])) for r in sel])

    gm = np.mean([abs(GCCOM[k] - theory(k)) for k in GCCOM])

    # At the paper's own grid size -- the like-for-like comparison.
    paper = [r for r in ok if r['nx'] == 128 and r['nz'] == 101]
    bypr = {}
    for r in paper:
        bypr.setdefault(round(r['ratio'], 1), []).append(r)
    if len(bypr) == 4:
        pick = [max(v, key=lambda z: z['mm']) for v in bypr.values()]
        print(f"\nAT THE PAPER'S GRID SIZE (128x101), mean |bias| over the four:")
        print(f"   ours {mean_bias(pick):.2f} deg      GCCOM (Fig 9) {gm:.2f} deg")

    # Convergence, per frequency, holding sponge AND window fixed.
    # With three or more grids, report the successive differences and a
    # geometric extrapolation.  Shrinking differences mean the scheme IS
    # converging; converging to a nonzero limit is a different failure from
    # diverging, and conflating the two cost a round here.
    grids = sorted({(r['nx'], r['nz']) for r in ok})
    # taus is set per frequency (T/30) by design, so it is NOT part of the
    # config key -- grouping on it would split one sweep into four tables.
    cfgs = sorted({(r['Lx'], r['lsl'], r['win']) for r in ok})
    if len(grids) > 1:
        notes = []
        for cfg in cfgs:
            sel = [r for r in ok if (r['Lx'], r['lsl'], r['win']) == cfg]
            gs = sorted({(r['nx'], r['nz']) for r in sel})
            if len(gs) < 2:
                continue
            Lx_m = cfg[0]*1000.0
            print(f"\nBIAS vs GRID  (Lx {cfg[0]:.0f} km, lsl {cfg[1]}, win [{cfg[2]}], taus = T/30)")
            print("%-5s " % "w/N" + " ".join("%-10s" % f"{g[0]}x{g[1]}" for g in gs)
                  + " %-18s %-16s" % ("successive d", "extrap e0 (p=1..2)"))
            for k in sorted({round(r['ratio'], 1) for r in sel}):
                v, cells = [], []
                for g in gs:
                    m = [r for r in sel if round(r['ratio'], 1) == k
                         and (r['nx'], r['nz']) == g]
                    if m:
                        b = m[0]['meas'] - theory(k); v.append(b)
                        cells.append(f"{b:+.2f}")
                    else:
                        cells.append("")
                d = [v[i+1]-v[i] for i in range(len(v)-1)]
                ds = ", ".join(f"{x:+.2f}" for x in d) if d else ""
                ext = ""
                # The geometric estimate below assumes a CONSTANT refinement
                # ratio.  The grids actually used (256/384/512/640) refine by
                # 1.50, 1.33, 1.25, so that assumption is wrong here and the
                # estimate is biased.  With three or more points, fit
                # e(h) = e0 + C h^p properly instead and report the spread
                # across fit choices, because a 3-parameter fit to 4 points is
                # near-interpolation and e0 is far less certain than the
                # residuals suggest.
                if len(v) >= 3:
                    hh = np.array([Lx_m/g[0] for g in gs[:len(v)]])
                    ee = np.array(v)
                    cands = []
                    for pf in (1.0, 1.5, 2.0):
                        A = np.vstack([np.ones_like(hh), hh**pf]).T
                        try:
                            a, _C = np.linalg.lstsq(A, ee, rcond=None)[0]
                            cands.append(a)
                        except Exception:
                            pass
                    if cands:
                        ext = (f"{np.mean(cands):+.2f}"
                               f" +/-{0.5*(max(cands)-min(cands)):.2f}")
                        notes.append((cfg, k, float(np.mean(cands)), 0.0))
                if ext:
                    pass
                elif len(d) >= 1 and abs(d[-1]) < 0.10:
                    # already flat: the last refinement moved it < 0.1 deg
                    ext = f"{v[-1]:+.2f}"
                    notes.append((cfg, k, v[-1], 0.0))
                elif len(d) >= 2 and d[-2] != 0:
                    q = d[-1]/d[-2]
                    if 0 < q < 1:            # geometric decay -> finite limit
                        lim = v[-1] + d[-1]*q/(1-q)
                        ext = f"{lim:+.2f}"
                        notes.append((cfg, k, lim, q))
                    elif q >= 1:
                        ext = "diverging"
                        notes.append((cfg, k, None, q))
                print("%-5.1f " % k + " ".join("%-10s" % c for c in cells)
                      + " %-18s %-16s" % (ds, ext))
        fin = [n for n in notes if n[2] is not None and abs(n[2]) > 0.3]
        div = [n for n in notes if n[2] is None]
        if div:
            print("\n   *** NOT CONVERGING ***")
            for cfg, k, _, q in div:
                print(f"   Lx {cfg[0]:.0f} km lsl {cfg[1]} win [{cfg[2]}] w/N={k}: "
                      f"differences are not shrinking (ratio {q:.2f})")
        if fin:
            print("\n   *** CONVERGING TO A NONZERO LIMIT ***")
            for cfg, k, lim, q in fin:
                print(f"   Lx {cfg[0]:.0f} km lsl {cfg[1]} win [{cfg[2]}] w/N={k}: -> {lim:+.2f} deg")
            print("   Differences ARE shrinking, so this is not a grid error --")
            print("   refining further will not remove it.  A limit that varies")
            print("   with the fit window means the beam is not straight at the")
            print("   theory angle over that window, which is a statement about")
            print("   the measurement convention, not about the discretization.")

    # Sponge sensitivity at fixed grid and window.  The answer must not depend
    # on the sponge.  If it does, the sponge is reaching the fit window and the
    # refinement series above is not measuring discretization error.
    print("\nSPONGE SENSITIVITY (same grid, same window, usable fits only)")
    worst = 0.0
    for k in sorted({round(r['ratio'], 1) for r in ok}):
        for g in grids:
          for Lxv in sorted({r['Lx'] for r in ok}):
            for w in sorted({r['win'] for r in ok}):
                sel = [r for r in ok if round(r['ratio'], 1) == k
                       and (r['nx'], r['nz']) == g and r['win'] == w
                       and r['Lx'] == Lxv]
                if len(sel) < 2:
                    continue
                sp = ", ".join(f"{r['lsl']}/{r['taus']}:{r['meas']-theory(k):+.2f}"
                               for r in sorted(sel, key=lambda z: z['lsl']))
                spread = max(r['meas'] for r in sel) - min(r['meas'] for r in sel)
                worst = max(worst, spread)
                mark = "   <- SPONGE-DEPENDENT" if spread > 0.3 else ""
                print(f"   w/N={k} Lx{Lxv:.0f} {g[0]}x{g[1]} win[{w}]  {sp}   "
                      f"spread {spread:.2f} deg{mark}")
    if worst > 0.3:
        print(f"\n   Worst spread {worst:.2f} deg from a parameter that should not")
        print("   affect the answer.  Until that is under ~0.2 deg, the "
              "refinement\n   series is not measuring the discretization.")


if __name__ == '__main__':
    main()
