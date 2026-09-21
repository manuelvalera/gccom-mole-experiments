#!/usr/bin/env python3
"""seiche_report.py -- read the seiche sweep and say what it means.

  python seiche_report.py seiche-logs

Groups the runs by what was varied and reports the observed convergence order
in each. The thing to look for is the SIGN of the order: a positive order
means the error shrinks under refinement (ordinary), a negative order means it
grows (the anomaly this sweep exists to explain).
"""
import os, re, sys, glob
import numpy as np


def scan(path):
    # utf-8-sig: PowerShell's Out-File -Encoding utf8 writes a BOM, which
    # defeats an anchored ^CMD match and made every run parse as None.
    txt = open(path, encoding='utf-8-sig', errors='ignore').read().replace('\r', '')
    m = re.search(r'relative error = ([+-][0-9.eE+-]+)', txt)
    if not m:
        return None
    cmd = re.search(r'^\s*CMD (.*)$', txt, re.M)
    cmd = cmd.group(1) if cmd else ''

    def opt(flag, default=None, cast=str):
        mm = re.search(r'--%s\s+(\S+)' % flag, cmd)
        return cast(mm.group(1)) if mm else default

    ms = re.search(r'--seiche\s+(\d+)\s+(\d+)', cmd)
    g = re.search(r'grid (\d+)x(\d+)\s+dx=([0-9.]+) m\s+dz=([0-9.]+) m', txt)
    if not (ms and g):
        return None
    return dict(
        I=int(ms.group(1)), J=int(ms.group(2)),
        nx=int(g.group(1)), nz=int(g.group(2)),
        dx=float(g.group(3)), dz=float(g.group(4)),
        err=float(m.group(1)),
        alpha=float(opt('alpha', 1e-4, float)),
        order=int(opt('order', 2, int)),
        bedmode=opt('bedmode', 'constraint'),
        spp=int(opt('spp', 120, int)),
        log=os.path.basename(path))


def order_of(hs, es):
    """Observed order p in err ~ h^p.  Negative p means the error GROWS."""
    hs, es = np.asarray(hs, float), np.asarray(es, float)
    if len(hs) < 2 or np.any(es == 0):
        return float('nan')
    return float(np.polyfit(np.log(hs), np.log(np.abs(es)), 1)[0])


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else 'seiche-logs'
    rows = [r for r in (scan(f) for f in sorted(glob.glob(os.path.join(d, '*.log')))) if r]
    if not rows:
        raise SystemExit(f"no completed seiche runs in {d}/")
    # Sections overlap deliberately (the same baseline config appears in more
    # than one), so collapse exact repeats before anything is averaged or
    # fitted -- otherwise a duplicated point silently gets double weight.
    seen, uniq = set(), []
    for r in rows:
        key = (r['I'], r['J'], r['nx'], r['nz'], r['alpha'], r['order'],
               r['bedmode'], round(r['err'], 12))
        if key in seen:
            continue
        seen.add(key); uniq.append(r)
    ndup = len(rows) - len(uniq)
    rows = uniq
    print(f"\n{len(rows)} distinct seiche runs in {d}/"
          + (f"  ({ndup} duplicate(s) collapsed)" if ndup else "") + "\n")

    print("%-7s %-10s %-9s %-7s %-8s %-12s %-9s" %
          ("mode", "grid", "dz (m)", "alpha", "k", "bedmode", "rel err"))
    print("-" * 74)
    for r in sorted(rows, key=lambda z: (z['I'], z['J'], z['nz'], z['alpha'], z['order'])):
        print("%-7s %-10s %-9.2f %-7.0e %-8d %-12s %+.4e" %
              (f"({r['I']},{r['J']})", f"{r['nx']}x{r['nz']}", r['dz'],
               r['alpha'], r['order'], r['bedmode'], r['err']))

    # ---- mode scan at fixed grid ------------------------------------
    base = [r for r in rows if r['alpha'] == 1e-4 and r['order'] == 2
            and r['bedmode'] == 'constraint']
    fixed = {}
    for r in base:
        fixed.setdefault((r['nx'], r['nz']), []).append(r)
    for g, rs in sorted(fixed.items()):
        if len({(r['I'], r['J']) for r in rs}) < 3:
            continue
        print(f"\nMODE SCAN at {g[0]}x{g[1]}  (p = J*pi/D0, so p*dz grows with J)")
        print("   %-8s %-10s %-12s %-12s" % ("mode", "p*dz", "rel err", "err/(p*dz)^2"))
        for r in sorted(rs, key=lambda z: (z['J'], z['I'])):
            pdz = (r['J']*np.pi/1000.0)*r['dz']
            print("   %-8s %-10.4f %+-12.4e %-12.4g" %
                  (f"({r['I']},{r['J']})", pdz, r['err'], r['err']/pdz**2))
        print("   If the last column is roughly constant the error is a")
        print("   vertical-symbol effect; if it varies wildly it is not.")

    # ---- refinement series -------------------------------------------
    series = {}
    for r in base:
        series.setdefault((r['I'], r['J']), []).append(r)
    for mode, rs in sorted(series.items()):
        rs = sorted(rs, key=lambda z: -z['dz'])      # coarse -> fine
        if len(rs) < 3:
            continue
        p = order_of([r['dz'] for r in rs], [r['err'] for r in rs])
        tag = "GROWS under refinement" if p < 0 else "shrinks (normal)"
        print(f"\nREFINEMENT, mode {mode}, alpha fixed at 1e-4")
        print("   " + "  ".join(f"{r['nx']}x{r['nz']}:{r['err']:+.3e}" for r in rs))
        print(f"   observed order p = {p:+.2f}   -> {tag}")

    # ---- alpha ---------------------------------------------------------
    al = [r for r in rows if r['order'] == 2 and r['bedmode'] == 'constraint'
          and (r['I'], r['J']) == (4, 2)]
    byg = {}
    for r in al:
        byg.setdefault((r['nx'], r['nz']), []).append(r)
    for g, rs in sorted(byg.items()):
        if len({r['alpha'] for r in rs}) < 3:
            continue
        rs = sorted(rs, key=lambda z: z['alpha'])
        spread = max(r['err'] for r in rs) - min(r['err'] for r in rs)
        print(f"\nALPHA SCAN at {g[0]}x{g[1]}")
        for r in rs:
            print("   alpha %-8.0e  %+.4e" % (r['alpha'], r['err']))
        print(f"   spread {spread:.3e}")
        if spread > 0.2*abs(np.mean([r['err'] for r in rs])):
            print("   -> the answer DEPENDS ON ALPHA. The lateral Robin rows are")
            print("      implicated: with no sponge they are the only thing")
            print("      holding the side walls, and alpha is fixed while the")
            print("      discrete normal derivative scales like 1/dz.")
        else:
            print("   -> alpha is not the driver; look elsewhere.")

    # ---- alpha scaled with 1/dz ---------------------------------------
    scaled = [r for r in rows if r['order'] == 2 and (r['I'], r['J']) == (4, 2)
              and abs(r['alpha']*r['dz'] - 1e-4*20.0) < 0.35*1e-4*20.0
              and r['alpha'] != 1e-4]
    if len(scaled) >= 3:
        scaled = sorted(scaled, key=lambda z: -z['dz'])   # coarse -> fine
        p = order_of([r['dz'] for r in scaled], [r['err'] for r in scaled])
        print("\nREFINEMENT with alpha SCALED as 1/dz")
        print("   " + "  ".join(f"{r['nx']}x{r['nz']}:{r['err']:+.3e}" for r in scaled))
        print(f"   observed order p = {p:+.2f}")
        print("   Compare with the fixed-alpha series above. If this one")
        print("   converges and that one does not, the boundary condition")
        print("   drifting with resolution is the whole story.")

    # ---- order and bedmode --------------------------------------------
    for key, label in (('order', 'OPERATOR ORDER'), ('bedmode', 'BEDMODE')):
        grp = {}
        for r in rows:
            if (r['I'], r['J']) != (4, 2) or r['alpha'] != 1e-4:
                continue
            grp.setdefault((r['nx'], r['nz']), {}).setdefault(r[key], r)
        for g, by in sorted(grp.items()):
            if len(by) < 2:
                continue
            print(f"\n{label} at {g[0]}x{g[1]}")
            for v, r in sorted(by.items(), key=lambda kv: str(kv[0])):
                print("   %-14s %+.4e" % (v, r['err']))

    print("\nReminder: at fixed mode, refining makes k*dz SMALLER, so any")
    print("interpolation-symbol error must shrink. An error that grows is")
    print("converging to the wrong operator, not converging slowly.")


if __name__ == '__main__':
    main()
