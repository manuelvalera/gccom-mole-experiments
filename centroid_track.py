#!/usr/bin/env python3
"""centroid_track.py -- is the beam curvature real, or an artifact of HOW we track it?

  python centroid_track.py refit-npz
  python centroid_track.py refit-npz --halfwidth 200 --thresh 0.3

WHY

curvature.py found the tracked beam bends at w/N=0.8: a sag of about -50 m
across the window at every grid from 256x101 to 640x251. Grid-independence
rules out a discretization artifact, but it does NOT rule out a MEASUREMENT
artifact, because the tracker itself is a convention:

    column-wise argmax of rms speed

picks the peak of the vertical profile at each x. If the beam's cross-profile
is asymmetric and evolves along the beam -- plausible for a beam born against
the bottom boundary, whose lower side is truncated near the source -- then the
peak drifts relative to the beam axis. That drift is a property of the
continuum field, so it is grid-independent too. The sag alone cannot tell the
two apart.

The depth_compare image also ruled out the other explanation I offered:
the reflected beam leaves the surface OUTWARD at |x| ~ 730 m and never enters
the tracked window, so the argmax is not hopping onto a reflection.

THE TEST

Track the beam a different way and see whether the curvature survives:
take slices PERPENDICULAR to the theory direction, and on each slice find the
energy-weighted centroid (weights rms^2, only above a fraction of the slice
maximum so the background does not pull it). That is the beam axis in the
sense an energy-flux argument would use.

    curvature survives with the centroid  -> the beam really bends
    curvature disappears                  -> it was the argmax convention

Either answer is publishable; they lead to different sentences in the paper.
"""
import os, sys, glob, argparse
import numpy as np
from scipy.interpolate import LinearNDInterpolator


def centres(XM, ZM, rms):
    n, m = ZM.shape[0] - 1, ZM.shape[1] - 1
    Xc = 0.25*(XM[:-1, :-1] + XM[1:, :-1] + XM[:-1, 1:] + XM[1:, 1:])
    Zc = 0.25*(ZM[:-1, :-1] + ZM[1:, :-1] + ZM[:-1, 1:] + ZM[1:, 1:])
    R = rms[1:-1, 1:-1]
    return Xc.ravel(), Zc.ravel(), R.ravel()


def argmax_track(XM, ZM, rms):
    n, m = ZM.shape[0] - 1, ZM.shape[1] - 1
    Xc = 0.25*(XM[:-1, :-1] + XM[1:, :-1] + XM[:-1, 1:] + XM[1:, 1:])
    Zc = 0.25*(ZM[:-1, :-1] + ZM[1:, :-1] + ZM[:-1, 1:] + ZM[1:, 1:])
    R = rms[1:-1, 1:-1]
    pts = []
    for j in range(R.shape[1]):
        col = R[:, j]
        if not np.isfinite(col).any() or np.nanmax(col) <= 0:
            continue
        k = int(np.nanargmax(col))
        pts.append((abs(Xc[k, j]), Zc[k, j]))
    return np.array(sorted(pts))


def centroid_track(interp, D0, theta, s_list, halfwidth, thresh, side):
    """Energy centroid on slices perpendicular to the theory direction."""
    t = np.array([side*np.cos(theta), np.sin(theta)])       # along beam
    nrm = np.array([-t[1], t[0]])                             # across beam
    src = np.array([0.0, -D0])
    eta = np.linspace(-halfwidth, halfwidth, 241)
    out = []
    for s in s_list:
        base = src + s*t
        P = base[None, :] + eta[:, None]*nrm[None, :]
        w = interp(P[:, 0], P[:, 1])
        good = np.isfinite(w)
        if good.sum() < 20:
            continue
        w = np.where(good, w, 0.0)**2
        if w.max() <= 0:
            continue
        w = np.where(w >= thresh*w.max(), w, 0.0)
        ec = float((w*eta).sum()/w.sum())
        p = base + ec*nrm
        # reject slices whose centroid sits on the window edge: the beam is
        # further off-axis than the slice reaches, and the value is meaningless
        if abs(ec) > 0.9*halfwidth:
            continue
        out.append((abs(p[0]), p[1]))
    return np.array(out)


def fit_angle(P, lo, hi):
    sel = (P[:, 0] >= lo) & (P[:, 0] <= hi)
    if sel.sum() < 6:
        return None
    x, z = P[sel, 0], P[sel, 1]
    sl = np.polyfit(x, z, 1)[0]
    return np.degrees(np.arctan(abs(sl)))


def sag(P, lo, hi):
    sel = (P[:, 0] >= lo) & (P[:, 0] <= hi)
    if sel.sum() < 10:
        return None, None, None
    x, z = P[sel, 0], P[sel, 1]
    c2 = np.polyfit(x, z, 2)
    c1 = np.polyfit(x, z, 1)
    r1 = np.sqrt(((z - np.polyval(c1, x))**2).mean())
    r2 = np.sqrt(((z - np.polyval(c2, x))**2).mean())
    return c2[0]*(hi - lo)**2, r1, r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir', nargs='?', default='refit-npz')
    ap.add_argument('--halfwidth', type=float, default=150.0,
                    help='half-width of each perpendicular slice, m')
    ap.add_argument('--thresh', type=float, default=0.5,
                    help='only rms^2 above this fraction of the slice max counts')
    g = ap.parse_args()

    files = sorted(glob.glob(os.path.join(g.dir, '*.npz')))
    if not files:
        raise SystemExit(f"no .npz in {g.dir}/")

    rows = {}
    for f in files:
        z = np.load(f)
        ratio = float(z['ratio']); D0 = float(np.abs(z['Z']).max())
        nz, nx = z['Z'].shape
        theta = np.arcsin(ratio); th_deg = np.degrees(theta)
        bounce = D0/np.tan(theta)
        x, zz, r = centres(z['X'], z['Z'], z['rms'])
        ok = np.isfinite(r)
        interp = LinearNDInterpolator(np.c_[x[ok], zz[ok]], r[ok])

        # along-beam distances covering 0.05-0.90 of a bounce horizontally
        s_list = np.linspace(0.05, 0.90, 170)*bounce/np.cos(theta)
        C = np.vstack([centroid_track(interp, D0, theta, s_list,
                                      g.halfwidth, g.thresh, sd)
                       for sd in (+1, -1)])
        A = argmax_track(z['X'], z['Z'], z['rms'])
        rows.setdefault((ratio, round(D0)), []).append(
            dict(grid=f"{nx-1}x{nz-1}", dz=D0/(nz-1), A=A, C=C,
                 th=th_deg, bounce=bounce))

    for (ratio, D0r), rs in sorted(rows.items()):
        rs = sorted(rs, key=lambda q: -q['dz'])
        th, b = rs[0]['th'], rs[0]['bounce']
        print("=" * 78)
        print(f"w/N = {ratio}   D0 = {D0r} m   theory {th:.2f} deg   bounce {b:.0f} m")
        print(f"slice half-width {g.halfwidth:.0f} m, threshold {g.thresh:.2f} of slice max")
        print("=" * 78)

        edges = np.array([0.10, 0.25, 0.40, 0.55, 0.70, 0.85])*b
        print("\nLOCAL angle bias by distance band (deg)")
        print("%-10s %-9s %s" % ("grid", "tracker",
              "  ".join("%4.2f-%4.2f" % (edges[i]/b, edges[i+1]/b)
                        for i in range(len(edges)-1))))
        for q in rs:
            for name, P in (("argmax", q['A']), ("centroid", q['C'])):
                cells = []
                for i in range(len(edges)-1):
                    a = fit_angle(P, edges[i], edges[i+1])
                    cells.append("   n/a   " if a is None else "%+8.2f " % (a - th))
                print("%-10s %-9s %s" % (q['grid'], name, " ".join(cells)))

        print("\nSAG across 0.10-0.75 of a bounce (m), and fit residuals")
        print("%-10s %-9s %-10s %-9s %-9s" % ("grid", "tracker", "sag", "lin rms", "quad rms"))
        for q in rs:
            for name, P in (("argmax", q['A']), ("centroid", q['C'])):
                sg, r1, r2 = sag(P, 0.10*b, 0.75*b)
                if sg is None:
                    print("%-10s %-9s insufficient points" % (q['grid'], name)); continue
                print("%-10s %-9s %+-10.1f %-9.2f %-9.2f" % (q['grid'], name, sg, r1, r2))

        print("\nWHOLE-WINDOW angle bias over 0.10-0.40 and 0.10-0.75 of a bounce")
        for q in rs:
            for name, P in (("argmax", q['A']), ("centroid", q['C'])):
                a1 = fit_angle(P, 0.10*b, 0.40*b); a2 = fit_angle(P, 0.10*b, 0.75*b)
                f = lambda a: "  n/a " if a is None else "%+6.2f" % (a - th)
                print("   %-10s %-9s  0.10-0.40: %s    0.10-0.75: %s"
                      % (q['grid'], name, f(a1), f(a2)))
        # extrapolation of the centroid series, so the output is self-contained
        if len(rs) >= 3:
            print("\nEXTRAPOLATED LIMIT of the centroid bias (e = A*dz^p + e0)")
            print("   %-10s %-17s %-17s %-18s" % ("window", "p=1 (resid)", "p=2 (resid)", "free p"))
            for lo, hi in ((0.10, 0.40), (0.10, 0.75)):
                hs, es = [], []
                for q in rs:
                    a = fit_angle(q['C'], lo*b, hi*b)
                    if a is not None:
                        hs.append(q['dz']); es.append(a - th)
                if len(hs) < 3:
                    continue
                hs, es = np.array(hs), np.array(es)
                cells = []
                for p in (1.0, 2.0):
                    M = np.vstack([hs**p, np.ones(len(hs))]).T
                    (A, e0), *_ = np.linalg.lstsq(M, es, rcond=None)
                    cells.append("%+6.2f (%.3f)" % (e0, np.abs(es - (A*hs**p + e0)).max()))
                best = None
                for p in np.linspace(0.3, 3.0, 271):
                    M = np.vstack([hs**p, np.ones(len(hs))]).T
                    (A, e0), *_ = np.linalg.lstsq(M, es, rcond=None)
                    ss = ((es - (A*hs**p + e0))**2).sum()
                    if best is None or ss < best[0]:
                        best = (ss, p, e0)
                print("   %-10s %-17s %-17s p=%.2f e0=%+.2f"
                      % ("%.2f-%.2f" % (lo, hi), cells[0], cells[1], best[1], best[2]))
            print("   (%d grids; with fewer than 5 the free-p fit is not meaningful)" % len(rs))

        print("\n   If the centroid rows are FLAT across bands while argmax drifts,")
        print("   the curvature was the tracker. If both drift the same way, the")
        print("   beam bends.\n")


if __name__ == '__main__':
    main()
