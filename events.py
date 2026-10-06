#!/usr/bin/env python3
"""events.py -- how asymmetric is each cold event, in numbers?

  python events.py mry-fig10/walter_sigma_a5_mooring.npz
  python events.py a.npz b.npz          (several records side by side)

Walter et al. (2012) distinguish two shapes by where the fast change is:

  canonical      an abrupt cold front, then warming over hours
  non-canonical  a sharp initial drop, cooling that continues for hours,
                 then an ABRUPT warm front (their Figure 2b: >= 1 degC
                 warming in about five minutes)

Reading that off a plot is subjective. For each event at the lowest sensor this
reports the time from onset to the minimum (cooling), from the minimum back to
ambient (warming), and the fastest cooling and warming rates. The asymmetry is
then a number: warming much faster than cooling is non-canonical, the reverse
canonical.

An event is where the anomaly falls below --thresh of the record's largest
excursion. Rates are taken from a lightly smoothed series so single-sample
spikes do not dominate.
"""
import argparse

import numpy as np


def analyse(path, thresh, smooth_s, alpha, T0):
    d = np.load(path)
    t = d['t']
    temp = T0 + d['b'][:, 0] / (9.80665 * alpha)       # lowest sensor
    # a run that diverged ends in NaN; one bad sample would poison the minimum
    # and every threshold derived from it, so keep only what came before
    # Cut at the first sample that is non-finite or physically impossible: a
    # diverging run grows through enormous finite values before it turns NaN,
    # and those would swamp every threshold derived from the record.
    bad = ~np.isfinite(temp) | (np.abs(temp - temp[0]) > 20.0)
    if bad.any():
        k = int(np.argmax(bad))
        t = t[:k]
        temp = temp[:k]
    T = float(d['T'])
    dt = float(np.median(np.diff(t)))
    k = max(int(round(smooth_s / dt)), 1)
    ker = np.ones(k) / k
    # pad with the edge values: zero-padding makes the first sample read as a
    # deep cold event and every real one falls below the threshold
    padded = np.pad(temp, (k // 2, k - 1 - k // 2), mode='edge')
    ts = np.convolve(padded, ker, mode='valid')
    base = np.median(ts[int(0.1 * len(ts)):int(0.3 * len(ts))])
    anom = ts - base
    level = thresh * anom.min()
    inside = anom < level
    edges = np.diff(inside.astype(int))
    starts = list(np.where(edges == 1)[0] + 1)
    ends = list(np.where(edges == -1)[0] + 1)
    if ends and starts and ends[0] < starts[0]:
        ends = ends[1:]
    rate = np.gradient(ts, t) * 60.0                      # degC per minute
    # widen each crossing to where the anomaly first leaves and returns to 10%
    # of its dip; an event whose first sharp drop recovers briefly before the
    # main cooling crosses the threshold twice, so overlapping windows are one
    # event and are merged before anything is measured
    wins = []
    for s0, e0 in zip(starts, ends):
        depth = anom[s0 + int(np.argmin(anom[s0:e0]))]
        lo = s0
        while lo > 0 and anom[lo] < 0.1 * depth:
            lo -= 1
        hi = e0
        while hi < len(anom) - 1 and anom[hi] < 0.1 * depth:
            hi += 1
        if wins and lo <= wins[-1][1]:
            wins[-1] = (wins[-1][0], max(hi, wins[-1][1]))
        else:
            wins.append((lo, hi))
    rows = []
    for lo, hi in wins:
        seg_min = lo + int(np.argmin(anom[lo:hi + 1]))
        depth = anom[seg_min]
        cool_h = (t[seg_min] - t[lo]) / 3600.0
        warm_h = (t[hi] - t[seg_min]) / 3600.0
        rows.append(dict(t=t[seg_min] / T, depth=-depth, cool_h=cool_h, warm_h=warm_h,
                         cool_rate=-rate[lo:seg_min + 1].min(),
                         warm_rate=rate[seg_min:hi + 1].max()))
    return rows, dt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('records', nargs='+')
    ap.add_argument('--thresh', type=float, default=0.3,
                    help='event threshold as a fraction of the deepest excursion')
    ap.add_argument('--smooth', type=float, default=120.0,
                    help='smoothing window in seconds before taking rates')
    ap.add_argument('--alpha', type=float, default=2.0e-4)
    ap.add_argument('--T0', type=float, default=11.0)
    g = ap.parse_args()

    for path in g.records:
        rows, dt = analyse(path, g.thresh, g.smooth, g.alpha, g.T0)
        print(f"\n{path}   (sampling {dt:.1f} s, smoothed over {g.smooth:.0f} s)")
        print(f"   {'t/T':>6s} {'dT':>6s} {'cool h':>7s} {'warm h':>7s} "
              f"{'max cool':>9s} {'max warm':>9s}  {'warm/cool rate':>14s}")
        for r in rows:
            ratio = r['warm_rate'] / max(r['cool_rate'], 1e-12)
            print(f"   {r['t']:6.2f} {r['depth']:6.3f} {r['cool_h']:7.2f} {r['warm_h']:7.2f} "
                  f"{r['cool_rate']:8.4f}/m {r['warm_rate']:8.4f}/m  {ratio:14.2f}")
        if rows:
            # spun-up events: after t/T = 3 if there are any, else the later half
            late = [r for r in rows if r['t'] >= 3.0] or rows[len(rows) // 2:]
            cr = float(np.median([r['cool_h'] for r in late]))
            wr = float(np.median([r['warm_h'] for r in late]))
            rr = float(np.median([r['warm_rate'] / max(r['cool_rate'], 1e-12) for r in late]))
            # Walter et al. separate the regimes by where the TIME goes: canonical
            # events cool fast and warm slowly; non-canonical ones open with a
            # sharp drop, keep cooling for hours, then warm quickly. A sharp first
            # drop makes the fastest cooling rate high in BOTH, so the rate ratio
            # alone misclassifies non-canonical events; the duration ratio does not.
            dr = cr / max(wr, 1e-12)
            if dr > 1.5:
                shape = 'non-canonical (cools for longer than it warms)'
            elif dr < 1 / 1.5:
                shape = 'canonical (warms for longer than it cools)'
            else:
                shape = 'roughly symmetric'
            print(f"   spun-up events ({len(late)}), median: cooling {cr:.2f} h, warming {wr:.2f} h "
                  f"-> duration ratio {dr:.2f}, fastest warm/cool rate {rr:.2f}  -> {shape}")
    print("\n   Walter et al.: bore period 6-20 h; warm front >= 1 degC in ~5 min "
          "(0.2 degC/min)")


if __name__ == '__main__':
    main()
