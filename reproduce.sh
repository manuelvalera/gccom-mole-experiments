#!/usr/bin/env bash
# reproduce.sh -- same as reproduce.ps1, for the Linux box.
#   export MOLE_SRC=/path/to/mole/src/matlab_octave
#   ./reproduce.sh          # everything
#   ./reproduce.sh 4        # just section 4
#
# Reference environment:
#   MOLE commit b0b9ff76aff40481036d30b63a3bb00243d921ea (2026-09-03)
#   Octave 8.4.0, Python 3.12.3, numpy 2.4.4, scipy 1.17.1, SuperLU (splu)
set -u
: "${MOLE_SRC:?set MOLE_SRC first}"
ONLY="${1:-all}"
sec(){ [ "$ONLY" = all ] || [ "$ONLY" = "$1" ]; }
PAT='r_x .*max|sponge:|consistency:|rms residual|max/mean|BEAM ANGLE|DIVERGED|static projection|VIOLATED|UNSTABLE|Raise --taus'
mkdir -p repro-logs; N=0
run(){ N=$((N+1)); log=$(printf "repro-logs/run_%03d.log" $N)
  echo "   python iwbcurv.py --mole \$MOLE_SRC --beat 9999 $*"
  python3 iwbcurv.py --mole "$MOLE_SRC" --beat 9999 "$@" >"$log" 2>&1; code=$?
  if grep -qE "$PAT" "$log" && [ $code -eq 0 ]; then grep -E "$PAT" "$log" | sed 's/^/   /'
  else echo "   !! no summary lines (exit=$code) -- tail of $log :"; tail -25 "$log" | sed 's/^/   | /'; fi; }

# preflight
for f in iwbcurv.py test_bedrows.py synth_tracker.py; do
  [ -f "$f" ] || { echo "MISSING $f"; exit 1; }; done
[ -d grids ] || { echo "MISSING grids/ -- run from your iwb-curv dir"; exit 1; }
[ -f "$MOLE_SRC/gridGen.m" ] || { echo "MOLE_SRC has no gridGen.m: $MOLE_SRC"; exit 1; }
echo "preflight OK"

if sec 1; then echo; echo "=== 1. bed constraint: acceptance test, 40 periods, no --bsc ==="
for r in 0.2 0.4 0.6 0.8; do echo "-- omega/N = $r"
  run --ratio $r --Lx 6000 --nx 256 --nz 101 --nper 40; done; fi

if sec 2; then echo; echo "=== 2. A/B against the old penalty, and the controls ==="
echo "-- old penalty (expect DIVERGED at t/T = 9.00)"
run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 40 --bedmode projection --bsc 12
echo "-- constraint, same grid"
run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 40
for m in post none; do echo "-- control: $m"
  run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 5 --bedmode $m; done
echo "-- equilibration invariance (must match --bedeq auto exactly)"
run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 40 --bedeq unit; fi

if sec 3; then echo; echo "=== 3. sponge: the contamination fix (worst case, omega/N = 0.8) ==="
for c in "0.10 100" "0.20 100" "0.20 50" "0.30 50" "0.30 25"; do set -- $c
  echo "-- lsl=$1 taus=$2"
  run --ratio 0.8 --Lx 6000 --nx 256 --nz 101 --nper 20 --lsl $1 --taus $2; done
echo "-- sponge stability guard (must refuse, not blow up)"
run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 2 --lsl 0.2 --taus 5; fi

if sec 4; then echo; echo "=== 4. corrected acceptance table, nz=101 and nz=201 ==="
# ratio Lx nx taus=T/30 spp ; 0.2 needs the bigger box (bounce = 4900 m)
for nz in 101 201; do
for c in "0.2 12000 512 149.6 60" "0.4 6000 256 74.8 120" \
         "0.6 6000 256 49.9 120"  "0.8 6000 256 37.4 120"; do set -- $c
  echo "-- omega/N = $1  nz=$nz"
  run --ratio $1 --Lx $2 --nx $3 --nz $nz --nper 20 --spp $5 \
      --lsl 0.2 --taus $4 --ab 20 --forcegrid; done; done; fi

if sec 5; then echo; echo "=== 5. RETRACTION control: vary r_x 10x at FIXED dz -- error must stay flat ==="
for a in 20 10 5 2; do echo "-- ab = $a m"
  run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 20 --spp 120 \
      --lsl 0.2 --taus 49.9 --ab $a; done; fi

if sec 6; then echo; echo "=== 6. convergence: refine dx and dz together (ab=10, r_x never binds) ==="
for c in "256 101" "384 151" "512 201"; do set -- $c
  echo "-- nx=$1 nz=$2"
  run --ratio 0.6 --Lx 6000 --nx $1 --nz $2 --ab 10 --nper 20 --spp 60 \
      --lsl 0.2 --taus 49.9; done; fi

if sec 7; then echo; echo "=== 7. does r_x > 1 break it? benchmark ridge, --forcegrid ==="
for nz in 101 151 201; do echo "-- nz=$nz"
  run --ratio 0.6 --Lx 6000 --nx 256 --nz $nz --ab 20 --nper 20 --spp 120 \
      --lsl 0.2 --taus 49.9 --forcegrid; done; fi

if sec 8; then echo; echo "=== 8. dead ends re-tested on THIS solver (should all be flat) ==="
echo "-- timestep: angle must be identical at spp 60/120/240"
for s in 60 120 240; do run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 20 --spp $s; done
echo "-- operator order: k = 2/4/6"
for k in 2 4 6; do run --ratio 0.6 --Lx 6000 --nx 256 --nz 101 --nper 20 --order $k; done; fi

if sec 9; then echo; echo "=== 9. offline tests (no Octave, no MOLE) ==="
python3 test_bedrows.py; python3 synth_tracker.py "$MOLE_SRC"; fi

if sec 10; then echo; echo "=== 10. save the best configuration for view.py ==="
run --ratio 0.6 --Lx 6000 --nx 256 --nz 201 --nper 20 --spp 120 \
    --lsl 0.2 --taus 49.9 --forcegrid --out iwb_r06_refined.npz
echo "then:  python3 view.py iwb_r06_refined.npz"; fi
