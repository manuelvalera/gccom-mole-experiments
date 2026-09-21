#!/usr/bin/env bash
# seiche_study.sh -- see seiche_study.ps1 for the full note on what is being
# chased and why. Same sections, same logs.
#   export MOLE_SRC=/path/to/mole/src/octave
#   ./seiche_study.sh        # all
#   ./seiche_study.sh 3      # one section
#   python3 seiche_report.py seiche-logs
set -u
: "${MOLE_SRC:?set MOLE_SRC first}"
ONLY="${1:-all}"
PAT='SEICHE|seiche mode|initial field projected|grid .*dx=|crossings|UNSTABLE|DIVERGED|\*\*\*'
mkdir -p seiche-logs; SEC=0; N=0
sec(){ if [ "$ONLY" = all ] || [ "$ONLY" = "$1" ]; then SEC=$1; N=0; echo; echo "=== $1. $2 ==="; return 0; fi; return 1; }
run(){ for a in "$@"; do [ -n "$a" ] || { echo "   !! blank argument, not running"; return; }; done
  N=$((N+1)); log=$(printf "seiche-logs/s%s_%03d.log" "$SEC" $N)
  echo "   python3 iwbcurv.py --mole \$MOLE_SRC --beat 9999 --ab 0 --bulge 0 $*"
  echo "CMD iwbcurv.py --mole $MOLE_SRC --beat 9999 --ab 0 --bulge 0 $*" > "$log"
  python3 iwbcurv.py --mole "$MOLE_SRC" --beat 9999 --ab 0 --bulge 0 "$@" >>"$log" 2>&1; c=$?
  if grep -qE "$PAT" "$log" && [ $c -eq 0 ]; then grep -E "$PAT" "$log" | sed "s/^/   /"
  else echo "   !! no result (exit=$c), tail of $log:"; tail -20 "$log" | sed "s/^/   | /"; fi; }

echo "=== 0. preflight ==="
[ -f iwbcurv.py ] || { echo "MISSING iwbcurv.py"; exit 1; }
[ -f "$MOLE_SRC/gridGen.m" ] || { echo "MOLE_SRC has no gridGen.m"; exit 1; }
run --Lx 6000 --nx 128 --nz 51 --seiche 4 2 --nper 6 --spp 400
[ "$ONLY" = 0 ] && exit 0

if sec 1 "mode scan at FIXED grid"; then
for m in "2 1" "4 1" "8 1" "2 2" "4 2" "2 4" "4 4" "2 8"; do set -- $m
  echo "-- mode ($1,$2)"
  run --Lx 6000 --nx 128 --nz 51 --seiche $1 $2 --nper 6 --spp 400; done; fi

if sec 2 "refinement at fixed mode"; then
for m in "4 2" "2 1"; do set -- $m; I=$1; J=$2
  for g in "128 51" "192 76" "256 101" "384 151"; do set -- $g
    echo "-- mode ($I,$J)  $1x$2"
    run --Lx 6000 --nx $1 --nz $2 --seiche $I $J --nper 6 --spp 400; done; done; fi

if sec 3 "THE SUSPECT: alpha, fixed and grid-scaled"; then
for a in 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1; do echo "-- alpha=$a  256x101"
  run --Lx 6000 --nx 256 --nz 101 --seiche 4 2 --nper 6 --spp 400 --alpha $a; done
for c in "128 51 1.0e-4" "192 76 1.5e-4" "256 101 2.0e-4" "384 151 3.0e-4"; do set -- $c
  echo "-- alpha scaled with 1/dz: $1x$2 alpha=$3"
  run --Lx 6000 --nx $1 --nz $2 --seiche 4 2 --nper 6 --spp 400 --alpha $3; done; fi

if sec 4 "operator order on the seiche"; then
for k in 2 4 6; do for g in "128 51" "256 101"; do set -- $g
  echo "-- k=$k  $1x$2"
  run --Lx 6000 --nx $1 --nz $2 --seiche 4 2 --nper 6 --spp 400 --order $k; done; done; fi

if sec 5 "separate dx from dz"; then
echo "-- dz fixed at nz=51, refine nx"
for nx in 96 128 192 256; do run --Lx 6000 --nx $nx --nz 51 --seiche 4 2 --nper 6 --spp 400; done
echo "-- dx fixed at nx=128, refine nz"
for nz in 41 51 76 101; do run --Lx 6000 --nx 128 --nz $nz --seiche 4 2 --nper 6 --spp 400; done; fi

if sec 6 "control: is the bed constraint involved?"; then
for bm in constraint post none; do echo "-- bedmode=$bm"
  run --Lx 6000 --nx 256 --nz 101 --seiche 4 2 --nper 6 --spp 400 --bedmode $bm; done; fi

echo; echo "logs in seiche-logs/; now: python3 seiche_report.py seiche-logs"
