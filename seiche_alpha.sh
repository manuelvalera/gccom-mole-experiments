#!/usr/bin/env bash
# seiche_alpha.sh -- see seiche_alpha.ps1 for the full note on why alpha is
# the suspect and why 1e-5 is a zero crossing rather than a fix.
set -u
: "${MOLE_SRC:?set MOLE_SRC first}"
ONLY="${1:-all}"
PAT='SEICHE|BEAM ANGLE|seiche mode|initial field projected|crossings|rms residual|max/mean|consistency|UNSTABLE|DIVERGED|\*\*\*'
mkdir -p alpha-logs; SEC=A; N=0
sec(){ if [ "$ONLY" = all ] || [ "$ONLY" = "$1" ]; then SEC=$1; N=0; echo; echo "=== $1. $2 ==="; return 0; fi; return 1; }
run(){ for a in "$@"; do [ -n "$a" ] || { echo "   !! blank argument"; return; }; done
  N=$((N+1)); log=$(printf "alpha-logs/s%s_%03d.log" "$SEC" $N)
  echo "   python3 iwbcurv.py --mole \$MOLE_SRC --beat 9999 $*"
  echo "CMD iwbcurv.py --mole $MOLE_SRC --beat 9999 $*" > "$log"
  python3 iwbcurv.py --mole "$MOLE_SRC" --beat 9999 "$@" >>"$log" 2>&1; c=$?
  if grep -qE "$PAT" "$log" && [ $c -eq 0 ]; then grep -E "$PAT" "$log" | sed "s/^/   /"
  else echo "   !! no result (exit=$c):"; tail -20 "$log" | sed "s/^/   | /"; fi; }

if sec A "alpha floor"; then
for a in 1e-9 1e-8 1e-7 1e-6 3e-6 1e-5; do echo "-- alpha=$a 256x101"
  run --ab 0 --bulge 0 --Lx 6000 --nx 256 --nz 101 --seiche 4 2 --nper 6 --spp 400 --alpha $a; done; fi

if sec B "TRUE convergence order at small alpha"; then
for g in "128 51" "192 76" "256 101" "384 151"; do set -- $g; echo "-- alpha=1e-6 $1x$2"
  run --ab 0 --bulge 0 --Lx 6000 --nx $1 --nz $2 --seiche 4 2 --nper 6 --spp 400 --alpha 1e-6; done
for g in "128 51" "256 101"; do set -- $g; echo "-- alpha=1e-6 mode (2,1) $1x$2"
  run --ab 0 --bulge 0 --Lx 6000 --nx $1 --nz $2 --seiche 2 1 --nper 6 --spp 400 --alpha 1e-6; done; fi

if sec C "THE ONE THAT MATTERS: do the beam angles move?"; then
for cfg in "0.6 49.87 133 1000" "0.8 37.40 75 562"; do set -- $cfg
  R=$1; TS=$2; W1=$3; W2=$4
  for a in 1e-4 1e-5 1e-6; do echo "-- ratio=$R alpha=$a 512x201"
    run --Lx 6000 --nx 512 --nz 201 --ratio $R --nper 20 --spp 120 \
        --lsl 0.20 --taus $TS --win $W1 $W2 --alpha $a --forcegrid; done; done; fi

if sec D "control: does a small alpha survive a long run?"; then
for a in 1e-4 1e-6; do echo "-- alpha=$a 40 periods"
  run --ab 0 --bulge 0 --Lx 6000 --nx 192 --nz 76 --seiche 4 2 --nper 40 --spp 400 --alpha $a; done; fi

echo; echo "logs in alpha-logs/; now: python3 seiche_report.py alpha-logs"
