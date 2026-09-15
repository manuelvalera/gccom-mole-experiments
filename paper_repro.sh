#!/usr/bin/env bash
# paper_repro.sh -- see paper_repro.ps1 for the full note on what matches
# Garcia et al. (2019) sec 3.3 and what does not.  Same sections.
#   export MOLE_SRC=/path/to/mole/src/matlab_octave
#   ./paper_repro.sh            # all
#   ./paper_repro.sh 2          # one section
#   python3 compare_paper.py paper-logs
set -u
: "${MOLE_SRC:?set MOLE_SRC first}"
ONLY="${1:-all}"
PAT='r_x .*max|sponge:|wander|consistency:|rms residual|max/mean|BEAM ANGLE|DIVERGED|static projection|VIOLATED|UNSTABLE|Raise --taus'
mkdir -p paper-logs; N=0; SEC=0
# section-prefixed log names: re-running one section must not clobber the rest
sec(){ if [ "$ONLY" = all ] || [ "$ONLY" = "$1" ]; then SEC=$1; N=0; return 0; fi; return 1; }
run(){ for a in "$@"; do [ -n "$a" ] || { echo "   !! empty argument -- a computed value came back blank; not running"; return; }; done
  N=$((N+1)); log=$(printf "paper-logs/s%s_%03d.log" "$SEC" $N)
  echo "   python3 iwbcurv.py --mole \$MOLE_SRC --beat 9999 $*"
  python3 iwbcurv.py --mole "$MOLE_SRC" --beat 9999 "$@" >"$log" 2>&1; code=$?
  if grep -qE "$PAT" "$log" && [ $code -eq 0 ]; then grep -E "$PAT" "$log" | sed 's/^/   /'
  else echo "   !! no summary lines (exit=$code) -- tail of $log :"; tail -25 "$log" | sed 's/^/   | /'; fi; }
# spp such that dt/taus <= 0.4 (the sponge relaxation is explicit)
spp(){ python3 -c "import math,sys;r=float(sys.argv[1]);t=float(sys.argv[2]);print(max(120,math.ceil(2*math.pi/(r*0.007)/(0.4*t))))" "$1" "$2"; }

echo "=== 0. preflight ==="
for f in iwbcurv.py view.py compare_paper.py; do [ -f "$f" ] || { echo "MISSING $f"; exit 1; }; done
[ -f grids/iwbridge/left.m ] || { echo "MISSING grids/iwbridge/ -- run from your iwb-curv dir"; exit 1; }
grep -qE 'bulge *= *0*\.?0*\*D0' grids/iwbridge/left.m && { echo "left.m has bulge = 0 -- NOT the curvilinear grid"; exit 1; }
grep -oE 'bulge *= *[0-9.]+\*D0' grids/iwbridge/left.m | sed 's/^/   /'
[ -f "$MOLE_SRC/gridGen.m" ] || { echo "MOLE_SRC has no gridGen.m: $MOLE_SRC"; exit 1; }
echo "   -- grid smoke test (must say FULLY CURVILINEAR) --"
run --Lx 3000 --nx 128 --nz 101 --ratio 0.6 --gridonly
[ "$ONLY" = 0 ] && exit 0

if sec 1; then echo; echo "=== 1. sponge: paper's Lsl=Lx/10,taus=100 vs alternatives (~25 min) ==="
for c in "0.10 100" "0.20 50" "0.30 30"; do set -- $c; L=$1; TS=$2
  for r in 0.2 0.4 0.6 0.8; do S=$(spp $r $TS)
    echo "-- lsl=$L taus=$TS  w/N=$r  spp=$S"
    run --Lx 3000 --nx 128 --nz 101 --ratio $r --nper 20 --spp $S \
        --lsl $L --taus $TS --win 200 500; done; done; fi

if sec 2; then echo; echo "=== 2. THE COMPARISON: paper config, calibrated sponge (~10 min) ==="
for r in 0.2 0.4 0.6 0.8; do S=$(spp $r 30)
  echo "-- w/N=$r  spp=$S"
  run --Lx 3000 --nx 128 --nz 101 --ratio $r --nper 20 --spp $S \
      --lsl 0.30 --taus 30 --win 200 500 --out paper-logs/pf$r.npz; done; fi

if sec 3; then echo; echo "=== 3. convergence: refine the paper config (~60 min) ==="
for c in "192 151" "256 201"; do set -- $c; NX=$1; NZ=$2
  for r in 0.2 0.4 0.6 0.8; do S=$(spp $r 30)
    echo "-- w/N=$r  nx=$NX nz=$NZ  spp=$S"
    run --Lx 3000 --nx $NX --nz $NZ --ratio $r --nper 20 --spp $S \
        --lsl 0.30 --taus 30 --win 200 500 --forcegrid; done; done; fi

if sec 4; then echo; echo "=== 4. Fig 8 reproduction ==="
echo "   paper Fig 8 is TOTAL u/u0 at t=19.5T (ubc=0), colour limit +/-0.15"
python3 view.py paper-logs/pf0.2.npz paper-logs/pf0.4.npz \
                paper-logs/pf0.6.npz paper-logs/pf0.8.npz \
                --total --clim 0.15 -o fig8_repro.png; fi

if sec 5; then echo; echo "=== 5. comparison table vs Fig 9 ==="
python3 compare_paper.py paper-logs; fi

if sec 6; then echo; echo "=== 6. WHY refinement diverges: window and sponge sensitivity at 256x201 (~35 min) ==="
# See paper_repro.ps1 section 6 for the full note.  (a) the paper's [200,500]
# window may sit in the near field; (b) Lsl=900 m may reach the window.
for r in 0.4 0.8; do S=$(spp $r 30)
  for w in "200 500" "300 700" "400 900"; do set -- $w
    echo "-- (a) window w/N=$r  [$1,$2]"
    run --Lx 3000 --nx 256 --nz 201 --ratio $r --nper 20 --spp $S \
        --lsl 0.30 --taus 30 --win $1 $2 --forcegrid; done
  echo "-- (a) window w/N=$r  auto"
  run --Lx 3000 --nx 256 --nz 201 --ratio $r --nper 20 --spp $S \
      --lsl 0.30 --taus 30 --forcegrid
  for c in "0.20 50" "0.40 30"; do set -- $c
    echo "-- (b) sponge w/N=$r  lsl=$1 taus=$2"
    run --Lx 3000 --nx 256 --nz 201 --ratio $r --nper 20 --spp $S \
        --lsl $1 --taus $2 --win 200 500 --forcegrid; done; done; fi

if sec 7; then echo; echo "=== 7. is the divergence a 3 km artifact? same test at Lx = 6 km (~90 min) ==="
# See paper_repro.ps1 section 7.  At 3 km the sponge needed to kill reflections
# reaches the fit window; at 6 km it does not.  dx,dz held to the 3 km values.
for c in "256 101" "384 151" "512 201"; do set -- $c; NX=$1; NZ=$2
  for r in 0.2 0.4 0.6 0.8; do
    TS=$(python3 -c "import math,sys;print(round(2*math.pi/(float(sys.argv[1])*0.007)/30,1))" $r)
    S=$(spp $r $TS)
    echo "-- w/N=$r  nx=$NX nz=$NZ  lsl=0.20 taus=$TS  spp=$S"
    run --Lx 6000 --nx $NX --nz $NZ --ratio $r --nper 20 --spp $S \
        --lsl 0.20 --taus $TS --win 200 500 --forcegrid; done; done
echo "-- sponge independence check at the finest grid"
for r in 0.4 0.8; do for l in 0.15 0.25; do
  TS=$(python3 -c "import math,sys;print(round(2*math.pi/(float(sys.argv[1])*0.007)/30,1))" $r)
  S=$(spp $r $TS)
  echo "-- w/N=$r  lsl=$l"
  run --Lx 6000 --nx 512 --nz 201 --ratio $r --nper 20 --spp $S \
      --lsl $l --taus $TS --win 200 500 --forcegrid; done; done; fi

if sec 8; then echo; echo "=== 8. window sweep at 6 km with a CLEAN sponge, windows as bounce fractions (~30 min) ==="
# See paper_repro.ps1 section 8.  Windows are fractions of the bounce so that
# different frequencies are measured at the same place along their own beams.
for r in 0.4 0.6 0.8; do
  TS=$(python3 -c "import math,sys;print(round(2*math.pi/(float(sys.argv[1])*0.007)/30,1))" $r)
  S=$(spp $r $TS)
  B=$(python3 -c "import math,sys;r=float(sys.argv[1]);print(round(1000/math.tan(math.asin(r))))" $r)
  for f in "0.10 0.30" "0.20 0.45" "0.30 0.60" "0.45 0.75"; do set -- $f
    W0=$(python3 -c "print(round($B*$1))"); W1=$(python3 -c "print(round($B*$2))")
    if [ "$W1" -gt 1800 ]; then echo "-- w/N=$r frac $1-$2 -> [$W0,$W1] m  SKIPPED (reaches the sponge)"; continue; fi
    echo "-- w/N=$r  bounce $B m  frac $1-$2  -> win [$W0,$W1] m"
    run --Lx 6000 --nx 512 --nz 201 --ratio $r --nper 20 --spp $S \
        --lsl 0.20 --taus $TS --win $W0 $W1 --forcegrid; done; done; fi

if sec 9; then echo; echo "=== 9. is the beam curvature grid-converged? (~25 min) ==="
# See paper_repro.ps1 section 9.  Same bounce fractions at 384x151; 512x201 is
# already in section 8.  Spread constant => physical; shrinking => numerical.
for r in 0.6 0.8; do
  TS=$(python3 -c "import math,sys;print(round(2*math.pi/(float(sys.argv[1])*0.007)/30,1))" $r)
  S=$(spp $r $TS)
  B=$(python3 -c "import math,sys;r=float(sys.argv[1]);print(round(1000/math.tan(math.asin(r))))" $r)
  for f in "0.10 0.30" "0.20 0.45" "0.30 0.60" "0.45 0.75"; do set -- $f
    W0=$(python3 -c "print(round($B*$1))"); W1=$(python3 -c "print(round($B*$2))")
    echo "-- w/N=$r  384x151  frac $1-$2  -> win [$W0,$W1] m"
    run --Lx 6000 --nx 384 --nz 151 --ratio $r --nper 20 --spp $S \
        --lsl 0.20 --taus $TS --win $W0 $W1 --forcegrid; done
  W0=$(python3 -c "print(round($B*0.10))"); W1=$(python3 -c "print(round($B*0.75))")
  echo "-- w/N=$r  512x201  full first bounce [$W0,$W1] -> npz"
  run --Lx 6000 --nx 512 --nz 201 --ratio $r --nper 20 --spp $S \
      --lsl 0.20 --taus $TS --win $W0 $W1 --forcegrid --out paper-logs/curve$r.npz; done
echo "   then:  python3 view.py paper-logs/curve0.6.npz paper-logs/curve0.8.npz -o curvature.png"; fi

echo; echo "$N runs, logs in paper-logs/"
echo "then:  python3 compare_paper.py paper-logs"
