# paper_repro.ps1 -- reproduce Garcia et al. (2019) JoCS 30:143-156, sec 3.3
# (field-scale internal wave beams over a Gaussian ridge), on the FULLY
# CURVILINEAR grid.
#
#   $env:MOLE_SRC = "C:\...\mole\src\matlab_octave"
#   .\paper_repro.ps1            # preflight, then all sections
#   .\paper_repro.ps1 -Only 2
#   python compare_paper.py paper-logs
#
# ---------------------------------------------------------------------------
# WHAT MATCHES THE PAPER
#   Lx = 3000 m, D0 = 1000 m, nz = 101      (paper: 128 x 6 x 101)
#   ridge D = D0 - ab exp(-x^2/2Lb^2), ab = 20 m, Lb = Lx/100 = 30 m   (Eq 26)
#   N = 0.007 s^-1 constant                                            (Eq 27)
#   ubc = u0 sin(wt), u0 = 0.01 m/s                                    (Eq 28)
#   sponge SL = -(u - ubc)/taus * exp(-4r/Lsl)                         (Eq 29)
#   20 tidal periods; beam angle from the max of the rms velocity over the
#   LAST TEN periods, fit over |x| in [200,500] m
#   omega/N = 0.2, 0.4, 0.6, 0.8
#
# WHAT DOES NOT, AND WHY
#   GRID.  The paper says "all simulations are conducted in sigma
#   coordinates".  This grid is fully curvilinear -- the xi-lines are not
#   vertical (73.8 m of top-to-bottom wander at Lx = 3 km) and the side walls
#   bulge by 0.15 D0.  That is deliberate and is the whole point of the
#   exercise, so the comparison is "same physics, harder grid", not a
#   like-for-like reproduction.
#
#   SPONGE.  The paper's Lsl = Lx/10, taus = 100 s does NOT work on this grid
#   at Lx = 3 km.  Measured, 20 periods, win [200,500]:
#       lsl 0.10 taus 100   w/N=0.6  rms  35.5 m  max/mean 3.3   err  +7.56
#       lsl 0.10 taus 100   w/N=0.8  rms 191.4 m  max/mean 2.6   err -24.02
#       lsl 0.20 taus  50   w/N=0.6  rms   3.6 m  max/mean 6.4   err  +0.38
#       lsl 0.20 taus  50   w/N=0.8  rms  58.0 m  max/mean 4.2   err  +8.32
#       lsl 0.30 taus  30   w/N=0.6  rms   4.3 m  max/mean 15.5  err  +0.78
#       lsl 0.30 taus  30   w/N=0.8  rms   6.3 m  max/mean 12.3  err  +0.46
#   Section 1 re-measures this so the deviation is evidenced, not asserted.
#   Whether the paper's weaker sponge suffices on a sigma grid and not here,
#   or whether their 3D domain with a lateral dimension bleeds energy
#   differently, is NOT established -- do not claim it does.
#
#   TIMESTEP.  The paper uses dt = 0.01 s.  Here dt is set per frequency so
#   that dt/taus <= 0.4, since the sponge relaxation is explicit and dt/taus
#   > 1 blows up.  The beam angle is dt-converged well before this.
# ---------------------------------------------------------------------------

param([string]$Only = "all")

$Pat = 'r_x .*max|sponge:|wander|consistency:|rms residual|max/mean|BEAM ANGLE|DIVERGED|static projection|VIOLATED|UNSTABLE|Raise --taus'
$LogDir = Join-Path (Get-Location) "paper-logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$script:N = 0
$script:SEC = "0"

# omega/N -> T (s), and spp giving dt/taus <= 0.4 at taus = 30 s
function Tper([double]$r) { return 2*[Math]::PI/($r*0.007) }
function Spp([double]$r, [double]$taus) {
    $need = [Math]::Ceiling((Tper $r)/(0.4*$taus))
    return [Math]::Max(120, $need)
}

function Sec([string]$n, [string]$t) {
    if ($Only -eq "all" -or $Only -eq $n) {
        Write-Host ""; Write-Host "=== $n. $t ===" -ForegroundColor Cyan
        # Log names carry the section number and the counter restarts per
        # section.  Without this, running -Only 7 after a full sweep rewrote
        # p001.. and destroyed the earlier sections' logs.
        $script:SEC = $n; $script:N = 0
        return $true
    }
    return $false
}

function Run([string[]]$A) {
    # A blank element here means an upstream expression failed and returned
    # nothing; python then sees "--spp --lsl" and dies on argparse.  Catch it
    # at the call site where the cause is visible.
    for ($i = 0; $i -lt $A.Count; $i++) {
        if ([string]::IsNullOrWhiteSpace($A[$i])) {
            $prev = if ($i -gt 0) { $A[$i-1] } else { "<first>" }
            Write-Host "   !! empty argument after '$prev' -- a computed value came back blank; not running" -ForegroundColor Red
            return
        }
    }
    $cmd = @("iwbcurv.py", "--mole", $env:MOLE_SRC, "--beat", "9999") + $A
    Write-Host ("   python " + ($cmd -join " ")) -ForegroundColor DarkGray
    $script:N++
    $log = Join-Path $LogDir ("s{0}_{1:d3}.log" -f $script:SEC, $script:N)
    & python @cmd 2>&1 | Out-File -FilePath $log -Encoding utf8
    $code = $LASTEXITCODE
    $hits = @(Select-String -Path $log -Pattern $Pat | ForEach-Object { $_.Line.TrimEnd() })
    if ($hits.Count -gt 0) { $hits | ForEach-Object { Write-Host "   $_" } }
    if ($hits.Count -eq 0 -or $code -ne 0) {
        Write-Host "   !! no summary lines (exit=$code) -- tail of $log :" -ForegroundColor Red
        Get-Content $log -Tail 25 | ForEach-Object { Write-Host "   | $_" -ForegroundColor DarkYellow }
    }
}

# ------------------------------------------------------------------ 0
Write-Host "=== 0. preflight ===" -ForegroundColor Cyan
$fail = $false
foreach ($f in "iwbcurv.py", "view.py", "compare_paper.py") {
    if (Test-Path $f) { Write-Host "   OK   $f" }
    else { Write-Host "   MISSING  $f" -ForegroundColor Red; $fail = $true }
}
if (-not (Test-Path "grids\iwbridge\left.m")) {
    Write-Host "   MISSING  grids\iwbridge\  (run from your iwb-curv dir)" -ForegroundColor Red; $fail = $true
} else {
    $b = Select-String -Path "grids\iwbridge\left.m" -Pattern 'bulge\s*=\s*([0-9.]+)\*D0'
    if ($b) {
        $bv = [double]$b.Matches[0].Groups[1].Value
        if ($bv -gt 0) { Write-Host "   OK   grids\iwbridge\  (bulge = $bv D0, fully curvilinear)" }
        else { Write-Host "   grids\iwbridge\left.m has bulge = 0 -- NOT the curvilinear grid" -ForegroundColor Red; $fail = $true }
    } else { Write-Host "   could not read bulge from left.m" -ForegroundColor Yellow }
}
if (-not $env:MOLE_SRC) { Write-Host "   MOLE_SRC not set" -ForegroundColor Red; $fail = $true }
elseif (-not (Test-Path (Join-Path $env:MOLE_SRC "gridGen.m"))) {
    Write-Host "   MOLE_SRC has no gridGen.m: $env:MOLE_SRC" -ForegroundColor Red; $fail = $true
} else { Write-Host "   OK   MOLE_SRC = $env:MOLE_SRC" }
Write-Host "   -- grid smoke test (must say FULLY CURVILINEAR) --"
Run @("--Lx","3000","--nx","128","--nz","101","--ratio","0.6","--gridonly")
if ($fail) { Write-Host "`npreflight failed" -ForegroundColor Red; exit 1 }
if ($Only -eq "0") { exit 0 }

# ------------------------------------------------------------------ 1
if (Sec "1" "sponge: the paper's Lsl=Lx/10, taus=100 vs alternatives (~25 min)") {
  $cfg = @( @("0.10","100"), @("0.20","50"), @("0.30","30") )
  foreach ($c in $cfg) {
    $lsl = $c[0]; $taus = [double]$c[1]
    foreach ($r in "0.2","0.4","0.6","0.8") {
      $spp = Spp ([double]$r) $taus
      Write-Host "-- lsl=$lsl taus=$taus  w/N=$r  spp=$spp"
      Run @("--Lx","3000","--nx","128","--nz","101","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl",$lsl,"--taus","$taus","--win","200","500")
    }
  }
}

# ------------------------------------------------------------------ 2
if (Sec "2" "THE COMPARISON: paper config, paper grid size, calibrated sponge (~10 min)") {
  foreach ($r in "0.2","0.4","0.6","0.8") {
    $spp = Spp ([double]$r) 30.0
    Write-Host "-- w/N=$r  spp=$spp"
    Run @("--Lx","3000","--nx","128","--nz","101","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.30","--taus","30","--win","200","500",
          "--out","paper-logs\pf$r.npz")
  }
}

# ------------------------------------------------------------------ 3
if (Sec "3" "convergence: refine the paper config (~60 min)") {
  $cfg = @( @("192","151"), @("256","201") )
  foreach ($c in $cfg) {
    $nx = $c[0]; $nz = $c[1]
    foreach ($r in "0.2","0.4","0.6","0.8") {
      $spp = Spp ([double]$r) 30.0
      Write-Host "-- w/N=$r  nx=$nx nz=$nz  spp=$spp"
      Run @("--Lx","3000","--nx",$nx,"--nz",$nz,"--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.30","--taus","30","--win","200","500",
            "--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 4
if (Sec "4" "Fig 8 reproduction") {
  Write-Host "   paper Fig 8 is TOTAL u/u0 at t=19.5T (ubc=0), colour limit +/-0.15"
  & python view.py paper-logs\pf0.2.npz paper-logs\pf0.4.npz `
                   paper-logs\pf0.6.npz paper-logs\pf0.8.npz `
                   --total --clim 0.15 -o fig8_repro.png
  Write-Host "   compare against Fig 8 on p.151 of the PDF"
}

# ------------------------------------------------------------------ 5
if (Sec "5" "comparison table vs Fig 9") { & python compare_paper.py paper-logs }

# ------------------------------------------------------------------ 6
if (Sec "6" "WHY refinement diverges: window and sponge sensitivity at 256x201 (~35 min)") {
  # Sections 2-3 show the bias crossing zero and GROWING as the grid refines,
  # for w/N = 0.4 and 0.8.  That is the opposite of the Lx=6 km / auto-window
  # behaviour (+1.99 -> +1.06 -> +0.37) and has to be explained before the
  # agreement at 128x101 means anything.
  #
  # Two candidates, and this section separates them:
  #   (a) THE WINDOW.  The paper fits |x| in [200,500] m at every frequency.
  #       At w/N=0.2 the bounce is 4900 m, so that window is the first 10% of
  #       a bounce -- near field, where HANDOFF sec 4 already recorded the
  #       ridge angle not yet settled onto the characteristic.  If this is it,
  #       moving the window outward at FIXED grid pulls the angle back toward
  #       theory, and the auto window agrees with the Lx=6 km result.
  #   (b) THE SPONGE.  Lsl = 0.30 Lx = 900 m in a +/-1500 m domain reaches
  #       within 600 m of the fit window.  If this is it, the answer at fixed
  #       grid depends on lsl/taus.  It should not.
  # If neither moves it, both are wrong and the beam itself is genuinely
  # steepening with resolution -- which would be the real finding.
  foreach ($r in "0.4","0.8") {
    $spp = Spp ([double]$r) 30.0
    foreach ($w in @("200","500"), @("300","700"), @("400","900")) {
      Write-Host "-- (a) window w/N=$r  [$($w[0]),$($w[1])]"
      Run @("--Lx","3000","--nx","256","--nz","201","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.30","--taus","30","--win",$w[0],$w[1],"--forcegrid")
    }
    Write-Host "-- (a) window w/N=$r  auto"
    Run @("--Lx","3000","--nx","256","--nz","201","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.30","--taus","30","--forcegrid")
    foreach ($c in @("0.20","50"), @("0.40","30")) {
      Write-Host "-- (b) sponge w/N=$r  lsl=$($c[0]) taus=$($c[1])"
      Run @("--Lx","3000","--nx","256","--nz","201","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl",$c[0],"--taus",$c[1],"--win","200","500","--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 7
if (Sec "7" "is the divergence a 3 km artifact? same test at Lx = 6 km (~90 min)") {
  # Section 6 found the answer at w/N=0.8, 256x201, win [200,500] moves 1.14
  # deg between sponge 0.30/30 and 0.40/30.  A parameter that should not
  # affect the answer moves it by more than the bias being measured, so the
  # 3 km refinement series is not measuring discretization error.
  #
  # The cause is geometric.  At Lx = 3 km this grid needs Lsl = 0.30-0.40 Lx
  # (900-1200 m) to suppress reflections, in a +/-1500 m box -- so the sponge
  # edge sits within 600 m of a [200,500] window.  Suppressing reflections and
  # keeping the sponge off the window are in conflict at 3 km.  At 6 km they
  # are not: lsl = 0.20 puts the sponge 2500 m from the window edge
  # (exp(-4*2500/1200) = 2e-4).
  #
  # Same paper window, same ridge, dx and dz held to the 3 km values by
  # doubling nx.  If the bias now converges toward zero and the sponge spread
  # drops below ~0.2 deg, the divergence was the 3 km domain and the paper's
  # own domain is the limitation.  If it still diverges with a clean sponge,
  # the divergence is real and belongs to the solver.
  $cfg = @( @("256","101"), @("384","151"), @("512","201") )
  foreach ($c in $cfg) {
    $nx = $c[0]; $nz = $c[1]
    foreach ($r in "0.2","0.4","0.6","0.8") {
      $taus = (Tper ([double]$r))/30.0
      $spp = Spp ([double]$r) $taus
      Write-Host "-- w/N=$r  nx=$nx nz=$nz  lsl=0.20 taus=$([Math]::Round($taus,1))  spp=$spp"
      Run @("--Lx","6000","--nx",$nx,"--nz",$nz,"--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win","200","500","--forcegrid")
    }
  }
  Write-Host "-- sponge independence check at the finest grid"
  foreach ($r in "0.4","0.8") {
    foreach ($l in "0.15","0.25") {
      $taus = (Tper ([double]$r))/30.0
      $spp = Spp ([double]$r) $taus
      Write-Host "-- w/N=$r  lsl=$l"
      Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl",$l,"--taus","$taus","--win","200","500","--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 8
if (Sec "8" "window sweep at 6 km with a CLEAN sponge, windows as bounce fractions (~30 min)") {
  # Section 6 did a window sweep, but at Lx = 3 km, where the sponge spread is
  # 1.14 deg at w/N=0.8 -- so it could not separate window from sponge.  And it
  # used fixed metre windows, which sit at wildly different fractions of the
  # bounce for each frequency (bounce = D0/tan(phi) is 4900 m at w/N=0.2 and
  # 750 m at 0.8), so "the same window" is not the same place along the beam.
  #
  # Here: Lx = 6 km, 512x201, lsl = 0.20 (sponge spread measured at 0.13-0.26
  # deg, so clean), and the windows are FRACTIONS OF THE BOUNCE so 0.4, 0.6 and
  # 0.8 are measured at the same place along their own beams.
  #
  # If the converged bias tracks bounce fraction and collapses across
  # frequencies when plotted that way, the residual is the measurement
  # convention -- the beam is not straight at the theory angle, and [200,500]
  # samples a different part of it at each frequency.  If the bias stays
  # frequency-dependent at matched bounce fraction, that hypothesis is dead
  # too and the residual belongs to the solver.
  #
  # w/N=0.2 is omitted: its bounce is 4900 m, so even 0.30 of a bounce reaches
  # past 1470 m and the outer windows would run into the sponge at 6 km.  It
  # needs Lx = 12 km, which is a separate (expensive) run.
  foreach ($r in "0.4","0.6","0.8") {
    $phi = [Math]::Asin([double]$r)
    $b   = 1000.0/[Math]::Tan($phi)
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    foreach ($f in @(0.10,0.30), @(0.20,0.45), @(0.30,0.60), @(0.45,0.75)) {
      $w0 = [Math]::Round($b*$f[0]); $w1 = [Math]::Round($b*$f[1])
      if ($w1 -gt 1800) {
        Write-Host "-- w/N=$r  bounce frac $($f[0])-$($f[1]) -> [$w0,$w1] m  SKIPPED (would reach the sponge)" -ForegroundColor Yellow
        continue
      }
      Write-Host "-- w/N=$r  bounce $([Math]::Round($b)) m  frac $($f[0])-$($f[1])  -> win [$w0,$w1] m"
      Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1","--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 9
if (Sec "9" "is the beam curvature grid-converged? (~25 min)") {
  # Section 8: at Lx=6 km, 512x201, clean sponge, the measured slope varies
  # across the first bounce by 0.26 deg at w/N=0.4 but 2.12 deg at 0.6 and
  # 3.26 deg at 0.8.  The beam is not straight at the higher frequencies, so
  # there is no single angle to compare with theory there.
  #
  # If that curvature is PHYSICAL, the window-to-window spread is the same at
  # 384x151 as at 512x201.  If it is NUMERICAL, the spread shrinks under
  # refinement.  512x201 is already done in section 8; this fills in 384x151
  # at the identical bounce fractions, and saves a field at the widest window
  # so the tracked maxima can be looked at directly.
  foreach ($r in "0.6","0.8") {
    $phi = [Math]::Asin([double]$r)
    $b   = 1000.0/[Math]::Tan($phi)
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    foreach ($f in @(0.10,0.30), @(0.20,0.45), @(0.30,0.60), @(0.45,0.75)) {
      $w0 = [Math]::Round($b*$f[0]); $w1 = [Math]::Round($b*$f[1])
      Write-Host "-- w/N=$r  384x151  frac $($f[0])-$($f[1])  -> win [$w0,$w1] m"
      Run @("--Lx","6000","--nx","384","--nz","151","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1","--forcegrid")
    }
    # a field at the widest usable window, for looking at the maxima directly
    $w0 = [Math]::Round($b*0.10); $w1 = [Math]::Round($b*0.75)
    Write-Host "-- w/N=$r  512x201  full first bounce [$w0,$w1] -> npz"
    Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1",
          "--forcegrid","--out","paper-logs\curve$r.npz")
  }
  Write-Host "   then:  python view.py paper-logs\curve0.6.npz paper-logs\curve0.8.npz -o curvature.png"
  Write-Host "   the white tracked maxima should bend away from the orange straight-line fit"
}

# ------------------------------------------------------------------ 10
if (Sec "10" "THE COMPARISON, done properly: fit the whole first bounce (~45 min)") {
  # Section 9 settled it.  The beam is STRAIGHT -- a line fits the tracked
  # maxima to 2.9 m over 867 m at w/N=0.6 -- but it carries a small oscillation
  # about that line.  A short fit window measures the local slope of the
  # wiggle; a long one averages it out.
  #
  # TWO effects, not one.  Converting each spread to an implied deviation
  # amplitude (spread in rad x mean window span):
  #     w/N=0.4  spans 458-688 m  spread 0.26 deg  ->  ~3 m
  #     w/N=0.6  spans 267-400 m  spread 2.12 deg  ->  ~13 m
  #     w/N=0.8  spans 150-225 m  spread 3.26 deg  ->  ~11 m
  # So 0.4's windows are longer AND its oscillation is ~4x smaller.  The
  # lever arm alone does not account for it; do not claim it does.
  #
  # So the well-posed measurement is a fit over the whole first bounce, and
  # that is what should be compared with Fig 9.  Already measured at 512x201:
  #     w/N=0.6  [133,1000]  180 pts  rms 2.9 m  ->  +0.09 deg
  #     w/N=0.8  [75,562]    104 pts  rms 4.3 m  ->  -0.69 deg
  # This section fills in 0.4 and 0.2 the same way.  0.2 needs Lx = 12 km
  # because its bounce is 4900 m; nx is doubled to hold dx, and the window
  # stops at 0.65 of a bounce to stay clear of the sponge.
  foreach ($c in @("0.4","6000","512"), @("0.2","12000","1024")) {
    $r = $c[0]; $lx = [double]$c[1]; $nx = $c[2]
    $b = 1000.0/[Math]::Tan([Math]::Asin([double]$r))
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    $hi = if ($r -eq "0.2") { 0.65 } else { 0.75 }
    $w0 = [Math]::Round($b*0.10); $w1 = [Math]::Round($b*$hi)
    $clear = $lx/2 - 0.20*$lx
    Write-Host "-- w/N=$r  Lx=$lx nx=$nx  bounce $([Math]::Round($b)) m  win [$w0,$w1] m  (sponge-free to $clear m)"
    Run @("--Lx","$lx","--nx",$nx,"--nz","201","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1",
          "--forcegrid","--out","paper-logs\bounce$r.npz")
  }
  # and re-save 0.6 / 0.8 so all four npz exist for one figure
  foreach ($c in @("0.6","133","1000"), @("0.8","75","562")) {
    $r = $c[0]
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    Write-Host "-- w/N=$r  re-save full-bounce field"
    Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.20","--taus","$taus","--win",$c[1],$c[2],
          "--forcegrid","--out","paper-logs\bounce$r.npz")
  }
  Write-Host "   then:  python view.py paper-logs\bounce0.2.npz paper-logs\bounce0.4.npz paper-logs\bounce0.6.npz paper-logs\bounce0.8.npz --total --clim 0.15 -o fig8_final.png"
}

# ------------------------------------------------------------------ 11
if (Sec "11" "near-field: push the INNER window edge out, outer edge fixed (~40 min)") {
  # The plot plus the two-grid check say the residual about the fitted line is
  # ONE ARCH, not an oscillation, and its amplitude scales as dx^1.4 -- so it
  # shrinks with refinement.  The innermost window (0.10-0.30 of a bounce) is
  # the exception: +1.55 at 384x151 and +1.63 at 512x201, essentially
  # unchanged.  That looks like the ridge near field, which is physical and
  # will not refine away.
  #
  # Test: hold the OUTER edge at 0.75 of a bounce and walk the INNER edge out.
  # If the bias falls as more near field is excluded and then flattens, the
  # plateau IS the far-field beam angle and the near field is the whole story.
  # If it keeps drifting with no plateau, it is not near field.
  #
  # Note the confound: moving the inner edge out also SHORTENS the window, and
  # a shorter window has a worse lever arm.  So a drift that continues past the
  # point where the span drops below ~300 m may be lever arm, not physics.
  # The span is printed for each run so this can be judged, not assumed.
  foreach ($r in "0.4","0.6","0.8") {
    $b = 1000.0/[Math]::Tan([Math]::Asin([double]$r))
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    $w1 = [Math]::Round($b*0.75)
    foreach ($f in 0.10, 0.20, 0.30, 0.40) {
      $w0 = [Math]::Round($b*$f)
      Write-Host "-- w/N=$r  inner $f of bounce -> win [$w0,$w1] m  (span $($w1-$w0) m)"
      Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1","--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 12
if (Sec "12" "fill in the coarse grids so every window has 3 points (~20 min)") {
  # Sections 8 and 9 gave the bounce-fraction windows at 512x201 (all three
  # frequencies) and at 384x151 (0.6 and 0.8 only).  With two grids per window
  # no extrapolation is possible, and two-point extrapolation is what led to
  # the claim that the full-bounce fit had converged -- it had not, it was
  # crossing zero on the way down.
  #
  # Measured 384x151 -> 512x201, the angle falls by 0.45-1.22 deg at every
  # window EXCEPT the innermost two, which are stationary (+0.08, -0.29).  So
  # the near field is resolved and the far field is not, which is backwards
  # from what near-field curvature would predict.
  #
  # This adds 256x101 for all three frequencies and 384x151 for 0.4, so every
  # window has three grids and compare_paper.py can extrapolate each one.
  foreach ($c in @("0.4","256","101"), @("0.6","256","101"), @("0.8","256","101"),
                 @("0.4","384","151")) {
    $r = $c[0]; $nx = $c[1]; $nz = $c[2]
    $b = 1000.0/[Math]::Tan([Math]::Asin([double]$r))
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    foreach ($f in @(0.10,0.30), @(0.20,0.45), @(0.30,0.60), @(0.45,0.75)) {
      $w0 = [Math]::Round($b*$f[0]); $w1 = [Math]::Round($b*$f[1])
      Write-Host "-- w/N=$r  ${nx}x${nz}  frac $($f[0])-$($f[1])  -> win [$w0,$w1] m"
      Run @("--Lx","6000","--nx",$nx,"--nz",$nz,"--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win","$w0","$w1","--forcegrid")
    }
  }
  # and the full-bounce window at all three grids for 0.6 and 0.8
  foreach ($c in @("0.6","133","1000"), @("0.8","75","562")) {
    $r = $c[0]
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    foreach ($g in @("256","101"), @("384","151")) {
      Write-Host "-- w/N=$r  $($g[0])x$($g[1])  full bounce [$($c[1]),$($c[2])]"
      Run @("--Lx","6000","--nx",$g[0],"--nz",$g[1],"--ratio",$r,"--nper","20",
            "--spp","$spp","--lsl","0.20","--taus","$taus","--win",$c[1],$c[2],"--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 13
if (Sec "13" "is the long-window agreement a limit or a crossing? one finer grid (~45 min)") {
  # Section 11: with the OUTER edge fixed at 0.75 of a bounce, pushing the
  # inner edge out changes the answer by 0.07 deg at w/N=0.4 and 0.28 at 0.6.
  # Near field is not the story; the outer edge is.  With windows reaching
  # 0.75 of a bounce at 512x201 the four frequencies give
  #     0.2 +0.43   0.4 +0.17   0.6 +0.09   0.8 -0.69
  # mean |bias| 0.35 deg, against 1.51 deg for GCCOM's Fig 9 points.
  #
  # BUT the grid series for those same windows is still drifting DOWN:
  #     w/N=0.6 [133,1000]  +1.86, +0.70, +0.09   extrap -0.59
  #     w/N=0.8 [75,562]    +1.64, +0.01, -0.69   extrap -1.22
  # so 512x201 may be a zero crossing, not a limit.  640x251 is a 1.25x
  # refinement on each axis.  If the next difference is much smaller than the
  # last (-0.61 at 0.6, -0.70 at 0.8), the drift is stopping and the agreement
  # is real.  If it continues at the same size, the values keep falling and
  # the agreement at 512x201 is a coincidence of resolution.
  foreach ($c in @("0.4","229","1718"), @("0.6","133","1000"), @("0.8","75","562")) {
    $r = $c[0]
    $taus = (Tper ([double]$r))/30.0
    $spp = Spp ([double]$r) $taus
    Write-Host "-- w/N=$r  640x251  win [$($c[1]),$($c[2])]"
    Run @("--Lx","6000","--nx","640","--nz","251","--ratio",$r,"--nper","20",
          "--spp","$spp","--lsl","0.20","--taus","$taus","--win",$c[1],$c[2],"--forcegrid")
  }
}

Write-Host ""
Write-Host "$script:N runs, logs in $LogDir" -ForegroundColor Cyan
Write-Host "then:  python compare_paper.py paper-logs" -ForegroundColor Cyan
