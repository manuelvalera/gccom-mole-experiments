# seiche_study.ps1 -- find why the seiche frequency error GROWS under refinement.
#
#   $env:MOLE_SRC = "C:\...\mole\src\octave"     # note: octave, not matlab_octave
#   .\seiche_study.ps1              # all sections
#   .\seiche_study.ps1 -Only 3
#   python seiche_report.py seiche-logs
#
# ---------------------------------------------------------------------------
# WHAT IS BEING CHASED
#
# A free standing internal wave in a flat rectangular box has an exact
# frequency, omega = N k / sqrt(k^2 + p^2), k = I pi / Lx, p = J pi / D0.
# No forcing, no sponge, no topography, no fit window, no tracker.  Measured
# at mode (4,2), spp=400:
#
#     grid        curvilinear      near-Cartesian
#     128x51        +2.91e-3          +5.04e-3
#     192x76        +3.95e-3          +6.02e-3
#     256x101       +4.32e-3          +6.36e-3
#
# The error GROWS as the grid refines, on both grids.  That is not slow
# convergence, it is convergence to the wrong answer.
#
# The buoyancy interpolation is NOT the explanation.  Its symbol is exactly
# cos^2(k dz/2), so N_eff = N cos(k dz/2) -- second order, negative, and
# vanishing.  Refining at FIXED mode makes k*dz smaller, so that error must
# shrink.  Something else grows.
#
# PRIME SUSPECT, stated before running: the lateral boundary.  With no sponge,
# the side walls are held only by the Robin rows from robinBC2D, alpha*phi +
# dphi/dn, with alpha = 1e-4 fixed.  The discrete normal derivative scales
# like 1/dz, so as the grid refines the derivative term grows relative to the
# fixed alpha and the EFFECTIVE boundary condition drifts from Robin toward
# pure Neumann.  A boundary condition that changes with resolution is exactly
# the kind of thing that converges to the wrong operator.
#
# If that is right, section 3 shows the error depending on alpha, and the
# refinement trend changing sign or flattening as alpha is scaled with the
# grid.  If alpha does nothing, the suspect is wrong and sections 1, 2 and 5
# narrow it instead.
# ---------------------------------------------------------------------------

param([string]$Only = "all")

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'SEICHE|seiche mode|initial field projected|grid .*dx=|crossings|UNSTABLE|DIVERGED|\*\*\*'
$LogDir = Join-Path (Get-Location) "seiche-logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$script:SEC = "0"; $script:N = 0

function Sec([string]$n, [string]$t) {
    if ($Only -eq "all" -or $Only -eq $n) {
        Write-Host ""; Write-Host "=== $n. $t ===" -ForegroundColor Cyan
        $script:SEC = $n; $script:N = 0
        return $true
    }
    return $false
}

function Run([string[]]$A) {
    foreach ($x in $A) {
        if ([string]::IsNullOrWhiteSpace($x)) {
            Write-Host "   !! blank argument, not running" -ForegroundColor Red; return
        }
    }
    $cmd = @("iwbcurv.py", "--mole", $env:MOLE_SRC, "--beat", "9999",
             "--ab", "0", "--bulge", "0") + $A
    $script:N++
    $log = Join-Path $LogDir ("s{0}_{1:d3}.log" -f $script:SEC, $script:N)
    Write-Host ("   python " + ($cmd -join " ")) -ForegroundColor DarkGray
    # the parser reads parameters back out of this echoed line
    "CMD " + ($cmd -join " ") | Out-File -FilePath $log -Encoding utf8
    & python @cmd 2>&1 | Out-File -FilePath $log -Encoding utf8 -Append
    $code = $LASTEXITCODE
    $hits = @(Select-String -Path $log -Pattern $Pat | ForEach-Object { $_.Line.TrimEnd() })
    if ($hits.Count -gt 0) { $hits | ForEach-Object { Write-Host "   $_" } }
    if ($hits.Count -eq 0 -or $code -ne 0) {
        Write-Host "   !! no result (exit=$code), tail of $log :" -ForegroundColor Red
        Get-Content $log -Tail 20 | ForEach-Object { Write-Host "   | $_" -ForegroundColor DarkYellow }
    }
}

# ------------------------------------------------------------------ 0
Write-Host "=== 0. preflight ===" -ForegroundColor Cyan
foreach ($f in "iwbcurv.py", "grids\iwbridge\left.m") {
    if (Test-Path $f) { Write-Host "   OK   $f" } else { Write-Host "   MISSING $f" -ForegroundColor Red; exit 1 }
}
if (Test-Path (Join-Path $env:MOLE_SRC "gridGen.m")) { Write-Host "   OK   MOLE_SRC" }
else { Write-Host "   MOLE_SRC has no gridGen.m" -ForegroundColor Red; exit 1 }
Write-Host "   -- one short seiche, must report a relative error --"
Run @("--Lx","6000","--nx","128","--nz","51","--seiche","4","2","--nper","6","--spp","400")
if ($Only -eq "0") { exit 0 }

# ------------------------------------------------------------------ 1
if (Sec "1" "mode scan at FIXED grid -- does the error track k*dz? (~3 min)") {
  # If the error is an interpolation/symbol effect it should grow with the
  # vertical mode number J (larger p = J pi / D0, so larger p*dz) and be
  # roughly flat in I.  If it is flat in both, it is not a symbol problem.
  foreach ($m in @("2","1"), @("4","1"), @("8","1"), @("2","2"), @("4","2"),
                  @("2","4"), @("4","4"), @("2","8")) {
    Write-Host "-- mode ($($m[0]),$($m[1]))"
    Run @("--Lx","6000","--nx","128","--nz","51","--seiche",$m[0],$m[1],
          "--nper","6","--spp","400")
  }
}

# ------------------------------------------------------------------ 2
if (Sec "2" "refinement at fixed mode, two modes (~6 min)") {
  # The headline anomaly, repeated for a second mode so it is not specific
  # to (4,2).
  foreach ($m in @("4","2"), @("2","1")) {
    foreach ($g in @("128","51"), @("192","76"), @("256","101"), @("384","151")) {
      Write-Host "-- mode ($($m[0]),$($m[1]))  $($g[0])x$($g[1])"
      Run @("--Lx","6000","--nx",$g[0],"--nz",$g[1],"--seiche",$m[0],$m[1],
            "--nper","6","--spp","400")
    }
  }
}

# ------------------------------------------------------------------ 3
if (Sec "3" "THE SUSPECT: alpha, fixed and grid-scaled (~8 min)") {
  # 3a: does the answer depend on alpha at all, at one grid?
  foreach ($a in "1e-6","1e-5","1e-4","1e-3","1e-2","1e-1") {
    Write-Host "-- alpha=$a  256x101"
    Run @("--Lx","6000","--nx","256","--nz","101","--seiche","4","2",
          "--nper","6","--spp","400","--alpha",$a)
  }
  # 3b: refinement at FIXED alpha (the current behaviour) vs alpha SCALED
  # with 1/dz, which keeps alpha and the discrete normal derivative in the
  # same ratio as the grid refines.  dz = 1000/(nz-1), so alpha ~ 1/dz means
  # alpha = 1e-4 * (nz-1)/50 normalised to the 128x51 case.
  $cfg = @( @("128","51","1.0e-4"), @("192","76","1.5e-4"),
            @("256","101","2.0e-4"), @("384","151","3.0e-4") )
  foreach ($c in $cfg) {
    Write-Host "-- alpha scaled with 1/dz: $($c[0])x$($c[1])  alpha=$($c[2])"
    Run @("--Lx","6000","--nx",$c[0],"--nz",$c[1],"--seiche","4","2",
          "--nper","6","--spp","400","--alpha",$c[2])
  }
}

# ------------------------------------------------------------------ 4
if (Sec "4" "operator order on the seiche (~3 min)") {
  # The beam-angle version of this test showed nothing, but the beam was
  # confounded by the ridge, the sponge and the window.  This one is not.
  foreach ($k in "2","4","6") {
    foreach ($g in @("128","51"), @("256","101")) {
      Write-Host "-- k=$k  $($g[0])x$($g[1])"
      Run @("--Lx","6000","--nx",$g[0],"--nz",$g[1],"--seiche","4","2",
            "--nper","6","--spp","400","--order",$k)
    }
  }
}

# ------------------------------------------------------------------ 5
if (Sec "5" "separate dx from dz (~5 min)") {
  # Refine one axis at a time.  If the error grows only when dz shrinks, the
  # vertical discretization or the top/bottom rows own it; only when dx
  # shrinks, the lateral rows do.
  Write-Host "-- dz fixed at nz=51, refine nx"
  foreach ($nx in "96","128","192","256") {
    Run @("--Lx","6000","--nx",$nx,"--nz","51","--seiche","4","2","--nper","6","--spp","400")
  }
  Write-Host "-- dx fixed at nx=128, refine nz"
  foreach ($nz in "41","51","76","101") {
    Run @("--Lx","6000","--nx","128","--nz",$nz,"--seiche","4","2","--nper","6","--spp","400")
  }
}

# ------------------------------------------------------------------ 6
if (Sec "6" "control: is the bed constraint involved? (~2 min)") {
  # With ab=0 the bed is flat, so the constraint reduces to w=0 there and
  # should be exactly satisfiable.  If bedmode changes the frequency error,
  # the constraint rows are contributing.
  foreach ($bm in "constraint","post","none") {
    Write-Host "-- bedmode=$bm  256x101"
    Run @("--Lx","6000","--nx","256","--nz","101","--seiche","4","2",
          "--nper","6","--spp","400","--bedmode",$bm)
  }
}

Write-Host ""
Write-Host "$script:N runs in the last section; logs in $LogDir" -ForegroundColor Cyan
Write-Host "now:  python seiche_report.py seiche-logs" -ForegroundColor Cyan
