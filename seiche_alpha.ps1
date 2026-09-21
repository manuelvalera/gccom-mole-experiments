# seiche_alpha.ps1 -- alpha is the driver. Find its floor, then find out what
# it did to the beam angles.
#
#   $env:MOLE_SRC = "...\mole\src\octave"
#   .\seiche_alpha.ps1
#   .\seiche_alpha.ps1 -Only C
#   python seiche_report.py alpha-logs
#
# ---------------------------------------------------------------------------
# WHERE THIS CAME FROM
#
# The seiche sweep found the frequency error is dominated by alpha, the
# coefficient in the lateral Robin rows, not by the discretization:
#
#     alpha    rel err      half-period scatter
#     1e-6     -4.08e-04    +/- 0.007 s
#     1e-5     +2.42e-05    +/- 0.089 s
#     1e-4     +4.32e-03    +/- 0.99  s      <- the default
#     1e-3     +4.36e-02    +/- 13.6  s
#
# alpha exists only to regularise the Poisson solve: a solid wall is pure
# Neumann (dphi/dn = 0), which leaves the system singular up to a constant,
# and alpha*phi is the smallest thing that removes the null space. So the
# physically correct wall is alpha -> 0, and the right value is the smallest
# one the factorisation tolerates.
#
# DO NOT read 1e-5 as the answer. The series is monotonic through zero between
# 1e-6 and 1e-5, so 1e-5 is a zero CROSSING -- a value that happens to cancel
# the residual error at this grid and this mode. Tuning alpha to make one
# number vanish is fitting, not fixing. The quantity that means something is
# the alpha -> 0 limit, which looks like about -4e-4 and is still unmeasured.
#
# Section A finds that limit and where conditioning gives out.
# Section B re-runs the refinement series there, for the TRUE convergence order.
# Section C is the one that matters for the paper: every beam angle quoted so
# far was computed at alpha = 1e-4.
# ---------------------------------------------------------------------------

param([string]$Only = "all")

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'SEICHE|BEAM ANGLE|seiche mode|initial field projected|crossings|rms residual|max/mean|consistency|UNSTABLE|DIVERGED|\*\*\*'
$LogDir = Join-Path (Get-Location) "alpha-logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$script:SEC = "A"; $script:N = 0

function Sec([string]$n, [string]$t) {
    if ($Only -eq "all" -or $Only -eq $n) {
        Write-Host ""; Write-Host "=== $n. $t ===" -ForegroundColor Cyan
        $script:SEC = $n; $script:N = 0; return $true
    }
    return $false
}

function Run([string[]]$A) {
    foreach ($x in $A) {
        if ([string]::IsNullOrWhiteSpace($x)) {
            Write-Host "   !! blank argument, not running" -ForegroundColor Red; return
        }
    }
    $cmd = @("iwbcurv.py", "--mole", $env:MOLE_SRC, "--beat", "9999") + $A
    $script:N++
    $log = Join-Path $LogDir ("s{0}_{1:d3}.log" -f $script:SEC, $script:N)
    Write-Host ("   python " + ($cmd -join " ")) -ForegroundColor DarkGray
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

# ------------------------------------------------------------------ A
if (Sec "A" "alpha floor: how small can it go, and what is the limit? (~6 min)") {
  # Watch TWO things in each line: the relative error, and the
  # "initial field projected" residuals. The floor is where div/bed start
  # climbing -- that is the factorisation losing the null space, and it is
  # the real constraint on how small alpha can be.
  foreach ($a in "1e-9","1e-8","1e-7","1e-6","3e-6","1e-5") {
    Write-Host "-- alpha=$a  256x101  mode (4,2)"
    Run @("--ab","0","--bulge","0","--Lx","6000","--nx","256","--nz","101",
          "--seiche","4","2","--nper","6","--spp","400","--alpha",$a)
  }
}

# ------------------------------------------------------------------ B
if (Sec "B" "TRUE convergence order at a small alpha (~8 min)") {
  # Same series as before, at alpha = 1e-6 rather than 1e-4. If the error
  # now SHRINKS under refinement, the growth was alpha contamination and the
  # discretization is fine. If it still grows at a floor-level alpha, there
  # is a second effect underneath and this is not finished.
  foreach ($g in @("128","51"), @("192","76"), @("256","101"), @("384","151")) {
    Write-Host "-- alpha=1e-6  $($g[0])x$($g[1])"
    Run @("--ab","0","--bulge","0","--Lx","6000","--nx",$g[0],"--nz",$g[1],
          "--seiche","4","2","--nper","6","--spp","400","--alpha","1e-6")
  }
  # and a second mode, so the order is not a property of (4,2) alone
  foreach ($g in @("128","51"), @("256","101")) {
    Write-Host "-- alpha=1e-6  mode (2,1)  $($g[0])x$($g[1])"
    Run @("--ab","0","--bulge","0","--Lx","6000","--nx",$g[0],"--nz",$g[1],
          "--seiche","2","1","--nper","6","--spp","400","--alpha","1e-6")
  }
}

# ------------------------------------------------------------------ C
if (Sec "C" "THE ONE THAT MATTERS: do the beam angles move? (~25 min)") {
  # Every Fig 9 number so far was computed at alpha = 1e-4, which the seiche
  # shows is the dominant error there. These runs are the forced problem WITH
  # the sponge, so alpha is not the only thing holding the walls and the
  # sensitivity may well be smaller -- but it has never been checked, and if
  # the angles move by more than ~0.1 deg the validation table has to be
  # requoted before anything is written up.
  #
  # theory: 0.6 -> 36.87 deg, window [133,1000], taus = T/30 = 49.87
  #         0.8 -> 53.13 deg, window [75,562],  taus = T/30 = 37.40
  foreach ($cfg in @("0.6","49.87","133","1000"), @("0.8","37.40","75","562")) {
    foreach ($a in "1e-4","1e-5","1e-6") {
      Write-Host "-- ratio=$($cfg[0])  alpha=$a  512x201"
      Run @("--Lx","6000","--nx","512","--nz","201","--ratio",$cfg[0],
            "--nper","20","--spp","120","--lsl","0.20","--taus",$cfg[1],
            "--win",$cfg[2],$cfg[3],"--alpha",$a,"--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ D
if (Sec "D" "control: does a small alpha survive a long run? (~6 min)") {
  # alpha regularises the solve. A value that is fine for 6 periods may drift
  # over 40, because the null-space component grows slowly. This is the check
  # that a smaller alpha is actually usable, not just accurate briefly.
  foreach ($a in "1e-4","1e-6") {
    Write-Host "-- alpha=$a  40 periods"
    Run @("--ab","0","--bulge","0","--Lx","6000","--nx","192","--nz","76",
          "--seiche","4","2","--nper","40","--spp","400","--alpha",$a)
  }
}

Write-Host ""
Write-Host "logs in $LogDir" -ForegroundColor Cyan
Write-Host "now:  python seiche_report.py alpha-logs" -ForegroundColor Cyan
Write-Host "then: python compare_paper.py alpha-logs   (for section C)" -ForegroundColor Cyan
