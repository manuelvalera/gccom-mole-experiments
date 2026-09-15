# reproduce.ps1 -- regenerate every table from this session.
#
#   $env:MOLE_SRC = "C:\path\to\mole\src\matlab_octave"
#   .\reproduce.ps1            # preflight, then everything
#   .\reproduce.ps1 -Only 4    # preflight, then section 4
#   .\reproduce.ps1 -Only 0    # preflight only
#
# Every run is written in full to .\repro-logs\run_NNN.log.  Summary lines are
# echoed; if a run produces NO summary lines or exits nonzero, the tail of its
# log is dumped instead.  (The first version of this script piped stderr into
# Select-String, so failing runs printed nothing at all and the sweep looked
# like it completed instantly.  Do not reintroduce that.)
#
# Reference environment:
#   MOLE commit b0b9ff76aff40481036d30b63a3bb00243d921ea (2026-09-03)
#   Octave 8.4.0, Python 3.12.3, numpy 2.4.4, scipy 1.17.1, SuperLU (splu)

param([string]$Only = "all")

$ErrorActionPreference = "Continue"
$Pat = 'r_x .*max|sponge:|consistency:|rms residual|max/mean|BEAM ANGLE|DIVERGED|static projection|VIOLATED|UNSTABLE|Raise --taus'
$LogDir = Join-Path (Get-Location) "repro-logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$script:N = 0

function Sec([string]$n, [string]$t) {
    if ($Only -eq "all" -or $Only -eq $n) {
        Write-Host ""; Write-Host "=== $n. $t ===" -ForegroundColor Cyan
        return $true
    }
    return $false
}

function Run([string[]]$A) {
    $cmd = @("iwbcurv.py", "--mole", $env:MOLE_SRC, "--beat", "9999") + $A
    Write-Host ("   python " + ($cmd -join " ")) -ForegroundColor DarkGray
    $script:N++
    $log = Join-Path $LogDir ("run_{0:d3}.log" -f $script:N)
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
foreach ($f in "iwbcurv.py", "test_bedrows.py", "synth_tracker.py") {
    if (Test-Path $f) { Write-Host "   OK   $f" }
    else { Write-Host "   MISSING  $f" -ForegroundColor Red; $fail = $true }
}
if (-not (Test-Path "grids")) { Write-Host "   MISSING  grids\  (run from your iwb-curv dir)" -ForegroundColor Red; $fail = $true }
else { Write-Host "   OK   grids\" }
if (-not $env:MOLE_SRC) { Write-Host "   MOLE_SRC not set" -ForegroundColor Red; $fail = $true }
elseif (-not (Test-Path (Join-Path $env:MOLE_SRC "gridGen.m"))) {
    Write-Host "   MOLE_SRC has no gridGen.m: $env:MOLE_SRC" -ForegroundColor Red
    Write-Host "   (Octave addpath does NOT error on a bad path -- it surfaces later as 'gridGen undefined')" -ForegroundColor Red
    $fail = $true
} else { Write-Host "   OK   MOLE_SRC = $env:MOLE_SRC" }
$pv = (& python -c "import sys,numpy,scipy,oct2py;print(sys.version.split()[0],numpy.__version__,scipy.__version__,oct2py.__version__)" 2>&1)
if ($LASTEXITCODE -ne 0) { Write-Host "   python/import failed: $pv" -ForegroundColor Red; $fail = $true }
else { Write-Host "   OK   python numpy scipy oct2py = $pv" }
Write-Host "   -- grid smoke test (seconds) --"
Run @("--Lx","6000","--nx","256","--nz","101","--ratio","0.6","--gridonly")
if ($fail) { Write-Host "`npreflight failed -- fix the above before running the sweep" -ForegroundColor Red; exit 1 }
if ($Only -eq "0") { exit 0 }

# ------------------------------------------------------------------ 1
if (Sec "1" "bed constraint: acceptance test, 40 periods, no --bsc") {
  foreach ($r in "0.2","0.4","0.6","0.8") {
    Write-Host "-- omega/N = $r"
    Run @("--ratio",$r,"--Lx","6000","--nx","256","--nz","101","--nper","40")
  }
}

# ------------------------------------------------------------------ 2
if (Sec "2" "A/B against the old penalty, and the controls") {
  Write-Host "-- old penalty (EXPECT: DIVERGED at t/T = 9.00)"
  Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","40","--bedmode","projection","--bsc","12")
  Write-Host "-- constraint, same grid"
  Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","40")
  foreach ($m in "post","none") {
    Write-Host "-- control: $m"
    Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","5","--bedmode",$m)
  }
  Write-Host "-- equilibration invariance (must match --bedeq auto exactly)"
  Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","40","--bedeq","unit")
}

# ------------------------------------------------------------------ 3
if (Sec "3" "sponge: the contamination fix (worst case, omega/N = 0.8)") {
  $cfg = @( @("0.10","100"), @("0.20","100"), @("0.20","50"), @("0.30","50"), @("0.30","25") )
  foreach ($c in $cfg) {
    $lsl = $c[0]; $taus = $c[1]
    Write-Host "-- lsl=$lsl taus=$taus"
    Run @("--ratio","0.8","--Lx","6000","--nx","256","--nz","101","--nper","20","--lsl",$lsl,"--taus",$taus)
  }
  Write-Host "-- sponge stability guard (EXPECT: refuses, prints *** UNSTABLE ***)"
  Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","2","--lsl","0.2","--taus","5")
}

# ------------------------------------------------------------------ 4
if (Sec "4" "corrected acceptance table, nz=101 and nz=201") {
  # ratio, Lx, nx, taus = T/30, spp ; 0.2 needs the bigger box (bounce = 4900 m)
  $cfg = @( @("0.2","12000","512","149.6","60"), @("0.4","6000","256","74.8","120"),
            @("0.6","6000","256","49.9","120"), @("0.8","6000","256","37.4","120") )
  foreach ($nz in "101","201") {
    foreach ($c in $cfg) {
      $r = $c[0]; $lx = $c[1]; $nx = $c[2]; $taus = $c[3]; $spp = $c[4]
      Write-Host "-- omega/N = $r  nz=$nz"
      Run @("--ratio",$r,"--Lx",$lx,"--nx",$nx,"--nz",$nz,"--nper","20","--spp",$spp,
            "--lsl","0.2","--taus",$taus,"--ab","20","--forcegrid")
    }
  }
}

# ------------------------------------------------------------------ 5
if (Sec "5" "RETRACTION control: vary r_x 10x at FIXED dz -- error must stay flat") {
  foreach ($a in "20","10","5","2") {
    Write-Host "-- ab = $a m"
    Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","20","--spp","120",
          "--lsl","0.2","--taus","49.9","--ab",$a)
  }
}

# ------------------------------------------------------------------ 6
if (Sec "6" "convergence: refine dx and dz together (ab=10, r_x never binds)") {
  $cfg = @( @("256","101"), @("384","151"), @("512","201") )
  foreach ($c in $cfg) {
    $nx = $c[0]; $nz = $c[1]
    Write-Host "-- nx=$nx nz=$nz"
    Run @("--ratio","0.6","--Lx","6000","--nx",$nx,"--nz",$nz,"--ab","10","--nper","20","--spp","60",
          "--lsl","0.2","--taus","49.9")
  }
}

# ------------------------------------------------------------------ 7
if (Sec "7" "does r_x > 1 break it? benchmark ridge, --forcegrid (EXPECT: VIOLATED, runs anyway)") {
  foreach ($nz in "101","151","201") {
    Write-Host "-- nz=$nz"
    Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz",$nz,"--ab","20","--nper","20","--spp","120",
          "--lsl","0.2","--taus","49.9","--forcegrid")
  }
}

# ------------------------------------------------------------------ 8
if (Sec "8" "dead ends re-tested on THIS solver (should all be flat)") {
  Write-Host "-- timestep: angle must be identical at spp 60/120/240"
  foreach ($s in "60","120","240") {
    Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","20","--spp",$s)
  }
  Write-Host "-- operator order: k = 2/4/6"
  foreach ($k in "2","4","6") {
    Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","101","--nper","20","--order",$k)
  }
}

# ------------------------------------------------------------------ 9
if (Sec "9" "offline tests (no Octave for the first, MOLE_SRC needed for the second)") {
  & python test_bedrows.py
  & python synth_tracker.py $env:MOLE_SRC
}

# ------------------------------------------------------------------ 10
if (Sec "10" "save the best configuration for view.py") {
  Run @("--ratio","0.6","--Lx","6000","--nx","256","--nz","201","--nper","20","--spp","120",
        "--lsl","0.2","--taus","49.9","--forcegrid","--out","iwb_r06_refined.npz")
  Write-Host "then:  python view.py iwb_r06_refined.npz"
}

Write-Host ""
Write-Host "$script:N runs, full logs in $LogDir" -ForegroundColor Cyan
