# speed_check.ps1 -- does the cache + PARDISO actually speed things up here,
# and does PARDISO give the same answers?
#
#   .\speed_check.ps1                 (~15 min; the first run pays Octave once)
#
# Prints per-step time for SuperLU and for PARDISO at several thread counts on
# the 896x351 problem, then checks PARDISO against the seiche value your
# machine already measured with SuperLU (256x101, mode 4,2, alpha 1e-6:
# relative error -4.0844e-04). A solver that changes that number is not faster,
# it is wrong.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'grid from cache|operators from cache|generating grid|splu|pardiso factorise|static projection|ms/step|relative error|BEAM ANGLE|DIVERGED|Error|\*\*\*'
$beam = @("--beat","9999","--Lx","6000","--nx","896","--nz","351","--ratio","0.8",
          "--spp","120","--lsl","0.20","--taus","37.40","--win","75","562","--forcegrid")
$oldThreads = $env:MKL_NUM_THREADS

function Go([string]$label, [string[]]$a) {
  Write-Host "-- $label" -ForegroundColor Cyan
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC @a 2>&1 |
    Select-String -Pattern $Pat | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   wall {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Write-Host "=== 0. warm the cache (Octave runs once, ~3 min) ===" -ForegroundColor Cyan
Go "cache fill, superlu, 1 period" ($beam + @("--nper","1","--solver","superlu"))

Write-Host "=== 1. per-step cost, 3 periods, cache warm ===" -ForegroundColor Cyan
Go "superlu" ($beam + @("--nper","3","--solver","superlu"))
foreach ($th in "4","8","16") {
  $env:MKL_NUM_THREADS = $th
  Go "pardiso, $th threads" ($beam + @("--nper","3","--solver","pardiso"))
}
$env:MKL_NUM_THREADS = $oldThreads

Write-Host "=== 2. correctness: seiche, 256x101, mode (4,2), alpha 1e-6 ===" -ForegroundColor Cyan
Write-Host "   must reproduce the SuperLU value -4.0844e-04" -ForegroundColor Cyan
$seiche = @("--beat","9999","--Lx","6000","--nx","256","--nz","101","--ab","0","--bulge","0",
            "--seiche","4","2","--nper","6","--spp","400","--alpha","1e-6")
Go "seiche, superlu" ($seiche + @("--solver","superlu"))
Go "seiche, pardiso" ($seiche + @("--solver","pardiso"))
