# mry_diag.ps1 -- why is the Monterey run growing?
#
#   copy iwbcurv.py .
#   .\mry_diag.ps1              (~5 min)
#
# THE SYMPTOM
#
# max|u| climbs through the run -- 2.86e-2, 3.07e-2, 4.33e-2 at periods 6, 8,
# 10 on 768x151, and 2.94e-2, 5.87e-2, 7.89e-2 on 1152x226. Faster on the finer
# grid, unchanged by a 3x smaller time step. That is a spatial instability,
# the same class as the one the buoyancy fix cured on the benchmark ridge --
# except --buoy energy is already on here.
#
# THE RUNS
#
# Free runs: no tide, a small random initial field, so any growth is the
# instability by itself. The probe box sits over the steep upper slope
# (x ~ 3.5-4.5 km), which is where the transect goes supercritical.
#
#   base        the configuration as run
#   scalarN     constant N instead of the measured profile -- does the
#               depth-varying stratification drive it?
#   N0          no stratification at all -- is it in the buoyancy coupling
#               again, or in the pressure/bed side?
#   logical     the OLD buoyancy coupling, for contrast
#   coarse/fine resolution dependence, which the forced runs already suggest
#
# NOTE. The first attempt at this was undermined by three instrumentation bugs,
# now fixed: --N 0 crashed in the theory-angle calculation; the blow-up
# threshold was 50x the NOISE amplitude, while projecting noise onto a 115:1
# grid starts four times above that, so runs were flagged as diverged when they
# were not; and the probe box measured height above the domain's deepest point,
# so over the sloping bed at x ~ 4 km it sat below the seabed and sampled the
# wrong place. Delete mry-diag\ before rerunning.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-diag"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$free = @("--beat", "9999", "--grids", "grids", "--gridname", "mryshelf",
          "--Lx", "10110", "--D0", "87.8", "--omega", "1.405e-4",
          "--taus", "1490", "--lsl", "0.15", "--forcegrid",
          "--solver", "cudss", "--device", "gpu",
          "--u0", "0", "--noise", "1e-3", "--probe", "--probebox", "4000", "1000", "60",
          "--nper", "20", "--spp", "600")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @free @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'GROWTH|DIVERGED|stiffness|probe box|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Go "base"      @("--nx", "768", "--nz", "151", "--nprofile", "mry_N.txt", "--buoy", "energy")
Go "transport" @("--nx", "768", "--nz", "151", "--nprofile", "mry_N.txt", "--buoy", "energy",
                 "--btmode", "transport")
Go "scalarN" @("--nx", "768", "--nz", "151", "--N", "4e-3", "--buoy", "energy")
Go "N0"      @("--nx", "768", "--nz", "151", "--N", "0", "--buoy", "energy")
Go "logical" @("--nx", "768", "--nz", "151", "--nprofile", "mry_N.txt", "--buoy", "logical")
Go "coarse"  @("--nx", "384", "--nz", "76", "--nprofile", "mry_N.txt", "--buoy", "energy")
Go "fine"    @("--nx", "1152", "--nz", "226", "--nprofile", "mry_N.txt", "--buoy", "energy")

Write-Host ""
Write-Host "summary:" -ForegroundColor Cyan
Get-ChildItem "$Out\*.log" | ForEach-Object {
    $r = Select-String -Path $_.FullName -Pattern 'GROWTH  ridge' | Select-Object -Last 1
    $d = Select-String -Path $_.FullName -Pattern 'DIVERGED at t/T=([0-9.]+)' | Select-Object -First 1
    $dtxt = if ($d) { "DIVERGED t/T=" + $d.Matches[0].Groups[1].Value } else { "ran to end" }
    $rtxt = if ($r) { ($r.Line -replace '^\s*GROWTH\s+ridge \|\|w\|\|\s+', '') } else { "" }
    Write-Host ("   {0,-14} {1,-22} {2}" -f $_.Name, $dtxt, $rtxt)
}
