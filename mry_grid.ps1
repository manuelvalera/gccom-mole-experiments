# mry_grid.ps1 -- is the Monterey growth a grid-quality problem?
#
#   copy iwbcurv.py .
#   copy mry_setup.py .
#   .\mry_grid.ps1              (~5 min on the GPU)
#
# WHAT IS ESTABLISHED
#
# Free runs (no tide, random start) on the shelf transect grow in ||u|| while
# ||w|| decays -- a horizontal flow with almost no vertical velocity, so not an
# internal wave. It:
#   needs stratification   (--N 0 decays; the profile grows)
#   is NOT the time step   (spp 600 and 2400 identical)
#   appears with refinement (192x51 and 256x51 stable; 384x76 grows; the forced
#                           runs grew faster at 768 and faster still at 1152)
#   survives --buoy energy  (the fix that cured the benchmark ridge)
#
# The suspect is the grid. The default stretching on a 115:1 domain gives cell
# angles from 6 to 167 degrees and a 15x spread in Jacobian, and the worst of it
# is at the shallow end where the local aspect ratio approaches 1000:1. The
# buoyancy exchange conserves energy in the area-weighted inner product, but the
# pressure projection is orthogonal in whatever inner product MOLE's D and G are
# adjoint in; on a near-uniform grid those coincide, and on this one they need
# not -- which would let the pressure do work on the flow.
#
# THE RUNS: same physics, four grids.
#
#   default    as built
#   gentle     milder stretching (--bt/--bb), so cells are closer to uniform
#   smooth600  bathymetry smoothed over 600 m instead of 200 m
#   trim20     transect cut where it shoals past 20 m, removing the worst cells
#
# If the growth tracks grid quality, this is a conditioning problem and the fix
# is the grid (or a better-conditioned projection). If it survives all four,
# it is in the discretization and needs the same treatment the bed coupling got.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-grid"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
python mry_setup.py --bathy narrow_bathy_100m.csv --out grids\mryshelf_s600 --smooth 600 | Out-Null
python mry_setup.py --bathy narrow_bathy_100m.csv --out grids\mryshelf_t20 --trim 20 | Out-Null

$free = @("--beat", "9999", "--grids", "grids", "--nprofile", "mry_N.txt",
          "--omega", "1.405e-4", "--taus", "1490", "--lsl", "0.15", "--forcegrid",
          "--solver", "cudss", "--device", "gpu", "--buoy", "energy",
          "--u0", "0", "--noise", "1e-3", "--probe", "--probebox", "4000", "1000", "60",
          "--nx", "384", "--nz", "76", "--nper", "20", "--spp", "600")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @free @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'grid angle|jacobian|probe t/T= +(1|10|20)\.0|GROWTH|DIVERGED' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Go "default"   @("--gridname", "mryshelf", "--Lx", "10110", "--D0", "87.8")
Go "gentle"    @("--gridname", "mryshelf", "--Lx", "10110", "--D0", "87.8", "--bt", "0.5", "--bb", "0.5")
Go "smooth600" @("--gridname", "mryshelf_s600", "--Lx", "10110", "--D0", "87.8")
Go "trim20"    @("--gridname", "mryshelf_t20")

Write-Host ""
Write-Host "compare the ||u|| column at periods 1, 10 and 20 in each run:" -ForegroundColor Cyan
Write-Host "growth there is the instability; ||w|| decays in all of them." -ForegroundColor Cyan
