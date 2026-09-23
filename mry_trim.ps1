# mry_trim.ps1 -- settle the shallow-end cut, then check the trimmed domain
# behaves under refinement.
#
#   .\mry_trim.ps1              (~8 min on the GPU)
#
# WHAT IS ESTABLISHED
#
# On the full transect, free runs grow in ||u|| (0.105 -> 0.47 over 20 periods)
# while ||w|| decays. Milder stretching does not help; smoothing the bathymetry
# does not help; trimming the transect where it shoals past 20 m does, and
# ||u|| then decays (0.105 -> 0.068). The Jacobian spread falls from 15x to
# 2.5x at the same time.
#
# The cause is structural: with nz fixed, dz = D(x)/nz, so 8 m of water gives
# cells 0.1 m tall and 26 m wide -- roughly 250:1 within a single cell, against
# 15:1 in the deep part. Refining makes those cells thinner rather than better,
# which is why the instability strengthened with resolution.
#
# SECTION A: how little can be trimmed? Each metre of trim costs transect, and
# the shallow end is where the shoaling internal tide would break, so cutting
# more than necessary is not free.
#
# SECTION B: with the cut chosen, does the trimmed domain converge? The full
# domain got worse with refinement; the trimmed one should not. These are
# FORCED runs with transport-conserving barotropic forcing, i.e. the
# configuration to use from here.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-trim"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
foreach ($d in 10, 15, 20) {
    python mry_setup.py --bathy narrow_bathy_100m.csv --out "grids\mryshelf_t$d" --trim $d | Out-Null
}

$common = @("--beat", "9999", "--grids", "grids", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--taus", "1490", "--lsl", "0.15", "--forcegrid",
            "--solver", "cudss", "--device", "gpu", "--buoy", "energy", "--spp", "600")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'jacobian|probe t/T= +(1|20)\.0|GROWTH|t/T= *(10|20)\.0 |max/mean|DIVERGED|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. free runs: how little trim is enough? ===" -ForegroundColor Cyan
foreach ($d in 10, 15, 20) {
    Go "free_t$d" @("--gridname", "mryshelf_t$d", "--nx", "384", "--nz", "76",
                    "--u0", "0", "--noise", "1e-3", "--probe",
                    "--probebox", "0", "20000", "500", "--nper", "20")
}

Write-Host "=== B. forced runs on the chosen cut, three resolutions ===" -ForegroundColor Cyan
foreach ($g in @("384", "76"), @("768", "151"), @("1152", "226")) {
    Go "forced_$($g[0])" @("--gridname", "mryshelf_t20", "--nx", $g[0], "--nz", $g[1],
                           "--btmode", "transport", "--ramp", "3", "--nper", "20",
                           "--out", (Join-Path $Out "forced_$($g[0]).npz"))
}

Write-Host ""
Write-Host "A: pick the smallest trim whose ||u|| decays." -ForegroundColor Cyan
Write-Host "B: max|u| should agree across the three grids, not grow with them." -ForegroundColor Cyan
Write-Host "then: python view.py $Out\forced_768.npz --rms --log -o mry_trimmed.png" -ForegroundColor Cyan
