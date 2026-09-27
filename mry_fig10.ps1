# mry_fig10.ps1 -- the three gaps: mooring site, resolution, diffusivity.
#
#   .\mry_fig10.ps1             (~35 min on the GPU)
#
# WHERE WE ARE
#
# On the measured transect, forced with a mode-1 wave at the offshore boundary,
# the thermistor traces show the wave arriving surface-intensified in deep
# water and becoming bottom-intensified as it shoals, with the 2-mab/6-mab
# amplitude ratio crossing 1 between the 80 and 75 m isobaths and reaching 1.68
# at 60 m. Ranges are 0.4-0.8 degC against the 0.5-1 degC Walter et al. observed.
# The cold excursions are sharp and bottom-intensified, the warm recoveries
# gentle -- the asymmetry their Figure 10 turns on.
#
# THREE THINGS THAT ARE NOT YET RIGHT
#
# 1 MOORING SITE. Their instruments sit at the 15 m isobath; ours lands inside
#   the shoreward sponge, so the series there is damped. Section A narrows that
#   sponge and checks the run stays clean.
# 2 RESOLUTION. They use dx = 5 m, dz = 1 m. We are at 26 m by 0.6 m. Section B
#   doubles the grid, which is now affordable because diffusion is sub-cycled
#   rather than setting the global time step.
# TIME STEP. Sub-cycling freed the DIFFUSIVE limit, not the advective one. At
# 768x151 the working step is dt = 18.6 s (spp 2400); at 74.5 s every run dies
# near t/T = 3, which is when the wave reaches the run-up -- so the constraint
# is the steepening front, not viscosity. dt therefore still scales with the
# grid: spp 2400 at 768, 4800 at 1536. What sub-cycling buys is that the
# diffusion no longer forces dt DOWN FURTHER on top of that.
#
# 3 DIFFUSIVITY. They use none. We needed kappa = 1e-4 or the buoyancy field
#   steepened to the grid scale at the run-up and the run died. Section C
#   lowers it on the finer grid, where the physical gradients are resolved
#   better and less is needed. If kappa can go to zero at 1536, the deviation
#   disappears rather than being documented.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-fig10"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mry20",
            "--Lx", "19900", "--D0", "87.8", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--lat", "36.8", "--forcegrid",
            "--solver", "cudss", "--device", "gpu", "--buoy", "energy",
            "--advect", "full", "--advscheme", "upwind2", "--nu", "1e-4",
            "--u0", "0", "--mode1", "3.7", "--nper", "8", "--ramp", "1",
            "--taus", "1500", "--lsl", "0.20", "--ulim", "5000")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'diffusion  :|viscosity  :|t/T= +8\.0|DIVERGED|time loop:' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. reach the 15 m isobath: narrow the shoreward sponge ===" -ForegroundColor Cyan
Go "spg_050" @("--nx", "768", "--nz", "151", "--spp", "2400", "--kappa", "1e-4",
               "--lslr", "0.050", "--frames", "fr_s050", "--framerate", "16")
Go "spg_020" @("--nx", "768", "--nz", "151", "--spp", "2400", "--kappa", "1e-4",
               "--lslr", "0.020", "--frames", "fr_s020", "--framerate", "16")
Go "spg_010" @("--nx", "768", "--nz", "151", "--spp", "2400", "--kappa", "1e-4",
               "--lslr", "0.010", "--frames", "fr_s010", "--framerate", "16")

Write-Host "=== B. resolution ===" -ForegroundColor Cyan
Go "fine" @("--nx", "1536", "--nz", "301", "--spp", "4800", "--kappa", "1e-4",
            "--lslr", "0.020", "--frames", "fr_fine", "--framerate", "16")

Write-Host "=== C. can the diffusivity come off on the finer grid? ===" -ForegroundColor Cyan
# five periods is enough to see whether it survives the run-up, and these are
# the slowest runs in the set
Go "fine_k25" @("--nx", "1536", "--nz", "301", "--spp", "4800", "--kappa", "2.5e-5",
                "--lslr", "0.020", "--nper", "5")
Go "fine_k0" @("--nx", "1536", "--nz", "301", "--spp", "4800", "--kappa", "0",
               "--lslr", "0.020", "--nper", "5")

Write-Host ""
Write-Host "A: whichever sponge survives, read the 15 m mooring from its frames:" -ForegroundColor Cyan
Write-Host "   python thermistors.py fr_s020 --isobath 15 --mab 2 4 6 --sponge 0.02 -o th15.png" -ForegroundColor Cyan
Write-Host "B: compare th15 from fr_fine against fr_s020 -- same physics, twice the grid." -ForegroundColor Cyan
Write-Host "C: if fine_k0 survives 8 periods, the diffusivity deviation is gone." -ForegroundColor Cyan
