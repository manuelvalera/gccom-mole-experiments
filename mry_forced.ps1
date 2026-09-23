# mry_forced.ps1 -- why do the forced runs blow up at t/T = 6.3?
#
#   copy iwbcurv.py .
#   .\mry_forced.ps1             (~10 min on the GPU)
#
# WHAT IS ESTABLISHED
#
# Free runs on the trimmed transect are stable at every trim depth tested
# (||u|| decays, 0.105 -> 0.069 over 20 periods). Forced runs diverge at
# t/T = 6.8, 6.3, 6.3 on 384, 768 and 1152 -- the SAME time at every
# resolution, which a grid instability would not do. The solver is linear, so
# the threshold scales with u0 and the divergence time is independent of it.
#
# The suspicion is resonance: the first baroclinic mode travels at roughly
# 0.14 m/s here, so an M2 wavelength is about 6 km in a 9.4 km domain, and if
# the sponges do not absorb that energy it simply accumulates.
#
#   SECTION A: let it run past the cutoff (--ulim 5000) and look at the SHAPE
#   of ||u||(t). Linear-in-time growth means accumulation with too little
#   damping; exponential means an instability.
#
#   SECTION B: vary the absorption. If growth is resonant, a wider or faster
#   sponge should cap it; if it is an instability, sponges will not help.
#
#   SECTION C: a domain long enough to be off-resonance (the untrimmed
#   transect is 10.1 km, the trimmed 9.4 km) plus a case with the sponge
#   covering a third of the domain.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-forced"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
python mry_setup.py --bathy narrow_bathy_100m.csv --out grids\mryshelf_t10 --trim 10 | Out-Null

$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mryshelf_t10",
            "--nprofile", "mry_N.txt", "--omega", "1.405e-4", "--forcegrid",
            "--solver", "cudss", "--device", "gpu", "--buoy", "energy",
            "--btmode", "transport", "--ramp", "3", "--nx", "384", "--nz", "76",
            "--spp", "600", "--probe", "--probebox", "0", "20000", "500",
            "--ulim", "5000")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'probe t/T= +(2|6|10|14|18)\.0|GROWTH|DIVERGED|max/mean' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. the shape of the growth ===" -ForegroundColor Cyan
Go "shape" @("--taus", "1490", "--lsl", "0.15", "--nper", "18")

Write-Host "=== B. absorption ===" -ForegroundColor Cyan
Go "sponge_wide"   @("--taus", "1490", "--lsl", "0.30", "--nper", "18")
Go "sponge_fast"   @("--taus", "450", "--lsl", "0.15", "--nper", "18")
Go "sponge_both"   @("--taus", "450", "--lsl", "0.30", "--nper", "18")

Write-Host "=== C. is it the domain length? ===" -ForegroundColor Cyan
Go "untrimmed_ref" @("--gridname", "mryshelf", "--Lx", "10110", "--D0", "87.8",
                     "--taus", "450", "--lsl", "0.30", "--nper", "18")

Write-Host ""
Write-Host "||u|| roughly linear in t -> accumulation; roughly exponential -> instability." -ForegroundColor Cyan
Write-Host "If the sponge cases saturate, the configuration needs more absorption." -ForegroundColor Cyan
