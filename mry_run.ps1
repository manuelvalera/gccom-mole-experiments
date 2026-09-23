# mry_run.ps1 -- first forced runs on the Monterey shelf transect.
#
#   python mry_setup.py --bathy narrow_bathy_100m.csv --out grids\mryshelf
#   .\mry_run.ps1
#   python view.py mry-npz\mry_10p.npz --rms --log -o mry_first.png
#
# CONFIGURATION, and why each number is what it is
#
#   --gridname mryshelf   the measured transect: 10.1 km long, 8 to 88 m deep
#   --nprofile mry_N.txt  N varies 2.1e-3 to 1.5e-2 over the column, a factor of
#                         seven; a scalar N cannot represent it
#   --omega 1.405e-4      M2, instead of the benchmark's omega = ratio * N
#   --spp 600             set by STIFFNESS, not by the forcing period: the
#                         buoyancy update needs N_max*dt < 1, and N_max*dt is
#                         1.14 here. At the benchmark's 120 steps it is 5.7 and
#                         the run blows up in a fifth of a period. Checked: 600
#                         and 2400 steps agree to 0.2%.
#   --taus 1490           T/30, the same ratio the benchmark used; the sponge
#                         relaxation is explicit, so dt/taus must stay below 1
#   --ramp 3              smooth start (see the spin-up work)
#   --buoy energy         the energy-consistent bed coupling; this transect has
#                         slopes up to 0.074, so the leak would be active
#
# STILL MISSING: rotation. At 36.8 N, f/omega = 0.62 and the beam slope is 22%
# shallower than what this computes. These runs are the plumbing check, not
# physics to interpret.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-npz"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mryshelf",
            "--Lx", "10110", "--D0", "87.8", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--taus", "1490", "--lsl", "0.15",
            "--ramp", "3", "--buoy", "energy", "--forcegrid",
            "--solver", "cudss", "--device", "gpu")

function Go([string]$tag, [string[]]$a) {
    $npz = Join-Path $Out "$tag.npz"
    $log = $npz -replace '\.npz$', '.log'
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $t0 = Get-Date
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a --out $npz 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'stiffness|lepticity|r_x|consistency: div|t/T=|ms/step|max/mean|DIVERGED|\*\*\*' |
      Select-Object -Last 8 | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Go "mry_10p"      @("--nx", "768", "--nz", "151", "--nper", "10", "--spp", "600")
Go "mry_10p_fine" @("--nx", "1152", "--nz", "226", "--nper", "10", "--spp", "600")
Go "mry_10p_dt"   @("--nx", "768", "--nz", "151", "--nper", "10", "--spp", "1800")

Write-Host ""
Write-Host "compare mry_10p with mry_10p_dt (time step) and mry_10p_fine (grid)." -ForegroundColor Cyan
Write-Host "then: python view.py $Out\mry_10p.npz --rms --log -o mry_first.png" -ForegroundColor Cyan
