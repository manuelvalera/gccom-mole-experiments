# mry_bore.ps1 -- the shoaling internal wave on the real Monterey slope.
#
#   .\mry_bore.ps1              (~15 min on the GPU)
#
# THE TARGET
#
# Walter, Woodson, Arthur, Fringer & Monismith (2012), JGR 117 C07017,
# Figure 10, non-canonical case: a mode-1 internal wave shoaling on the actual
# southern Monterey Bay slope, which they model with SUNTANS on a 20 km x 80 m
# transect at dx 5 m, dz 1 m, forced at the offshore boundary at M2, with an
# eddy viscosity of 1e-4 m^2/s and no scalar diffusivity. They run six tidal
# periods (2.7e5 steps at dt = 1 s).
#
# Their distinction is the internal Iribarren number xi = s / sqrt(a/lambda).
# For this stratification mode1.py gives c = 0.209 m/s and lambda = 9.4 km, and
# on the real slope (s ~ 0.04) xi ~ 2 needs an offshore amplitude a ~ 3.7 m --
# the non-canonical regime, where the wave does not have room to break and
# instead surges up and down the slope.
#
# SETUP
#
#   grid     the measured transect extended offshore to 20 km at 87.8 m, cut
#            where it shoals past 10 m (the cells there are 1000:1 and host an
#            instability; see ROADMAP)
#   forcing  mode-1 wave at the offshore boundary, --mode1, no barotropic tide
#   physics  full nonlinear advection, rotation at 36.8 N, nu = 1e-4
#
# SCALAR DIFFUSIVITY. Walter et al. use none, letting their scheme set it at
# dx = 5 m and dt = 1 s. At our resolution the buoyancy field steepens at the
# run-up until max|b| jumps two orders of magnitude in fifty steps and the run
# dies -- while max|u| is still falling, so it is the b equation alone.
# --kappa 1e-4, the same value as the viscosity, carries it through. This is a
# deviation from their setup and belongs in any writeup; the alternative is
# their resolution, which needs the implicit viscous term first.
#
# TIME STEP. Explicit diffusion needs nu*dt*lambda_max below 2, and lambda_max
# lives in the thinnest cells: 36 1/m^2 at 192x51 but 321 at 768x151. At
# spp 600 that gives 2.4 and the run dies at step 40 whatever the amplitude.
# spp 2400 brings it to 0.6. Walter et al. never met this because their dt is
# 1 s against our 19 s -- the M2 period sets our clock, so the same physical
# viscosity is four orders of magnitude harder to carry explicitly.
#
# The shoreward sponge is 0.05 of the domain on purpose: wider absorbs the
# shoaling wave that is the subject of the experiment. At a = 3.7 m that was
# enough to blow up before the buoyancy ghost cells were fixed; with upwind2 it
# now runs.
#
# A  amplitude: 1 m (weak, should stay linear), 3.7 m (xi ~ 2), 7 m (stronger)
# B  scheme at the working amplitude. NOTE: on this grid the LIMITED schemes
#    (minmod, vanleer) go unstable on the shoreward slope around t/T = 3, while
#    upwind1 and upwind2 run clean. Both limiters read the downwind cell and
#    switch order cell by cell, and on a strongly terrain-following grid the
#    along-xi gradient of b carries the isopycnal tilt, so the switching varies
#    face to face. upwind2 is the working choice for now; the limiter question
#    has to be settled before the bore case, since a bore needs monotonicity.
#
# Watch max|u| and the animations: a wave that steepens and surges is the
# result; one that simply reflects is not.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-bore"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
python mry_setup.py --bathy narrow_bathy_100m.csv --out grids\mry20 --extend 20 --trim 10 | Out-Null

$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mry20",
            "--Lx", "19900", "--D0", "87.8", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--lat", "36.8", "--forcegrid",
            "--solver", "cudss", "--device", "gpu", "--buoy", "energy",
            "--advect", "full", "--nu", "1e-4", "--kappa", "1e-4", "--u0", "0",
            "--nx", "768", "--nz", "151", "--spp", "2400", "--nper", "8",
            "--ramp", "1", "--taus", "1500", "--lsl", "0.20", "--lslr", "0.05",
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
    $out | Select-String -Pattern 'mode 1:|viscosity :|t/T=|time loop:|max/mean|DIVERGED|Traceback|Error' |
      Select-Object -Last 10 | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. amplitude ===" -ForegroundColor Cyan
Go "a1.0" @("--mode1", "1.0", "--advscheme", "upwind2", "--out", (Join-Path $Out "a1.0.npz"))
Go "a3.7" @("--mode1", "3.7", "--advscheme", "upwind2", "--out", (Join-Path $Out "a3.7.npz"))
Go "a7.0" @("--mode1", "7.0", "--advscheme", "upwind2", "--out", (Join-Path $Out "a7.0.npz"))

Write-Host "=== B. scheme at xi ~ 2 ===" -ForegroundColor Cyan
Go "a3.7_mm" @("--mode1", "3.7", "--advscheme", "minmod", "--out", (Join-Path $Out "a3.7_mm.npz"))

Write-Host "=== C. frames for the animation ===" -ForegroundColor Cyan
Go "frames" @("--mode1", "3.7", "--advscheme", "upwind2", "--frames", "fr_bore", "--framerate", "16")

Write-Host ""
Write-Host "python animate.py fr_bore -o bore.gif --field b --fps 12" -ForegroundColor Cyan
Write-Host "python view.py $Out\a3.7.npz --rms --log -o bore_rms.png" -ForegroundColor Cyan
