# beam_nonlinear.ps1 -- the internal-wave beam, nonlinearly.
#
#   .\beam_nonlinear.ps1        (~10 min on the GPU)
#   python centroid_track.py beam-nl-npz
#
# THE QUESTION
#
# The 2021 mimetic GCCOM reproduced the lock exchange but could not sustain an
# internal-wave beam: section 4.2.3 of the dissertation reports beams forming
# at the right angle and then dissipating, in both 2-D and 3-D, and attributes
# it to the first-order one-sided advection (d2D/d3D are [-1, 1]/dx).
#
# This repository already passes that case LINEARLY -- beams at the theoretical
# angle, converging under refinement -- which isolates the rest of the
# formulation as sound. The remaining question is whether the advection scheme
# is what kills them. The lock release just measured the cost of first order on
# a case with a known answer: Fr 0.674 against 0.699 for second order.
#
# Here the same comparison is run on the beam itself:
#
#   linear     the verified reference, no advection at all
#   upwind1    first order, i.e. the 2021 scheme
#   upwind2    second order
#   minmod     second order away from extrema, first order at them
#
# The numerical diffusion of first-order upwind is |u| dx / 2, which at these
# beam velocities and this grid is ~3e-2 m^2/s -- four orders of magnitude above
# anything physical. If that is what dissipated the beams, upwind1 loses them
# and upwind2 keeps them.
#
# What to look at afterwards: max/mean in each log (a coherent beam sits at
# 12-25, reflections at ~5), and the angles from centroid_track.py.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "beam-nl-npz"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "9999", "--Lx", "6000", "--nx", "512", "--nz", "201",
            "--ratio", "0.8", "--nper", "25", "--spp", "120", "--lsl", "0.20",
            "--taus", "37.40", "--win", "75", "562", "--alpha", "1e-6",
            "--ramp", "3", "--buoy", "energy", "--forcegrid",
            "--solver", "cudss", "--device", "gpu")

function Go([string]$tag, [string[]]$a) {
    $npz = Join-Path $Out "f0.8_$tag.npz"
    $log = $npz -replace '\.npz$', '.log'
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a --out $npz 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'max/mean|BEAM ANGLE|ms/step|DIVERGED' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Go "linear"  @("--advect", "none")
Go "upwind1" @("--advect", "full", "--advscheme", "upwind1")
Go "upwind2" @("--advect", "full", "--advscheme", "upwind2")
Go "minmod"  @("--advect", "full", "--advscheme", "minmod")

Write-Host ""
Write-Host "max/mean far below the linear case = the beam is being dissipated." -ForegroundColor Cyan
Write-Host "then: python centroid_track.py $Out" -ForegroundColor Cyan
