# mry_rotate.ps1 -- the Monterey transect with rotation, which is where the
# physics finally is.
#
#   copy iwbcurv.py .
#   .\mry_rotate.ps1            (~10 min on the GPU)
#
# ROTATION IS VERIFIED
#
# The seiche now carries the rotating dispersion relation,
#   omega^2 = (N^2 kx^2 + f^2 pz^2) / (kx^2 + pz^2),
# and the solver matches it at every latitude with the SAME relative error as
# without rotation (-1.81e-3 at 128x51, from 0 to 85 degrees), converging at
# second order (2.10, 2.21). The (u, v) pair is advanced by an exact rotation
# through f*dt, so the Coriolis terms add no error of their own.
#
# At Monterey f/omega_M2 = 0.62, so the beam slope is 22% shallower than
# everything computed so far. Every run below therefore differs from the
# earlier ones physically, not just in detail.
#
#   A  with and without rotation, same everything else: how much does it move?
#   B  grid convergence with rotation on
#   C  a long run for the animation

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-rot"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mryshelf",
            "--Lx", "10110", "--D0", "87.8", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--forcegrid", "--solver", "cudss", "--device", "gpu",
            "--buoy", "energy", "--spp", "600", "--ulim", "5000",
            "--lsl", "0.30", "--taus", "450", "--btmode", "transport", "--ramp", "3",
            "--probe", "--probebox", "0", "20000", "500")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'rotation:|theory:|probe t/T= +16\.0|max/mean|DIVERGED' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. does rotation change the answer? ===" -ForegroundColor Cyan
Go "norot"  @("--nx", "768", "--nz", "151", "--nper", "16",
              "--out", (Join-Path $Out "norot.npz"))
Go "rot"    @("--nx", "768", "--nz", "151", "--nper", "16", "--lat", "36.8",
              "--out", (Join-Path $Out "rot.npz"))

Write-Host "=== B. convergence with rotation ===" -ForegroundColor Cyan
foreach ($g in @("384", "76"), @("1152", "226")) {
    Go "rot_$($g[0])" @("--nx", $g[0], "--nz", $g[1], "--nper", "16", "--lat", "36.8",
                        "--out", (Join-Path $Out "rot_$($g[0]).npz"))
}

Write-Host "=== C. frames for the animation ===" -ForegroundColor Cyan
Go "frames" @("--nx", "768", "--nz", "151", "--nper", "12", "--lat", "36.8",
              "--frames", "frames_mry_rot", "--framerate", "8")

Write-Host ""
Write-Host "A: compare rms|w| and max/mean between norot and rot." -ForegroundColor Cyan
Write-Host "B: rms|w| should agree across 384, 768 and 1152 (area-weighted now)." -ForegroundColor Cyan
Write-Host "C: python animate.py frames_mry_rot -o monterey_rot.gif --field speed --fps 10" -ForegroundColor Cyan
