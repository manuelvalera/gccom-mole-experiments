# depth_test.ps1 -- is the beam curvature caused by approaching the reflection?
#
#   .\depth_test.ps1
#   python curvature.py depth-npz
#
# THE HYPOTHESIS
#
# At w/N=0.8 the sag across the fit window is a steady -50 m at every grid from
# 256x101 to 640x251, so it is not a discretization artifact. But the beam
# descends tan(53.13) * 638 = 850 m over 0.85 of a bounce in a 1000 m domain,
# which puts the outer bands close to the bottom reflection. Where the incident
# and reflected beams overlap, the column-wise maximum stops tracking the
# incident beam and migrates -- which would look exactly like curvature.
#
# THE TEST
#
# Double the depth. One bounce is D0 / tan(theta), so D0 = 2000 m doubles the
# bounce to 1500 m at w/N=0.8 while leaving the physics, the ridge and the
# beam angle unchanged. dz is held at 10 m (nz 201 -> 401) so the beam is
# resolved identically.
#
#   If the curvature is a reflection artifact, the local-angle profile plotted
#   against BOUNCE FRACTION is unchanged, and against absolute distance it
#   stretches. The sag c*L^2 roughly quadruples, since L doubles.
#
#   If the curvature is a property of the beam itself, it tracks absolute
#   distance instead, and the sag at a given |x| in metres is the same.
#
# Both depths are run at the same dx and dz so nothing else moves.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "depth-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }

# D0, nz, taus (= T/30, unchanged by depth), window in metres
foreach ($cfg in @("1000","201","37.40"), @("2000","401","37.40")) {
  $npz = Join-Path $Out ("d{0}_0.8.npz" -f $cfg[0])
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
  Write-Host "-- D0=$($cfg[0]) m  nz=$($cfg[1])  (dz = 10 m either way)"
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx 512 --nz $cfg[1] --D0 $cfg[0] --ratio 0.8 `
      --nper 20 --spp 120 --lsl 0.20 --taus $cfg[2] `
      --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
    Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|\*\*\*|r_x|Jacobian' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  if (-not (Test-Path $npz)) {
    Write-Host "   !! no npz -- check --D0 is accepted" -ForegroundColor Red }
}
Write-Host ""
Write-Host "curvature.py reads D0 from the grid, so the bounce is rescaled" -ForegroundColor Cyan
Write-Host "automatically. Compare section 1 between the two depths." -ForegroundColor Cyan
Write-Host "now: python curvature.py depth-npz" -ForegroundColor Cyan
