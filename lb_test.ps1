# lb_test.ps1 -- is the first-order beam-angle error the resolution of a
# NARROW beam?
#
#   .\lb_test.ps1
#   python ab_compare.py 30=ab5-npz 60=lb60-npz      (labels are ridge WIDTHS here)
#
# WHAT THE RIDGE-HEIGHT TEST LEFT
#
# Near-field limits (centroid, 0.10-0.40 of a bounce) at ab = 5, 10, 20 m are
# -0.69, -0.77, -1.05 deg, and the window dependence shrinks with the ridge:
# at ab = 5 the near and long windows agree (-0.69 vs -0.72). So the curvature
# and part of the offset scale with ridge height, but two things do NOT:
#
#   1. a residual near -0.7 deg that survives as the ridge shrinks, and
#   2. the first-order error slope: (bias at 512 - bias at 896) / (dz 5 - 2.857)
#      is 0.39 deg/m at ab = 5 and 0.37 deg/m at ab = 20.
#
# An O(h) error whose size does not depend on ridge height cannot come from
# discretizing the ridge's HEIGHT. It could come from its WIDTH: the beam is
# born with a cross-beam scale set by Lb = 30 m, which the grids resolve with
# only a handful of cells, while the seiche mode that converged cleanly at
# second order spans dozens.
#
# THE TEST
#
# Double the ridge width, and double its height with it so the maximum slope
# is unchanged (ab/Lb = 5/30 = 10/60) -- the topography stays equally gentle,
# only wider. If the first-order slope roughly halves and the limit moves, the
# narrow beam's resolution is the culprit. If neither moves, it is not.
#
# The earlier Lb = 120 attempt failed because the wide ridge radiated too
# weakly to track (rms residual 46 m, max/mean 5.2). Lb = 60 is a smaller step;
# check the solver's max/mean stays above ~8 before trusting the numbers.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "lb60-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
foreach ($g in @("512","201"), @("640","251"), @("896","351")) {
  $npz = Join-Path $Out ("f0.8_{0}x{1}.npz" -f $g[0], $g[1])
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
  Write-Host "-- Lb=60 ab=10  $($g[0])x$($g[1])"
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx $g[0] --nz $g[1] --ratio 0.8 --ab 10 --lb 60 `
      --nper 20 --spp 120 --lsl 0.20 --taus 37.40 `
      --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
    Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
    Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}
Write-Host ""
Write-Host "now: python ab_compare.py 30=ab5-npz 60=lb60-npz" -ForegroundColor Cyan
