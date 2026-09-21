# finer.ps1 -- two more grids at w/N=0.8, to decide between first and second
# order for the beam angle.
#
#   .\finer.ps1
#   python centroid_track.py refit-npz
#
# WHY
#
# With the centroid tracker, the near-field bias at w/N=0.8 (0.10-0.40 of a
# bounce) is +2.76, +1.45, +0.82, +0.46 on 256/384/512/640. That series fits
# two stories almost equally well, and they end in different places:
#
#     p=1   limit -1.09 deg   residual 0.02
#     p=2   limit +0.13 deg   residual 0.13
#
# The seiche is second order, but the seiche is a flat box. The beam also has
# topography whose SAMPLING changes with the grid (crest 19.0 -> 19.8 m over
# 256 -> 640), and that is a separate convergence process the seiche cannot
# see. So p=2 cannot be assumed; it has to be measured.
#
# The two stories predict different things at finer grids:
#
#     grid       dz (m)    p=1 predicts    p=2 predicts
#     768x301    3.33        +0.19           +0.43
#     896x351    2.86        +0.00           +0.35
#
# A 0.25-0.35 deg separation is well above the ~0.1 deg measurement noise, so
# two points settle it. These are the largest runs so far -- expect 15-20 and
# 25-35 minutes, and several GB of memory for the factorisation. The script
# skips any npz that already exists, so it can be stopped and resumed.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "refit-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }

foreach ($g in @("768","301"), @("896","351")) {
  $npz = Join-Path $Out ("f0.8_{0}x{1}.npz" -f $g[0], $g[1])
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
  Write-Host "-- ratio=0.8  $($g[0])x$($g[1])  -> $npz"
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx $g[0] --nz $g[1] --ratio 0.8 `
      --nper 20 --spp 120 --lsl 0.20 --taus 37.40 `
      --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
    Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|crest|splu|DIVERGED|\*\*\*|MemoryError' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
  if (-not (Test-Path $npz)) {
    Write-Host "   !! no npz -- check memory; the factorisation is the limit" -ForegroundColor Red }
}
Write-Host ""
Write-Host "now: python centroid_track.py refit-npz" -ForegroundColor Cyan
Write-Host "compare the 0.10-0.40 centroid row at 768 and 896 with +0.19/+0.00 (p=1)" -ForegroundColor Cyan
Write-Host "and +0.43/+0.35 (p=2)." -ForegroundColor Cyan
