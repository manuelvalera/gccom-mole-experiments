# refit_fields.ps1 -- rerun the grid series saving the FIELDS, so the fit
# window can be varied afterwards for free.
#
#   .\refit_fields.ps1
#   python refit.py refit-npz
#
# The alpha sweep ruled the lateral boundary out of the beam-angle offset:
# 1e-4 and 1e-6 differ by 0.01-0.10 deg at every grid while the series shape is
# unchanged. And the seiche showed the core discretization is second order and
# converges to zero. So the remaining suspects are the ridge near field, the
# sponge, and the fit window -- and the window is pure post-processing, which
# means it can be tested on saved fields at no simulation cost.
#
# These runs are at alpha=1e-6, which is the floor-level value the seiche
# established, so the saved fields are the best available.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "refit-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }

foreach ($cfg in @("0.6","49.87","133","1000"), @("0.8","37.40","75","562")) {
  foreach ($g in @("256","101"), @("384","151"), @("512","201"), @("640","251")) {
    $npz = Join-Path $Out ("f{0}_{1}x{2}.npz" -f $cfg[0], $g[0], $g[1])
    if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
    Write-Host "-- ratio=$($cfg[0])  $($g[0])x$($g[1])  -> $npz"
    & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
        --Lx 6000 --nx $g[0] --nz $g[1] --ratio $cfg[0] `
        --nper 20 --spp 120 --lsl 0.20 --taus $cfg[1] `
        --win $cfg[2] $cfg[3] --alpha 1e-6 --forcegrid --out $npz 2>&1 |
      Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    if (-not (Test-Path $npz)) {
      Write-Host "   !! no npz written -- check --out is supported" -ForegroundColor Red }
  }
}
Write-Host ""
Write-Host "now: python refit.py refit-npz" -ForegroundColor Cyan
