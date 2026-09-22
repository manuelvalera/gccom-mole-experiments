# ramp_test.ps1 -- does a smooth start shorten the spin-up, on the grid where
# spin-up is worst?
#
#   .\ramp_test.ps1                       (~35 min with PARDISO, 4 threads)
#   python time_compare.py time896        (unramped: 20, 40, 80 periods)
#   python time_compare.py time896r       (ramped:   25, 40 periods)
#
# The near-field centroid bias moved +0.14 deg between 20 and 40 periods at
# 512x201 but +0.52 deg at 896x351, so the spin-up transient lasts longer on
# finer grids -- and it biased the six-grid convergence series. The tide was
# switched on abruptly; --ramp 3 brings it on smoothly over three periods.
#
# Needed: an unramped 80-period run as the settled reference (the 20- and
# 40-period ones already exist in time896), and ramped runs at 25 and 40.
# If ramped-40 (or even ramped-25) matches unramped-80 in the cent .10-.40
# column, the settled six-grid series can be run at that length instead of 80.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
if (-not $env:MKL_NUM_THREADS) { $env:MKL_NUM_THREADS = "4" }
foreach ($d in "time896", "time896r") {
  if (-not (Test-Path $d)) { New-Item -ItemType Directory -Path $d | Out-Null } }

function Go([string]$npz, [string[]]$extra) {
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; return }
  Write-Host "-- $npz  $($extra -join ' ')" -ForegroundColor Cyan
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx 896 --nz 351 --ratio 0.8 --ab 20 --spp 120 `
      --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --forcegrid `
      --solver pardiso --out $npz @extra 2>&1 |
    Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
    Select-String -Pattern 't/T=|ms/step|BEAM ANGLE|max/mean|DIVERGED|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Go "time896\ab20_p80.npz"  @("--nper", "80", "--ramp", "0")
Go "time896r\ab20_p25.npz" @("--nper", "25", "--ramp", "3")
Go "time896r\ab20_p40.npz" @("--nper", "40", "--ramp", "3")

Write-Host ""
Write-Host "now: python time_compare.py time896   and   python time_compare.py time896r" -ForegroundColor Cyan
