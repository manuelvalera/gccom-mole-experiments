# time_test.ps1 -- are the beam angles converged in TIME?
#
#   .\time_test.ps1
#   Select-String -Path time-npz\*.log -Pattern "t/T="          # max|u| history
#   python time_compare.py
#
# WHY
#
# The rms field is accumulated over the LAST 10 PERIODS of each run
# (iwbcurv.py: `if t / T > g.nper - 10`). So a 20-period run averages periods
# 10-20, a 40-period run 30-40. Every validation number so far is a 20-period
# run.
#
# instab2 section B: ab = 20 on 896x351 moved from 51.72 deg at 20 periods to
# 53.63 at 40 (solver argmax), while max/mean rose 14.7 -> 16.3. Two readings:
#
#   SPIN-UP  -- the field has not reached its periodic state by periods 10-20.
#               The angle settles as nper grows and max|u| plateaus.
#   GROWTH   -- a slow version of the ab = 40 instability. max|u| keeps rising
#               and the angle keeps drifting, eventually diverging.
#
# Physical time is what matters for spin-up, so this runs on the cheap 512x201
# grid. The 20-period fields already exist (refit-npz for ab = 20, ab5-npz for
# ab = 5) and are copied in so the series is complete.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "time-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
$seed = @{ "20" = "refit-npz\f0.8_512x201.npz"; "5" = "ab5-npz\f0.8_512x201.npz" }
foreach ($ab in "20","5") {
  $dst = Join-Path $Out ("ab{0}_p20.npz" -f $ab)
  if (-not (Test-Path $dst) -and (Test-Path $seed[$ab])) { Copy-Item $seed[$ab] $dst }
  $lengths = if ($ab -eq "20") { @("40","60","80") } else { @("40","80") }
  foreach ($np in $lengths) {
    $npz = Join-Path $Out ("ab{0}_p{1}.npz" -f $ab, $np)
    if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
    Write-Host "-- ab=$ab  512x201  nper=$np"
    $t0 = Get-Date
    & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
        --Lx 6000 --nx 512 --nz 201 --ratio 0.8 --ab $ab `
        --nper $np --spp 120 --lsl 0.20 --taus 37.40 `
        --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
      Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
      Select-String -Pattern 't/T=|BEAM ANGLE|max/mean|DIVERGED|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
  }
}
Write-Host ""
Write-Host "now: python time_compare.py" -ForegroundColor Cyan
