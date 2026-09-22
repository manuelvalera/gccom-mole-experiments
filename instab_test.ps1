# instab_test.ps1 -- is the shallow-beam family at ab = 40 a GROWING mode?
#
#   .\instab_test.ps1
#   python view.py instab-npz\p10.npz ab40-npz\f0.8_896x351.npz instab-npz\p40.npz --rms --log --xlim 1600 -o instab.png
#
# WHAT THE IMAGE SHOWED
#
# At ab = 40 m on the fine grids a second family of beams radiates from the
# ridge at ~22 deg, alongside the forced beam at 53 deg. At 640x251 the tracker
# hops between the two; at 896x351 the shallow family dominates (max/mean 85,
# against 11-19 for every clean run).
#
# The solver is linear and time-invariant: its only driver is the tide at
# w/N = 0.8, and such a system can only respond at the forcing frequency. A
# beam at 22 deg belongs to a wave at about sin(22 deg) = 0.37 N -- a frequency
# nobody is forcing -- so it is a FREE mode of the discrete system. For it to
# overwhelm the forced beam in 20 periods it must be growing.
#
# THE TEST
#
# Same grid, same ridge, three run lengths. If it is an instability, the shallow
# family and the field's max/mean grow with run length (and a long enough run
# may trip the divergence detector). If it is a startup transient that simply
# decays slowly, it FADES with run length. The 20-period run already exists in
# ab40-npz, so only 10 and 40 are needed.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "instab-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
foreach ($np in "10","40") {
  $npz = Join-Path $Out ("p{0}.npz" -f $np)
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
  Write-Host "-- ab=40  896x351  nper=$np"
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx 896 --nz 351 --ratio 0.8 --ab 40 `
      --nper $np --spp 120 --lsl 0.20 --taus 37.40 `
      --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
    Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
    Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|UNSTABLE|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}
Write-Host ""
Write-Host "compare max/mean: 10 periods (above), 20 periods = 84.8, 40 periods (above)" -ForegroundColor Cyan
Write-Host "then: python view.py instab-npz\p10.npz ab40-npz\f0.8_896x351.npz instab-npz\p40.npz --rms --log --xlim 1600 -o instab.png" -ForegroundColor Cyan
