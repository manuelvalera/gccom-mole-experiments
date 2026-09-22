# reproducer.ps1 -- measure the instability directly, and find a cheap case
# that reproduces it.
#
#   copy iwbcurv.py .   (the new one: --noise, --probe)
#   .\reproducer.ps1                    (~25 min, PARDISO at 4 threads)
#
# ESTABLISHED
#
# ab = 20 on 896x351, unramped: max|u| roughly doubles every 16 periods and the
# 80-period run is ruined. The same run at 512x201 is stable to 80 periods. At
# ab = 40 the growth is faster and the time step makes no difference, so it is
# a property of the discretised equations near the ridge.
#
# FREE RUNS
#
# --u0 0 --noise 1e-3: no tide, a small random initial field, same grid, same
# sponge. With nothing forcing it, any sustained growth is an instability, and
# --probe fits its e-folding time instead of leaving it to be read off by eye.
#
# THE REPRODUCER
#
# The instability appears with fine resolution AROUND THE RIDGE. A domain a
# quarter as wide (Lx 1500 m, nx 224) keeps exactly the same dx near the ridge
# (1500/224 = 6000/896) and the same vertical grid, at a quarter of the cost.
# If its growth rate matches the full domain, every further experiment -- the
# eigenvalue analysis, operator swaps, stretching, order -- can run on it in
# minutes instead of 20.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
if (-not $env:MKL_NUM_THREADS) { $env:MKL_NUM_THREADS = "4" }
$Out = "repro-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }

function Free([string]$tag, [string[]]$grid) {
  $npz = Join-Path $Out "$tag.npz"
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; return }
  Write-Host "-- $tag   $($grid -join ' ')" -ForegroundColor Cyan
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 --ratio 0.8 --spp 120 `
      --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --forcegrid `
      --solver pardiso --u0 0 --noise 1e-3 --probe --nper 40 --out $npz @grid 2>&1 |
    Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
    Select-String -Pattern 'random initial|probe t/T=\s+(10|20|30|40)\.0|GROWTH|DIVERGED|ms/step|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Free "full_896_ab20"    @("--Lx","6000","--nx","896","--nz","351","--ab","20")
Free "narrow_224_ab20"  @("--Lx","1500","--nx","224","--nz","351","--ab","20")
Free "full_512_ab20"    @("--Lx","6000","--nx","512","--nz","201","--ab","20")
Free "narrow_224_ab40"  @("--Lx","1500","--nx","224","--nz","351","--ab","40")

Write-Host ""
Write-Host "summary:" -ForegroundColor Cyan
Select-String -Path "$Out\*.log" -Pattern "GROWTH" | ForEach-Object {
  Write-Host ("   {0,-22} {1}" -f (Split-Path $_.Path -Leaf), $_.Line.Trim()) }
Write-Host ""
Write-Host "then: python view.py $Out\full_896_ab20.npz $Out\narrow_224_ab20.npz --log --xlim 700 -o repro.png" -ForegroundColor Cyan
