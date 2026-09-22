# instab2.ps1 -- what kind of instability is it, and are the "clean" runs
# actually stable?
#
#   .\instab2.ps1
#
# ESTABLISHED
#
# At ab = 40 m on 896x351 the field grows without bound: max/mean 37.7 at 10
# periods, 84.8 at 20, and the 40-period run DIVERGED at t/T = 20.33 with
# max|u| = 60 u0. A family of beams at ~22 deg -- a frequency nobody forces --
# is the unstable mode radiating.
#
# SECTION A: time step or operator?
#
# dt = T/spp is fixed by the forcing period, NOT by the grid. If the growth
# depends on dt, some stability limit tightens as cells shrink or skew near the
# ridge, and a smaller dt suppresses it -- a nuisance with a known cure. If the
# growth is the same at every dt, the SPATIAL operator itself has a growing
# eigenmode: no time step will fix it, and the defect is in how the curvilinear
# terms are discretised near steep topography. That is the more serious result,
# and the one that would matter for Monterey.
#
#   spp 60 / 120 / 240 at 20 periods: compare max/mean with the 84.8 at spp=120.
#
# SECTION B: are the clean runs clean?
#
# Every validation number so far was taken at 20 periods with ab <= 20 m, where
# max/mean sat at 11-19 and looked healthy. But a SLOW instability would look
# exactly like that at 20 periods. Doubling the run length at ab = 20 on the
# finest grid is the check: max/mean should stay near 14.7 and the beam angle
# should not move. If it does, the validation table needs a stability caveat.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = Join-Path (Get-Location) "instab-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }

function Go([string]$tag, [string[]]$extra) {
  $npz = Join-Path $Out ("{0}.npz" -f $tag)
  if (Test-Path $npz) { Write-Host "   have $npz, skipping"; return }
  Write-Host "-- $tag"
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
      --Lx 6000 --nx 896 --nz 351 --ratio 0.8 `
      --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --forcegrid --out $npz @extra 2>&1 |
    Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
    Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|UNSTABLE|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Write-Host "=== A. time step (ab = 40, 20 periods; spp=120 gave max/mean 84.8) ===" -ForegroundColor Cyan
Go "ab40_spp60"  @("--ab","40","--nper","20","--spp","60")
Go "ab40_spp240" @("--ab","40","--nper","20","--spp","240")

Write-Host "=== B. are the clean runs clean? (ab = 20, 40 periods; 20 periods gave max/mean 14.7) ===" -ForegroundColor Cyan
Go "ab20_p40" @("--ab","20","--nper","40","--spp","120")

Write-Host ""
Write-Host "A: max/mean falls at spp 240 and rises at 60 -> time-step limit (cure: smaller dt)" -ForegroundColor Cyan
Write-Host "   roughly the same at all three           -> growing mode of the spatial operator" -ForegroundColor Cyan
Write-Host "B: max/mean ~14.7 and angle unchanged       -> the validation runs are stable" -ForegroundColor Cyan
