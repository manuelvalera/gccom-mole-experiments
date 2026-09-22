# ab_test.ps1 -- does the -1.0 deg beam-angle limit depend on ridge height?
#
#   .\ab_test.ps1
#   python ab_compare.py 10=ab10-npz 20=refit-npz 40=ab40-npz
#
# WHY
#
# The limit survives fixing the crest (-1.10 vs -1.08), so crest sampling is
# not it. One candidate not yet tested is the DIAGNOSTIC: the tracker runs on
# rms of TOTAL speed, which includes the barotropic tide. Beam and tide
# oscillate at the same frequency with a phase difference, so rms^2 carries a
# cross term that varies across the beam and can shift its apparent centre --
# a shift that does not refine away, i.e. it lands in the limit.
#
# The solver is linear, so the beam amplitude scales with ridge height ab while
# the barotropic tide does not. In a linear beam the angle cannot depend on ab.
#
#   limit moves with ab   -> the tide is contaminating the diagnostic;
#                            the fix is to track baroclinic rms instead
#   limit independent     -> the offset is intrinsic to the computed beam
#
# ab = 10 and 40 bracket the existing ab = 20 runs in refit-npz. At ab = 40 the
# steepest ridge slope is 40/30 * exp(-1/2) = 0.81 against a beam slope of 1.33
# at w/N = 0.8, so the ridge stays subcritical and the physics is unchanged.
# Same three grids at every ab so the comparison is like for like.

param([string[]]$Heights = @("10","40"))

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
foreach ($ab in $Heights) {
  $Out = Join-Path (Get-Location) ("ab{0}-npz" -f $ab)
  if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
  foreach ($g in @("512","201"), @("640","251"), @("896","351")) {
    $npz = Join-Path $Out ("f0.8_{0}x{1}.npz" -f $g[0], $g[1])
    if (Test-Path $npz) { Write-Host "   have $npz, skipping"; continue }
    Write-Host "-- ab=$ab  $($g[0])x$($g[1])"
    $t0 = Get-Date
    & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
        --Lx 6000 --nx $g[0] --nz $g[1] --ratio 0.8 --ab $ab `
        --nper 20 --spp 120 --lsl 0.20 --taus 37.40 `
        --win 75 562 --alpha 1e-6 --forcegrid --out $npz 2>&1 |
      Tee-Object -FilePath ($npz -replace '\.npz$','.log') |
      Select-String -Pattern 'BEAM ANGLE|rms residual|max/mean|DIVERGED|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
  }
}
Write-Host ""
Write-Host "now: python ab_compare.py 10=ab10-npz 20=refit-npz 40=ab40-npz" -ForegroundColor Cyan
