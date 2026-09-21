# beam_alpha.ps1 -- redo the beam-angle convergence series at a floor-level
# alpha, and find out whether the first-order fit was alpha too.
#
#   .\beam_alpha.ps1
#   python compare_paper.py beam-alpha-logs
#
# ---------------------------------------------------------------------------
# WHY
#
# The seiche is second order once alpha is at the floor:
#
#     grid       alpha=1e-4    alpha=1e-6
#     128x51     +2.91e-03     -1.81e-03
#     192x76     +3.95e-03     -7.72e-04
#     256x101    +4.32e-03     -4.08e-04
#     384x151    +4.58e-03     -1.49e-04
#                GROWS, p=-0.5  SHRINKS, p=2.27
#
# At alpha=1e-4 it grew under refinement and fitted p=1 with a nonzero limit.
# At alpha=1e-6 it converges cleanly at second order to zero (e0 = +5.9e-05,
# residual 7.8e-08 on the p=2 fit).
#
# The beam-angle series that produced "first order, e0 = -1.63 deg" was run
# entirely at alpha=1e-4. At 512x201 alpha is only worth 0.07 deg, because the
# sponge also holds the walls there -- but that has never been measured at
# 256x101, where the Robin rows carry relatively more of the load.
#
# If the coarse end drops when alpha does, the first-order fit dissolves the
# same way the seiche's did, and the convergence claim in the paper changes
# from "first order to a nonzero limit" to "second order to zero".
#
# If the series is UNCHANGED, then the beam angle really does have a
# grid-independent offset that the seiche does not, and the cause is something
# the flat-box problem cannot see -- the ridge, the sponge, or the window.
# That is a different and still-open question, but at least a sharp one.
# ---------------------------------------------------------------------------

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'BEAM ANGLE|rms residual|max/mean|consistency|UNSTABLE|DIVERGED|\*\*\*'
$LogDir = Join-Path (Get-Location) "beam-alpha-logs"
if (-not (Test-Path $LogDir)) { New-Item -ItemType Directory -Path $LogDir | Out-Null }
$script:N = 0

function Run([string[]]$A) {
    foreach ($x in $A) {
        if ([string]::IsNullOrWhiteSpace($x)) {
            Write-Host "   !! blank argument" -ForegroundColor Red; return } }
    $cmd = @("iwbcurv.py","--mole",$env:MOLE_SRC,"--beat","9999") + $A
    $script:N++
    $log = Join-Path $LogDir ("b{0:d3}.log" -f $script:N)
    Write-Host ("   python " + ($cmd -join " ")) -ForegroundColor DarkGray
    "CMD " + ($cmd -join " ") | Out-File -FilePath $log -Encoding utf8
    & python @cmd 2>&1 | Out-File -FilePath $log -Encoding utf8 -Append
    $code = $LASTEXITCODE
    $hits = @(Select-String -Path $log -Pattern $Pat | ForEach-Object { $_.Line.TrimEnd() })
    if ($hits.Count -gt 0) { $hits | ForEach-Object { Write-Host "   $_" } }
    if ($hits.Count -eq 0 -or $code -ne 0) {
        Write-Host "   !! no result (exit=$code):" -ForegroundColor Red
        Get-Content $log -Tail 20 | ForEach-Object { Write-Host "   | $_" -ForegroundColor DarkYellow } }
}

# omega/N = 0.6: theory 36.87, window [133,1000], taus = T/30 = 49.87
# omega/N = 0.8: theory 53.13, window [75,562],  taus = T/30 = 37.40
# alpha 1e-4 is rerun at each grid so the comparison is like-for-like on this
# machine, not against numbers from an earlier session.
foreach ($cfg in @("0.6","49.87","133","1000"), @("0.8","37.40","75","562")) {
  foreach ($g in @("256","101"), @("384","151"), @("512","201"), @("640","251")) {
    foreach ($a in "1e-4","1e-6") {
      Write-Host "-- ratio=$($cfg[0])  $($g[0])x$($g[1])  alpha=$a"
      Run @("--Lx","6000","--nx",$g[0],"--nz",$g[1],"--ratio",$cfg[0],
            "--nper","20","--spp","120","--lsl","0.20","--taus",$cfg[1],
            "--win",$cfg[2],$cfg[3],"--alpha",$a,"--forcegrid")
    }
  }
}

Write-Host ""
Write-Host "logs in $LogDir" -ForegroundColor Cyan
Write-Host "now: python compare_paper.py beam-alpha-logs" -ForegroundColor Cyan
