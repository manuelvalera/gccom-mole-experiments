# mry_sponge.ps1 -- fix the absorption, then find out what the transect actually does.
#
#   .\mry_sponge.ps1            (~12 min on the GPU)
#
# WHAT IS ESTABLISHED
#
# The steady amplitude depends almost entirely on the sponge: ||u|| settles at
# 13.4, 4.15, 5.08, 1.65 and 0.80 across five sponge/domain combinations. So
# the earlier numbers measured the sponge, not the shelf. A wide, slowly
# relaxed sponge (lsl 0.30, taus 450) reaches a flat steady state by t/T = 10;
# a narrow fast one is worse than a narrow slow one, because an abrupt damping
# profile reflects instead of absorbing.
#
# The untrimmed transect is also steady with that sponge, which matters: the
# shallow end is where a shoaling internal tide would steepen and break, so
# keeping it is worth something if it is numerically sound.
#
#   A  sponge independence. Vary width and rate about the chosen setting. If
#      the answer moves by more than a few percent, absorption is still
#      setting it and no amplitude is trustworthy.
#   B  is the shallow end really safe now? Free runs on the untrimmed domain
#      with the good sponge, at two resolutions -- these grew before.
#   C  grid convergence of the forced solution, which is the actual validation
#      question and has not been answerable until now.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "mry-sponge"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "9999", "--grids", "grids", "--gridname", "mryshelf",
            "--Lx", "10110", "--D0", "87.8", "--nprofile", "mry_N.txt",
            "--omega", "1.405e-4", "--forcegrid", "--solver", "cudss", "--device", "gpu",
            "--buoy", "energy", "--spp", "600", "--ulim", "5000",
            "--probe", "--probebox", "0", "20000", "500")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'probe t/T= +(10|16)\.0|GROWTH  domain|DIVERGED|max/mean|max\|u\|' |
      Select-Object -Last 5 | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. is the answer independent of the sponge? ===" -ForegroundColor Cyan
foreach ($s in @("0.30", "450"), @("0.40", "450"), @("0.30", "300"), @("0.40", "300")) {
    Go "indep_l$($s[0])_t$($s[1])" @("--lsl", $s[0], "--taus", $s[1], "--nx", "384", "--nz", "76",
                                     "--btmode", "transport", "--ramp", "3", "--nper", "16")
}

Write-Host "=== B. shallow end with proper absorption ===" -ForegroundColor Cyan
foreach ($g in @("384", "76"), @("768", "151")) {
    Go "free_$($g[0])" @("--lsl", "0.30", "--taus", "450", "--nx", $g[0], "--nz", $g[1],
                         "--u0", "0", "--noise", "1e-3", "--nper", "20")
}

Write-Host "=== C. grid convergence of the forced solution ===" -ForegroundColor Cyan
foreach ($g in @("384", "76"), @("768", "151"), @("1152", "226")) {
    Go "conv_$($g[0])" @("--lsl", "0.30", "--taus", "450", "--nx", $g[0], "--nz", $g[1],
                         "--btmode", "transport", "--ramp", "3", "--nper", "16",
                         "--out", (Join-Path $Out "conv_$($g[0]).npz"))
}

Write-Host ""
Write-Host "A: ||u|| at t/T=16 within a few percent across all four -> absorption is adequate." -ForegroundColor Cyan
Write-Host "B: ||u|| decaying at both resolutions -> the untrimmed domain is usable." -ForegroundColor Cyan
Write-Host "C: ||u|| agreeing across resolutions -> a converged solution, at last." -ForegroundColor Cyan
