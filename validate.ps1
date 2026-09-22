# validate.ps1 -- the full re-measurement, with every fix in place.
#
#   copy iwbcurv.py .
#   .\validate.ps1                                  (~15 min on the GPU)
#   python seiche_report.py val-logs
#   python centroid_track.py val-npz
#
# WHY NOW
#
# Every earlier beam number was taken before three things were understood:
#   - the start-up transient (fixed by --ramp 3; settled by ~15 periods)
#   - the buoyancy coupling that leaks energy over a sloping bed
#     (fixed by --buoy energy; the ab=40 replica stopped diverging)
#   - alpha (fixed; default now 1e-6)
# and the GPU loop has made a 20-period run at 896x351 take 20 seconds instead
# of six minutes, so the whole series can be redone at settled run lengths.
#
# SECTIONS
#   A  stability: free runs on the narrow replica, with and without the fix,
#      long enough for the old instability to have shown itself
#   B  seiche convergence, to reconfirm second order with the final code
#   C  the beam series at w/N = 0.8 and 0.6, six grids, ramped, settled,
#      averaged over periods 15-25
#
# Everything runs on the GPU. Drop "--device gpu --solver cudss" to use the CPU.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$gpu = @("--solver", "cudss", "--device", "gpu")
$npz = "val-npz"
$logs = "val-logs"
foreach ($d in $npz, $logs) {
    if (-not (Test-Path $d)) {
        New-Item -ItemType Directory -Path $d | Out-Null
    }
}
$script:N = 0

function Run([string]$tag, [string[]]$a) {
    $script:N++
    $log = Join-Path $logs "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $t0 = Get-Date
    # Tee-Object appends in PowerShell's default (UTF-16) encoding, which mixes
    # with the UTF-8 first line and makes the log unreadable to the Python
    # parsers. Capture, then write the whole log as UTF-8 in one go.
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @gpu @a 2>&1
    $header = "CMD iwbcurv.py $($a -join ' ')"
    ,$header + $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'GROWTH|DIVERGED|BEAM ANGLE|relative error|ms/step|max/mean|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

$narrow = @("--beat", "9999", "--Lx", "1500", "--nx", "224", "--nz", "351", "--ratio", "0.8",
            "--spp", "120", "--lsl", "0.20", "--taus", "37.40", "--win", "75", "562",
            "--alpha", "1e-6", "--forcegrid", "--u0", "0", "--noise", "1e-3", "--probe")

Write-Host "=== A. stability of the buoyancy fix (free runs, narrow replica) ===" -ForegroundColor Cyan
Run "A_ab40_energy_p100" ($narrow + @("--ab", "40", "--nper", "100", "--buoy", "energy"))
Run "A_ab40_logical_p25" ($narrow + @("--ab", "40", "--nper", "25", "--buoy", "logical"))
Run "A_ab40_N0_p100"     ($narrow + @("--ab", "40", "--nper", "100", "--N", "0"))
Run "A_ab20_energy_p150" ($narrow + @("--ab", "20", "--nper", "150", "--buoy", "energy"))

Write-Host "=== B. seiche convergence with the final code ===" -ForegroundColor Cyan
foreach ($g in @("128", "51"), @("192", "76"), @("256", "101"), @("384", "151")) {
    Run "B_seiche_$($g[0])x$($g[1])" @("--beat", "9999", "--Lx", "6000", "--nx", $g[0], "--nz", $g[1],
        "--ab", "0", "--bulge", "0", "--seiche", "4", "2", "--nper", "6", "--spp", "400",
        "--alpha", "1e-6", "--buoy", "energy")
}

Write-Host "=== C. beam series, ramped and settled, buoyancy fix on ===" -ForegroundColor Cyan
foreach ($cfg in @("0.8", "37.40", "75", "562"), @("0.6", "49.87", "133", "1000")) {
    foreach ($g in @("256", "101"), @("384", "151"), @("512", "201"),
                   @("640", "251"), @("768", "301"), @("896", "351")) {
        $tag = "C_f$($cfg[0])_$($g[0])x$($g[1])"
        Run $tag @("--beat", "9999", "--Lx", "6000", "--nx", $g[0], "--nz", $g[1],
            "--ratio", $cfg[0], "--nper", "25", "--spp", "120", "--lsl", "0.20",
            "--taus", $cfg[1], "--win", $cfg[2], $cfg[3], "--alpha", "1e-6",
            "--ramp", "3", "--buoy", "energy", "--forcegrid",
            "--out", (Join-Path $npz "f$($cfg[0])_$($g[0])x$($g[1]).npz"))
    }
}

Write-Host ""
Write-Host "now:" -ForegroundColor Cyan
Write-Host "   python seiche_report.py $logs" -ForegroundColor Cyan
Write-Host "   python centroid_track.py $npz" -ForegroundColor Cyan
Write-Host "   Select-String -Path $logs\A_*.log -Pattern 'GROWTH  ridge'" -ForegroundColor Cyan
