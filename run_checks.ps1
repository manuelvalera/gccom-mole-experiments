# run_checks.ps1 -- convergence / confound checks for the 3-D cone at omega/N = 0.6
#
#   .\run_checks.ps1
#   .\run_checks.ps1 -Mole "C:/path/to/mole/src/matlab_octave" -Threads 8
#
# Writes a full log plus a one-line-per-run summary.  Safe to leave unattended:
# each run is independent, a crash or out-of-memory kill only loses that run.

param(
    [string]$Mole    = "C:/Users/dap21/gccom-mole-experiments-2026/mole/src/matlab_octave",
    [string]$Script  = ".\cone3d.py",
    [int]   $Threads = 8,
    [double]$Ratio   = 0.6
)

$env:MKL_NUM_THREADS  = "$Threads"
$env:OMP_NUM_THREADS  = "$Threads"
$env:CONE3D_THREADS   = "$Threads"
$env:MKL_DYNAMIC      = "FALSE"

$stamp   = Get-Date -Format "yyyyMMdd_HHmmss"
$log     = "cone_checks_$stamp.log"
$summary = "cone_checks_$stamp.txt"

# label, a, sig, spg, nper
$runs = @(
    @{n="baseline            "; a=64; sig=0.070; spg=0.05; nper=6 },
    @{n="sponge pulled back  "; a=64; sig=0.070; spg=0.03; nper=6 },
    @{n="more cells on source"; a=80; sig=0.070; spg=0.05; nper=6 },
    @{n="more travel         "; a=80; sig=0.056; spg=0.05; nper=6 },
    @{n="both levers         "; a=96; sig=0.056; spg=0.03; nper=6 },
    @{n="longer run          "; a=64; sig=0.070; spg=0.05; nper=12}
)

$theory = [math]::Round([math]::Asin($Ratio) * 180 / [math]::PI, 2)

"cone3d convergence checks   omega/N=$Ratio   theory $theory deg   threads=$Threads" |
    Tee-Object -FilePath $summary
"started $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')" | Tee-Object -FilePath $summary -Append
("-" * 92) | Tee-Object -FilePath $summary -Append
("{0,-20} {1,4} {2,6} {3,6} {4,5} {5,10} {6,9} {7,8}" -f `
    "run","a","sig","spg","nper","measured","error","minutes") |
    Tee-Object -FilePath $summary -Append

foreach ($r in $runs) {
    $label = $r.n
    $args  = @("cart","--mole",$Mole,"--a",$r.a,"--nper",$r.nper,
               "--sig",$r.sig,"--spg",$r.spg,"--dtfac","2.0","--spp","60",
               "--ratio",$Ratio)

    "`n`n=================================================================" | Add-Content $log
    "RUN: $label   a=$($r.a) sig=$($r.sig) spg=$($r.spg) nper=$($r.nper)"   | Add-Content $log
    "started $(Get-Date -Format 'HH:mm:ss')"                                 | Add-Content $log
    "=================================================================`n"   | Add-Content $log

    Write-Host "[$(Get-Date -Format 'HH:mm:ss')] running: $label (a=$($r.a)) ..." -NoNewline
    $t0  = Get-Date
    $out = & python $Script @args 2>&1
    $out | Add-Content $log
    $mins = [math]::Round(((Get-Date) - $t0).TotalMinutes, 1)

    $line = $out | Select-String "CONE HALF-ANGLE" | Select-Object -Last 1
    if ($line) {
        $m = [regex]::Match($line.ToString(),
             'measured\s+([-\d.]+).*?error\s+([-+\d.]+)\s+deg\s+\(([-+\d.]+%)\)')
        $meas = $m.Groups[1].Value
        $err  = "$($m.Groups[2].Value) ($($m.Groups[3].Value))"
    } else {
        $meas = "FAILED"
        $err  = "see log"
    }
    Write-Host "  -> $meas deg   ($mins min)"

    ("{0,-20} {1,4} {2,6} {3,6} {4,5} {5,10} {6,9} {7,8}" -f `
        $label.Trim(), $r.a, $r.sig, $r.spg, $r.nper, $meas, $err, $mins) |
        Tee-Object -FilePath $summary -Append
}

("-" * 92)                                              | Tee-Object -FilePath $summary -Append
"finished $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"    | Tee-Object -FilePath $summary -Append
"theory = $theory deg.  If every run sits near 45, the offset is real."  |
    Tee-Object -FilePath $summary -Append
"full output: $log"                                     | Tee-Object -FilePath $summary -Append

Write-Host "`nSummary written to $summary"
Write-Host "Full log written to $log"
Get-Content $summary
