# ablate.ps1 -- which part of the scheme is unstable?
#
#   copy iwbcurv.py .     (new: near-ridge probe, growth report on divergence)
#   .\ablate.ps1          (~30 min)
#
# THE REPRODUCER
#
# Narrow domain (Lx 1500 m, nx 224, nz 351) -- same resolution around the ridge
# as 896x351 at a quarter of the cost. Free run: no tide, small random initial
# field. At ab = 40 it DIVERGED at t/T = 20.2 in 2 minutes, matching the
# full-domain run. Each line below changes ONE thing from that baseline.
#
# HOW TO READ IT
#
#   diverges like baseline  -> that piece is not the cause
#   stops diverging         -> that piece is the cause, or part of it
#
#   --N 0          no stratification. Kills buoyancy completely; what is left is
#                  the projection, the bed constraint and the sponge. If it
#                  still grows, the instability is in the pressure/velocity/bed
#                  coupling. If it stops, it is in the buoyancy coupling -- the
#                  logical-space interpolators Ic_z and Idf_w on skewed cells.
#   --alpha 1e-4   the lateral Robin coefficient.
#   --order 4      the interior operators.
#   --lsl 0.10     half the sponge width.
#   --bt/--bb 0.1  almost no horizontal stretching, so far less skew near the
#                  ridge -- the grid-geometry suspect.
#   nx 128 nz 201  the 512x201 resolution around the ridge: does ab = 40 need
#                  the fine grid to go unstable, as in the full domain?
#
# (The lagged pressure gradient gp is NOT on the list: it is a pure gradient,
# which an exact projection removes completely, so it cannot affect velocity.)
#
# The last run is 150 periods at ab = 20: long enough for a mode with an
# e-folding time of ~23 periods to rise out of the noise, which settles whether
# ab = 20 carries the same instability.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
if (-not $env:MKL_NUM_THREADS) { $env:MKL_NUM_THREADS = "4" }
$Out = "ablate-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
$narrow = @("--Lx","1500","--nx","224","--nz","351")

function Free([string]$tag, [string[]]$a) {
  $npz = Join-Path $Out "$tag.npz"; $log = $npz -replace '\.npz$','.log'
  if (Test-Path $log) { Write-Host "   have $log, skipping"; return }
  Write-Host "-- $tag   $($a -join ' ')" -ForegroundColor Cyan
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 --ratio 0.8 --spp 120 `
      --taus 37.40 --win 75 562 --alpha 1e-6 --forcegrid --solver pardiso `
      --u0 0 --noise 1e-3 --probe --out $npz @a 2>&1 |
    Tee-Object -FilePath $log |
    Select-String -Pattern 'GROWTH|DIVERGED|\*\*\*' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

$b = $narrow + @("--ab","40","--lsl","0.20","--nper","25")
Free "base_ab40"        $b
Free "N0_ab40"          ($b + @("--N","0"))
Free "alpha4_ab40"      ($b + @("--alpha","1e-4"))
Free "order4_ab40"      ($b + @("--order","4"))
Free "lsl10_ab40"       ($narrow + @("--ab","40","--lsl","0.10","--nper","25"))
Free "nostretch_ab40"   ($b + @("--bt","0.1","--bb","0.1"))
Free "coarse_ab40"      @("--Lx","1500","--nx","128","--nz","201","--ab","40","--lsl","0.20","--nper","25")
Free "long_ab20"        ($narrow + @("--ab","20","--lsl","0.20","--nper","150"))

Write-Host ""
Write-Host "summary (ridge norm is the one that matters):" -ForegroundColor Cyan
Get-ChildItem "$Out\*.log" | ForEach-Object {
  $f = $_.Name
  $d = Select-String -Path $_.FullName -Pattern 'DIVERGED at t/T=([0-9.]+)' | Select-Object -First 1
  $r = Select-String -Path $_.FullName -Pattern 'GROWTH  ridge' | Select-Object -Last 1
  $dtxt = if ($d) { "DIVERGED t/T=" + $d.Matches[0].Groups[1].Value } else { "ran to end" }
  $rtxt = if ($r) { ($r.Line -replace '^\s*GROWTH\s+ridge \|\|w\|\|\s+','') } else { "" }
  Write-Host ("   {0,-20} {1,-22} {2}" -f $f, $dtxt, $rtxt)
}
