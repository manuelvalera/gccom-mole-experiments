# fix_test.ps1 -- does the energy-consistent buoyancy coupling remove the
# instability for good, and what does it do to the beam?
#
#   copy iwbcurv.py .     (new: --buoy energy)
#   .\fix_test.ps1        (~45 min, PARDISO at 4 threads)
#
# THE DIAGNOSIS
#
# Ablation: with --N 0 the instability disappears, so it lives in the w <-> b
# exchange. MOLE's interpolator pair is one-way at the bed: the bed face takes
# its buoyancy from the boundary point, while the first interior cell takes
# half of its forcing from the bed face. On a flat bottom w_bed = 0 and nothing
# happens (the seiche never saw it); over a slope w_bed = slope * u and energy
# leaks in at the ridge. --buoy energy rebuilds the w-forcing as the exact
# adjoint of Idf_w in the physical-area inner product (row sums exactly 1,
# adjoint to ~1e-15).
#
# Verified in a single-core test: the ab = 40 replica that diverged at
# t/T = 20.2 ran to 25 periods, with the near-ridge growth rate down from
# +0.185 to +0.006 per period. These runs check whether that residual is real.
#
#   A  ab = 40 replica, 100 periods       residual growth: does it keep going?
#   B  ab = 20 replica, 150 periods       the original coupling DIVERGED at 79
#   C  ab = 40 replica, --N 0, 100 periods  the no-buoyancy baseline to compare A with
#   D  896x351, ab = 20, abrupt start, 80 periods, forced
#                                         the run that blew up (x2 per ~16 periods);
#                                         should now settle, and gives a beam angle

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
if (-not $env:MKL_NUM_THREADS) { $env:MKL_NUM_THREADS = "4" }
$Out = "fix-npz"
if (-not (Test-Path $Out)) { New-Item -ItemType Directory -Path $Out | Out-Null }
$common = @("--beat","9999","--ratio","0.8","--spp","120","--taus","37.40","--win","75","562",
            "--alpha","1e-6","--forcegrid","--solver","pardiso","--lsl","0.20")
$narrow = @("--Lx","1500","--nx","224","--nz","351")
$free   = @("--u0","0","--noise","1e-3","--probe")

function Go([string]$tag, [string[]]$a) {
  $npz = Join-Path $Out "$tag.npz"; $log = $npz -replace '\.npz$','.log'
  if (Test-Path $log) { Write-Host "   have $log, skipping"; return }
  Write-Host "-- $tag" -ForegroundColor Cyan
  $t0 = Get-Date
  & python iwbcurv.py --mole $env:MOLE_SRC @common --out $npz @a 2>&1 |
    Tee-Object -FilePath $log |
    Select-String -Pattern 'buoyancy coupling|GROWTH|DIVERGED|t/T= +(16|32|48|64|80)\.0 +max|BEAM ANGLE|max/mean' |
    ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
  Write-Host ("   took {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Go "A_energy_ab40_p100"  ($narrow + $free + @("--ab","40","--nper","100","--buoy","energy"))
Go "B_energy_ab20_p150"  ($narrow + $free + @("--ab","20","--nper","150","--buoy","energy"))
Go "C_N0_ab40_p100"      ($narrow + $free + @("--ab","40","--nper","100","--N","0"))
Go "D_energy_896_ab20_p80" @("--Lx","6000","--nx","896","--nz","351","--ab","20","--nper","80","--buoy","energy")

Write-Host ""
Write-Host "A vs C: if their ridge rates match, the fix is complete." -ForegroundColor Cyan
Write-Host "B: ran to 150 = the ab = 20 instability is gone (it diverged at 79 before)." -ForegroundColor Cyan
Write-Host "D: max|u| flat to 80 periods = fixed where it mattered; note the beam angle." -ForegroundColor Cyan
