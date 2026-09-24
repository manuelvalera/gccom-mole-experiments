# lock_test.ps1 -- lock release, against the published mimetic-GCCOM numbers.
#
#   copy iwbcurv.py .
#   .\lock_test.ps1              (~6 min on the GPU)
#
# The configuration is the one in reference2021/LockRelease3Dv7fullymimeticmodular.m:
# 0.8 m long, 0.0656 m deep, g' = 0.0224 m/s^2, so sqrt(g'H) = 0.0383 m/s and an
# energy-conserving front travels at half of that. Reported as
# Fr = u_front / sqrt(g' H / 2), where the inviscid value is 1/sqrt(2) = 0.7071
# and the 2021 mimetic GCCOM reported 0.705 at 401x6x101.
#
# NOTE ON alpha. It regularises a matrix whose entries scale as 1/dx^2, so it is
# scale-dependent, not a constant: 1e-6 is right for the 6 km benchmark and
# breaks the projection outright at 0.8 m, where ~1 gives residuals of 1e-14.
#
# A  resolution: does Fr approach the inviscid value as the grid refines?
# B  scheme: minmod against upwind2 measures how much the limiter's diffusion
#    costs -- the same question that decides whether a beam survives.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "lock-logs"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}
$common = @("--beat", "199", "--Lx", "0.8", "--D0", "0.0656", "--ab", "0", "--bulge", "0",
            "--N", "0", "--lock", "0.0224", "--tend", "10", "--alpha", "1",
            "--advect", "full", "--buoy", "energy", "--forcegrid",
            "--solver", "cudss", "--device", "gpu")

function Go([string]$tag, [string[]]$a) {
    $log = Join-Path $Out "$tag.log"
    if (Test-Path $log) {
        Write-Host "   have $log, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $out = & python iwbcurv.py --mole $env:MOLE_SRC @common @a 2>&1
    $out | Out-File -FilePath $log -Encoding utf8
    $out | Select-String -Pattern 'LOCK RELEASE|Fr =|static projection|DIVERGED' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
}

Write-Host "=== A. resolution, minmod ===" -ForegroundColor Cyan
Go "res_101"  @("--nx", "101", "--nz", "26", "--dt", "0.008", "--advscheme", "minmod")
Go "res_201"  @("--nx", "201", "--nz", "51", "--dt", "0.004", "--advscheme", "minmod")
Go "res_401"  @("--nx", "401", "--nz", "101", "--dt", "0.002", "--advscheme", "minmod")
Go "res_801"  @("--nx", "801", "--nz", "201", "--dt", "0.001", "--advscheme", "minmod")

Write-Host "=== B. scheme at 401x101 ===" -ForegroundColor Cyan
Go "sch_upwind1" @("--nx", "401", "--nz", "101", "--dt", "0.002", "--advscheme", "upwind1")
Go "sch_upwind2" @("--nx", "401", "--nz", "101", "--dt", "0.002", "--advscheme", "upwind2")

Write-Host ""
Write-Host "A: Fr should rise toward 0.707 as the grid refines." -ForegroundColor Cyan
Write-Host "B: upwind1 is the 2021 scheme -- the gap to upwind2 is what first-order" -ForegroundColor Cyan
Write-Host "   advection costs, which is what dissipated the beams." -ForegroundColor Cyan
