# gpu_bench.ps1 -- what the GPU loop is worth on a realistic run length.
#
#   copy iwbcurv.py .
#   .\gpu_bench.ps1            (~6 min)
#
# The earlier check used 3 periods, which overstates the non-solve cost: every
# step falls inside the final 10-period averaging window, and the five
# consistency checks are a fixed per-run cost spread over only 360 steps. This
# runs 20 periods, where the averaging covers half the steps, so the ms/step
# split is the one that applies to production runs.
#
# Also compares against the CPU loop on the same case, and checks the fields
# agree to round-off, as gpu_check did.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'time loop on the GPU|cudss factorise|splu|ms/step|BEAM ANGLE|DIVERGED|\*\*\*'
$case = @("--beat","9999","--Lx","6000","--nx","896","--nz","351","--ratio","0.8","--nper","20",
          "--spp","120","--lsl","0.20","--taus","37.40","--win","75","562","--alpha","1e-6",
          "--ramp","3","--buoy","energy","--forcegrid")

function Go([string]$tag, [string[]]$extra) {
    Write-Host "-- $tag" -ForegroundColor Cyan
    $t0 = Get-Date
    & python iwbcurv.py --mole $env:MOLE_SRC @case @extra 2>&1 |
      Select-String -Pattern $Pat |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   wall {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

Go "GPU loop, cuDSS" @("--solver", "cudss", "--device", "gpu", "--out", "bench_gpu.npz")
Go "CPU loop, cuDSS" @("--solver", "cudss", "--device", "cpu", "--out", "bench_cpu.npz")
Go "CPU loop, SuperLU" @("--solver", "superlu", "--device", "cpu")

& python -c @"
import numpy as np
a = np.load('bench_cpu.npz')
b = np.load('bench_gpu.npz')
for k in ('rms', 'snap_u', 'snap_w'):
    if k in a and k in b:
        x = a[k]
        y = b[k]
        rel = np.abs(x - y).max() / max(np.abs(x).max(), 1e-300)
        print(f'   {k:7s} GPU vs CPU: {rel:.1e}   ' + ('OK' if rel < 1e-9 else 'MISMATCH'))
"@
