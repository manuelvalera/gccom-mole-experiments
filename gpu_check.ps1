# gpu_check.ps1 -- is the all-GPU time loop correct, and how fast is it?
#
#   copy iwbcurv.py .
#   .\gpu_check.ps1            (~3 min)
#
# The CPU path is bit-for-bit unchanged by this port (checked). The GPU path
# is new, so it has to pass three things:
#
#   1. the seiche still reproduces -4.0844e-04
#   2. a short 896x351 beam run gives the 42.44 deg fingerprint
#   3. its whole saved field matches the CPU+cuDSS run of the same case to
#      round-off -- not bit-identical, because GPU sparse products add in a
#      different order, but differences should sit near 1e-12, not 1e-3
#
# Watch the card in a second window if you like:
#   nvidia-smi --query-gpu=temperature.gpu,power.draw,utilization.gpu --format=csv -l 5

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'time loop on the GPU|cudss factorise|static projection|consistency|ms/step|relative error|BEAM ANGLE|Error|Traceback|\*\*\*'
$beam = @("--beat","9999","--Lx","6000","--nx","896","--nz","351","--ratio","0.8","--nper","3",
          "--spp","120","--lsl","0.20","--taus","37.40","--win","75","562","--alpha","1e-6",
          "--forcegrid","--solver","cudss")

Write-Host "=== 1. seiche with the device-resident solver: must print -4.0844e-04 ===" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 --Lx 6000 --nx 256 --nz 101 --ab 0 --bulge 0 `
    --seiche 4 2 --nper 6 --spp 400 --alpha 1e-6 --solver cudss --device gpu 2>&1 |
  Select-String -Pattern $Pat | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }

Write-Host "=== 2. beam 896x351, 3 periods, GPU loop: must print 42.44 deg ===" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC @beam --device gpu --out gpu_check_gpu.npz 2>&1 |
  Select-String -Pattern $Pat | Select-Object -Last 5 | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }

Write-Host "=== 3. the same run with the CPU loop, for a field-by-field comparison ===" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC @beam --device cpu --out gpu_check_cpu.npz 2>&1 |
  Select-String -Pattern 'ms/step|BEAM ANGLE' | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }

& python -c @"
import numpy as np
a = np.load('gpu_check_cpu.npz'); b = np.load('gpu_check_gpu.npz')
for k in ('rms', 'snap_u', 'snap_w'):
    if k in a and k in b:
        x, y = a[k], b[k]
        rel = np.abs(x - y).max() / max(np.abs(x).max(), 1e-300)
        print(f'   {k:7s} max relative difference GPU vs CPU: {rel:.1e}   ' +
              ('OK' if rel < 1e-9 else 'MISMATCH -- do not use the GPU loop'))
"@
