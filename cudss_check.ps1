# cudss_check.ps1 -- is --solver cudss correct inside the model, and how fast?
#
#   copy iwbcurv.py .
#   .\cudss_check.ps1          (~3 min)
#
# Same bar PARDISO had to clear. The seiche must reproduce the SuperLU value
# exactly, and a short beam run must give the same angle as SuperLU did on the
# same grid (42.44 deg at 896x351, 3 periods -- not a physical angle, just a
# fingerprint of the whole field) while reporting the per-step split.
#
# To watch the GPU while it runs, in a second window:
#   nvidia-smi --query-gpu=temperature.gpu,power.draw,utilization.gpu --format=csv -l 5

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Pat = 'cudss factorise|splu|static projection|consistency|ms/step|relative error|BEAM ANGLE|Error|Traceback|\*\*\*'

Write-Host "=== 1. seiche 256x101, mode (4,2), alpha 1e-6: must print -4.0844e-04 ===" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 --Lx 6000 --nx 256 --nz 101 --ab 0 --bulge 0 `
    --seiche 4 2 --nper 6 --spp 400 --alpha 1e-6 --solver cudss 2>&1 |
  Select-String -Pattern $Pat | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }

Write-Host "=== 2. beam 896x351, 3 periods: must print 42.44 deg, like SuperLU ===" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 --Lx 6000 --nx 896 --nz 351 --ratio 0.8 `
    --nper 3 --spp 120 --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --forcegrid --solver cudss 2>&1 |
  Select-String -Pattern $Pat | Select-Object -Last 4 | ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }

Write-Host ""
Write-Host "Both match -> cudss is safe to use. Compare the ms/step line with the SuperLU" -ForegroundColor Cyan
Write-Host "run's 132 ms (109 solve + 24 everything else)." -ForegroundColor Cyan
