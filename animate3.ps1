# animate3.ps1 -- the three validated cases, with nonlinear advection on.
#
#   copy iwbcurv.py .
#   copy animate.py .
#   .\animate3.ps1              (~5 min on the GPU)
#
# Each case is the configuration it was verified with, run again with
# --advect full so the animations show the nonlinear solver rather than the
# linear one. Frames carry u, w and buoyancy; the renderer picks the field that
# makes each case legible.
#
#   lock release   buoyancy: the interface collapses into a gravity current and
#                  Kelvin-Helmholtz billows roll up along it. Fr ~ 0.69 against
#                  theory 0.7071 and the 2021 mimetic GCCOM's 0.705.
#   seiche         vertical velocity: a standing internal wave in a flat box,
#                  the case with an exact frequency. The seiche loop now carries
#                  the advective terms too. Run at --seicheamp 1e-3: the linear
#                  solver is amplitude-independent, but 1 m/s in a 1 km box is
#                  violently nonlinear and blows up. At 1e-3 the measured
#                  frequency matches the linear value to 0.04%, which is the
#                  linear limit of the nonlinear code -- a verification in its
#                  own right.
#   beam           baroclinic speed: the internal-wave beam over a ridge, the
#                  experiment the 2021 model could not sustain. 25 periods,
#                  steady, 53.55 deg against 53.13 theory.
#
# GIFs need nothing extra; for .mp4 ffmpeg must be on PATH.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }

Write-Host "-- lock release (0.8 m box, g' = 0.0224)" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
    --Lx 0.8 --D0 0.0656 --ab 0 --bulge 0 --nx 401 --nz 101 --N 0 `
    --lock 0.0224 --tend 10 --dt 0.002 --alpha 1 `
    --advect full --advscheme upwind2 --buoy energy --forcegrid `
    --solver cudss --device gpu --frames fr_lock --framerate 40 |
  Select-String -Pattern 'LOCK RELEASE|Fr =|frames ->' |
  ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
python animate.py fr_lock -o anim_lock.gif --field b --fps 12

Write-Host "-- seiche, mode (4,2), finite amplitude" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
    --Lx 6000 --nx 256 --nz 101 --ab 0 --bulge 0 --seiche 4 2 `
    --nper 3 --spp 400 --alpha 1e-6 --advect full --advscheme upwind2 `
    --buoy energy --seicheamp 1e-3 --frames fr_seiche --framerate 30 |
  Select-String -Pattern 'relative error|frames ->|DIVERGED' |
  ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
python animate.py fr_seiche -o anim_seiche.gif --field w --fps 12

Write-Host "-- internal-wave beam over a ridge, nonlinear" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
    --Lx 6000 --nx 512 --nz 201 --ratio 0.8 --nper 8 --spp 120 `
    --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --ramp 3 `
    --buoy energy --forcegrid --advect full --advscheme upwind2 `
    --solver cudss --device gpu --frames fr_beam --framerate 12 |
  Select-String -Pattern 'BEAM ANGLE|max/mean|frames ->' |
  ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
python animate.py fr_beam -o anim_beam.gif --field speed --fps 12

Write-Host ""
Write-Host "anim_lock.gif  anim_seiche.gif  anim_beam.gif" -ForegroundColor Cyan
Write-Host "the beam also reads well as signed w: python animate.py fr_beam -o anim_beam_w.gif --field w" -ForegroundColor Cyan
