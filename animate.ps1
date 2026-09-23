# animate.ps1 -- watch a run develop, from the first step to the formed beams.
#
#   copy iwbcurv.py .
#   copy animate.py .
#   .\animate.ps1
#
# Two animations: the Garcia benchmark ridge, where the beams appear out of the
# start-up transient over about three periods, and the Monterey shelf transect.
#
# The solver writes one snapshot every 1/framerate of a period into a
# directory; animate.py renders them on the real curvilinear grid with the
# colour scale fixed across the run, so brightening is the field growing rather
# than the scale moving. The depth-mean flow is removed by default -- the
# barotropic tide is several times larger than the beams and hides them.
#
# Disk: a frame is about 0.1 MB at 256x101 and 1.5 MB at 896x351. 12 per period
# for 6 periods is 72 frames, so keep --framerate modest on the big grids.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }

Write-Host "-- benchmark ridge, 6 periods" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
    --Lx 6000 --nx 512 --nz 201 --ratio 0.8 --nper 6 --spp 120 `
    --lsl 0.20 --taus 37.40 --win 75 562 --alpha 1e-6 --ramp 3 --buoy energy `
    --forcegrid --solver cudss --device gpu `
    --frames frames_ridge --framerate 12 |
  Select-String -Pattern 'frames ->|ms/step|BEAM ANGLE' |
  ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
python animate.py frames_ridge -o ridge.gif --field speed --fps 12

Write-Host "-- Monterey shelf transect (trimmed), 6 periods" -ForegroundColor Cyan
& python iwbcurv.py --mole $env:MOLE_SRC --beat 9999 `
    --grids grids --gridname mryshelf_t20 --nx 768 --nz 151 `
    --nprofile mry_N.txt --omega 1.405e-4 --taus 1490 --lsl 0.15 `
    --btmode transport --ramp 3 --buoy energy --forcegrid `
    --solver cudss --device gpu --nper 6 --spp 600 `
    --frames frames_mry --framerate 8 |
  Select-String -Pattern 'frames ->|ms/step|max/mean|DIVERGED' |
  ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
python animate.py frames_mry -o monterey.gif --field speed --fps 10

Write-Host ""
Write-Host "also try: python animate.py frames_ridge -o ridge_w.gif --field w" -ForegroundColor Cyan
Write-Host "  (signed vertical velocity, which shows the phase moving through the beam)" -ForegroundColor Cyan
