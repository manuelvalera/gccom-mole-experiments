# reflect_test.ps1 -- is the far-field bending the approaching reflection?
#
#   copy iwbcurv.py .
#   .\reflect_test.ps1              (~6 min on the GPU)
#   python curvature.py reflect-npz
#
# WHAT IS ESTABLISHED
#
# With the ramp, the buoyancy fix and alpha = 1e-6, the near field (0.10-0.40
# of a bounce) converges to theory at both frequencies (+0.12 deg at w/N = 0.8,
# +0.37 at 0.6). But the measured angle still drifts as the fit window moves
# outward -- 2.18 deg of spread at 0.8 -- and the outermost band is the most
# negative on every grid. At 0.70-0.85 of a bounce the beam is within about a
# beam width of the surface, so it overlaps its own reflection.
#
# THE TEST
#
# Move the reflection and keep everything else. Deepening the domain lengthens
# the bounce (D0 / tan(theta)); scaling Lx with it keeps the domain the same
# number of bounces wide and the sponge the same distance away, which is what
# spoiled the earlier attempt at this. The ridge (ab, Lb) and the resolution
# (dx, dz) are held fixed, so the near field is physically identical.
#
#   bands line up by BOUNCE FRACTION      -> the reflection sets the bending
#   bands line up by ABSOLUTE DISTANCE    -> it comes from the ridge instead
#
# Run at a moderate resolution (dx ~ 13 m, dz ~ 5.7 m) so the deep case stays
# well inside GPU memory; the comparison is between the two depths, not against
# the converged answer.

if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
$Out = "reflect-npz"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out | Out-Null
}

function Go([string]$tag, [string[]]$a) {
    $npz = Join-Path $Out "$tag.npz"
    if (Test-Path $npz) {
        Write-Host "   have $npz, skipping"
        return
    }
    Write-Host "-- $tag" -ForegroundColor Cyan
    $t0 = Get-Date
    $out = & python iwbcurv.py --mole $env:MOLE_SRC --solver cudss --device gpu `
        --beat 9999 --ratio 0.8 --nper 25 --spp 120 --lsl 0.20 --taus 37.40 `
        --alpha 1e-6 --ramp 3 --buoy energy --forcegrid --ab 20 --lb 30 `
        --out $npz @a 2>&1
    $out | Out-File -FilePath ($npz -replace '\.npz$', '.log') -Encoding utf8
    $out | Select-String -Pattern 'BEAM ANGLE|max/mean|ms/step|resident|DIVERGED|MemoryError|\*\*\*' |
      ForEach-Object { Write-Host "   $($_.Line.TrimEnd())" }
    Write-Host ("   {0:n1} min" -f ((Get-Date) - $t0).TotalMinutes)
}

# same dx and dz in both: Lx/nx = 6000/448 = 9000/672, D0/nz = 1000/176 = 1500/264
Go "D1000" @("--Lx", "6000", "--nx", "448", "--nz", "177", "--D0", "1000", "--win", "75", "562")
Go "D1500" @("--Lx", "9000", "--nx", "672", "--nz", "265", "--D0", "1500", "--win", "112", "843")

Write-Host ""
Write-Host "curvature.py reads D0 from each file, so its bands are already in" -ForegroundColor Cyan
Write-Host "bounce fractions: compare the two tables row by row." -ForegroundColor Cyan
Write-Host "   python curvature.py $Out" -ForegroundColor Cyan
