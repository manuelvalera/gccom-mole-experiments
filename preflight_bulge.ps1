# preflight_bulge.ps1 -- confirm left.m/right.m actually honour --bulge.
#
# The solver's --seiche guard checks the FLAG, not the grid. If left.m and
# right.m do not read the BULGE global, --bulge 0 is accepted while the walls
# keep their 150 m sinusoidal bulge, the domain is not a rectangle, and the
# analytic frequency omega = N k / sqrt(k^2+p^2) does not apply to it.
# That is a ~3% offset masquerading as a numerical defect.
if (-not $env:MOLE_SRC) { Write-Error "set `$env:MOLE_SRC"; exit 1 }
foreach ($f in "grids\iwbridge\left.m", "grids\iwbridge\right.m") {
    if (Select-String -Path $f -Pattern "global BULGE" -Quiet) {
        Write-Host "   OK      $f reads BULGE"
    } else {
        Write-Host "   STALE   $f does NOT read BULGE -- replace it" -ForegroundColor Red
    }
}
Write-Host "   -- wander must differ between the two runs if BULGE is live --"
foreach ($b in "0.15", "0") {
    $o = & python iwbcurv.py --mole $env:MOLE_SRC --Lx 6000 --nx 64 --nz 33 `
                  --ab 0 --bulge $b --gridonly 2>&1 | Out-String
    if ($o -match "side-wall bulge[^\r\n]*") { $m = $Matches[0] }
    elseif ($o -match "wander[^\r\n]*") { $m = $Matches[0] }
    else { $m = "(no geometry line found)" }
    Write-Host ("   --bulge {0,-5} {1}" -f $b, $m)
}
