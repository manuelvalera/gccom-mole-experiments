# animate_web.ps1 -- rebuild the three animations at a size sensible for a README.
#
#   .\animate3.ps1        (first: produces fr_lock, fr_seiche, fr_beam)
#   .\animate_web.ps1
#
# A full-resolution GIF of a 400-frame run runs to tens of MB, and every release
# of this repository is archived to Zenodo, so the web versions are cut down:
# every second or third frame, lower dpi, and the beam cropped to the part that
# matters. Aim for a couple of MB each.

$Out = "docs\animations"
if (-not (Test-Path $Out)) {
    New-Item -ItemType Directory -Path $Out -Force | Out-Null
}

python animate.py fr_lock   -o "$Out\lock_release.gif" --field b     --fps 12 --stride 2 --dpi 80
python animate.py fr_seiche -o "$Out\seiche.gif"       --field w     --fps 12 --stride 2 --dpi 80
python animate.py fr_beam   -o "$Out\beam.gif"         --field speed --fps 12 --stride 2 --dpi 80 --xlim 2500

Write-Host ""
Get-ChildItem "$Out\*.gif" | ForEach-Object {
    Write-Host ("   {0,-22} {1,7:n2} MB" -f $_.Name, ($_.Length / 1MB))
}
Write-Host ""
Write-Host "Anything over ~3 MB: raise --stride or lower --dpi and rerun." -ForegroundColor Cyan
