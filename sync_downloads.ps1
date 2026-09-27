# sync_downloads.ps1 -- move the newest downloaded copy of each tracked file
# into the repo, overwriting, and clean up the browser's " (1)" duplicates.
#
#   .\sync_downloads.ps1                 # sync everything it recognises
#   .\sync_downloads.ps1 iwbcurv.py      # just these
#   .\sync_downloads.ps1 -WhatIf         # show what would move, change nothing
#
# The browser will not overwrite: a second download of iwbcurv.py becomes
# "iwbcurv (1).py". This takes the most recent of those, drops it on the repo
# copy, and deletes the stragglers.

[CmdletBinding(SupportsShouldProcess = $true)]
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$Only,
    [string]$Downloads = (Join-Path $env:USERPROFILE 'Downloads'),
    [int]$WithinHours = 24
)

# Files this repo expects to receive. Anything not listed is left alone.
$tracked = @(
    'iwbcurv.py', 'animate.py', 'advect_test.py', 'phase_speed.py',
    'mry_setup.py', 'mode1.py', 'view.py', 'centroid_track.py', 'refit.py',
    'curvature.py', 'solver_bench.py', 'style_check.py',
    'ROADMAP.md', 'README.md', 'STYLE.md'
) + (Get-ChildItem -Path . -Filter '*.ps1' | ForEach-Object { $_.Name })

if ($Only) {
    $tracked = $Only
}

$cutoff = (Get-Date).AddHours(-$WithinHours)
$moved = 0
$skipped = 0

foreach ($name in ($tracked | Select-Object -Unique)) {
    $stem = [System.IO.Path]::GetFileNameWithoutExtension($name)
    $ext = [System.IO.Path]::GetExtension($name)
    # matches both "iwbcurv.py" and "iwbcurv (3).py"
    $pattern = "$stem*$ext"
    $cands = Get-ChildItem -Path $Downloads -Filter $pattern -File -ErrorAction SilentlyContinue |
        Where-Object {
            $_.LastWriteTime -gt $cutoff -and
            ($_.BaseName -eq $stem -or $_.BaseName -match "^$([regex]::Escape($stem)) \(\d+\)$")
        }
    if (-not $cands) {
        continue
    }
    $newest = $cands | Sort-Object LastWriteTime -Descending | Select-Object -First 1
    $dest = Join-Path (Get-Location) $name

    if ((Test-Path $dest) -and
        (Get-FileHash $newest.FullName).Hash -eq (Get-FileHash $dest).Hash) {
        Write-Host ("   {0,-22} unchanged" -f $name) -ForegroundColor DarkGray
        $cands | Remove-Item -Force -ErrorAction SilentlyContinue
        $skipped++
        continue
    }

    if ($PSCmdlet.ShouldProcess($name, "overwrite from $($newest.Name)")) {
        Copy-Item $newest.FullName $dest -Force
        $cands | Remove-Item -Force -ErrorAction SilentlyContinue
        Write-Host ("   {0,-22} <- {1}  ({2:HH:mm})" -f $name, $newest.Name, $newest.LastWriteTime) -ForegroundColor Green
        $moved++
    }
}

Write-Host ""
Write-Host "$moved updated, $skipped already current" -ForegroundColor Cyan
if ($moved -gt 0) {
    Write-Host "check before committing:  git status --short" -ForegroundColor Cyan
}
