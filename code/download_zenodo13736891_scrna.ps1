param(
    [string]$ManifestPath = "outputs\NPC_SecAct_Zenodo13736891_scrna_manifest_v0.1.csv",
    [string]$Destination = "work\zenodo13736891_scrna_archives",
    [int]$ThrottleLimit = 3
)

$ErrorActionPreference = "Stop"
$destinationPath = [System.IO.Path]::GetFullPath($Destination)
New-Item -ItemType Directory -Force $destinationPath | Out-Null

$rows = Import-Csv $ManifestPath
$downloadRows = foreach ($row in $rows) {
    $target = Join-Path $destinationPath $row.filename
    $isComplete = $false
    if (Test-Path -LiteralPath $target) {
        $file = Get-Item -LiteralPath $target
        if ($file.Length -eq [int64]$row.size_bytes) {
            $hash = (Get-FileHash -Algorithm MD5 -LiteralPath $target).Hash.ToLowerInvariant()
            $isComplete = $hash -eq $row.md5
        }
    }
    if ($isComplete) {
        [pscustomobject]@{ row = $row; skip = $true }
    } else {
        [pscustomobject]@{ row = $row; skip = $false }
    }
}

$completed = @($downloadRows | Where-Object skip | ForEach-Object {
    [pscustomobject]@{
        filename = $_.row.filename
        status = "verified_existing"
        size_bytes = [int64]$_.row.size_bytes
        md5 = $_.row.md5
    }
})

$pending = @($downloadRows | Where-Object { -not $_.skip } | ForEach-Object { $_.row })
$downloaded = @($pending | ForEach-Object -Parallel {
    $row = $_
    $target = Join-Path $using:destinationPath $row.filename
    $resolve = "zenodo.org:443:188.184.103.118"
    Write-Host "START $($row.filename)"
    & curl.exe --silent --show-error --fail --location `
        --retry 10 --retry-all-errors --connect-timeout 20 `
        --continue-at - --resolve $resolve `
        $row.download_url --output $target
    if ($LASTEXITCODE -ne 0) {
        throw "curl failed for $($row.filename) with exit code $LASTEXITCODE"
    }

    $file = Get-Item -LiteralPath $target
    $hash = (Get-FileHash -Algorithm MD5 -LiteralPath $target).Hash.ToLowerInvariant()
    if ($file.Length -ne [int64]$row.size_bytes -or $hash -ne $row.md5) {
        throw "Integrity check failed for $($row.filename): size=$($file.Length), md5=$hash"
    }
    Write-Host "VERIFIED $($row.filename)"
    [pscustomobject]@{
        filename = $row.filename
        status = "downloaded_verified"
        size_bytes = $file.Length
        md5 = $hash
    }
} -ThrottleLimit $ThrottleLimit)

$status = @($completed + $downloaded) | Sort-Object filename
$statusPath = Join-Path $destinationPath "download_status.tsv"
$status | Export-Csv -Delimiter "`t" -NoTypeInformation -Encoding utf8 $statusPath
$status | Format-Table -AutoSize
Write-Host "STATUS $statusPath"
