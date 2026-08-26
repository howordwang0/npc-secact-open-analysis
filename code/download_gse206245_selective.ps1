param(
    [string]$FileManifest = "work/data_audit/GSE206245_filelist.txt",
    [string]$OutputDirectory = "work/gse206245_selective",
    [int]$ThrottleLimit = 2
)

$ErrorActionPreference = "Stop"
$workspace = (Get-Location).Path
$manifestPath = [IO.Path]::GetFullPath((Join-Path $workspace $FileManifest))
$outputPath = [IO.Path]::GetFullPath((Join-Path $workspace $OutputDirectory))
New-Item -ItemType Directory -Force -Path $outputPath | Out-Null

$requiredSuffixes = @(
    "_barcodes.tsv.gz",
    "_features.tsv.gz",
    "_matrix.mtx.gz",
    "_scalefactors_json.json.gz",
    "_tissue_lowres_image.png.gz",
    "_tissue_positions_list.csv.gz"
)
$gsmIds = 6248646..6248657 | ForEach-Object { "GSM$_" }

$rows = Get-Content -LiteralPath $manifestPath |
    Select-Object -Skip 1 |
    ConvertFrom-Csv -Delimiter "`t" -Header "record_type", "Name", "Time", "Size", "Type" |
    Where-Object {
        $currentName = $_.Name
        $hasRequiredSuffix = @(
            $requiredSuffixes | Where-Object { $currentName.EndsWith($_) }
        ).Count -gt 0
        $_.record_type -eq "File" -and
        ($gsmIds -contains (($currentName -split "_")[0])) -and
        $hasRequiredSuffix
    }

$expectedCount = $gsmIds.Count * $requiredSuffixes.Count
if ($rows.Count -ne $expectedCount) {
    throw "Expected $expectedCount required files in the GEO manifest, found $($rows.Count)"
}

$results = $rows | ForEach-Object -ThrottleLimit $ThrottleLimit -Parallel {
    $ErrorActionPreference = "Continue"
    $PSNativeCommandUseErrorActionPreference = $false
    $name = $_.Name
    $expectedBytes = [int64]$_.Size
    $target = Join-Path $using:outputPath $name
    $gsm = ($name -split "_")[0]
    $url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6248nnn/$gsm/suppl/$name"

    if (Test-Path -LiteralPath $target) {
        $existingBytes = (Get-Item -LiteralPath $target).Length
        if ($existingBytes -eq $expectedBytes) {
            return [pscustomobject]@{
                file = $name
                expected_bytes = $expectedBytes
                actual_bytes = $existingBytes
                status = "already_complete"
            }
        }
        if ($existingBytes -gt $expectedBytes) {
            throw "$name exceeds the GEO byte count ($existingBytes > $expectedBytes)"
        }
    }

    $downloaded = $false
    for ($attempt = 1; $attempt -le 20; $attempt++) {
        & curl.exe --ipv4 --silent --show-error --fail --location --http1.1 --ssl-no-revoke `
            --retry 3 --retry-delay 5 --retry-all-errors `
            --speed-limit 1024 --speed-time 90 --continue-at - `
            --output $target $url 2>&1 | Out-Null
        $curlExitCode = $LASTEXITCODE
        $actualAttemptBytes = if (Test-Path -LiteralPath $target) {
            (Get-Item -LiteralPath $target).Length
        } else {
            0
        }
        if ($curlExitCode -eq 0 -and $actualAttemptBytes -eq $expectedBytes) {
            $downloaded = $true
            break
        }
        Start-Sleep -Seconds ([Math]::Min(5 * $attempt, 30))
    }
    if (-not $downloaded) {
        throw "Repeated download failure for $name"
    }

    $actualBytes = (Get-Item -LiteralPath $target).Length
    if ($actualBytes -ne $expectedBytes) {
        throw "$name failed byte validation ($actualBytes != $expectedBytes)"
    }
    [pscustomobject]@{
        file = $name
        expected_bytes = $expectedBytes
        actual_bytes = $actualBytes
        status = "downloaded_and_validated"
    }
}

$statusPath = Join-Path $outputPath "download_status.tsv"
$results | Sort-Object file | Export-Csv -LiteralPath $statusPath -Delimiter "`t" -NoTypeInformation
$validated = ($results | Where-Object { $_.actual_bytes -eq $_.expected_bytes }).Count
if ($validated -ne $expectedCount) {
    throw "Only $validated of $expectedCount required files validated"
}
Write-Output "Validated $validated required GSE206245 files in $outputPath"
