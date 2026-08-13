Set-StrictMode -Version Latest

function Get-TransactionalFileSha256 {
    [CmdletBinding()]
    param([Parameter(Mandatory = $true)][string]$Path)

    $stream = [System.IO.File]::Open(
        [System.IO.Path]::GetFullPath($Path),
        [System.IO.FileMode]::Open,
        [System.IO.FileAccess]::Read,
        [System.IO.FileShare]::Read)
    $sha = [System.Security.Cryptography.SHA256]::Create()
    try {
        return (($sha.ComputeHash($stream) | ForEach-Object { $_.ToString('x2') }) -join '')
    }
    finally {
        $sha.Dispose()
        $stream.Dispose()
    }
}

function Get-TransactionalDirectoryIdentity {
    [CmdletBinding()]
    param([Parameter(Mandatory = $true)][string]$Directory)

    $root = [System.IO.Path]::GetFullPath($Directory)
    if (-not (Test-Path -LiteralPath $root -PathType Container)) {
        throw "Publish staging directory does not exist: $root"
    }
    $rootItem = Get-Item -LiteralPath $root -Force
    if (($rootItem.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
        throw "Publish staging root must not be a reparse point: $root"
    }
    $entries = [System.Collections.Generic.List[System.IO.FileSystemInfo]]::new()
    $pending = [System.Collections.Generic.Stack[string]]::new()
    $pending.Push($root)
    while ($pending.Count -gt 0) {
        $directory = $pending.Pop()
        foreach ($entry in @(Get-ChildItem -LiteralPath $directory -Force)) {
            if (($entry.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
                throw "Publish staging content must not contain reparse points: $($entry.FullName)"
            }
            $entries.Add($entry)
            if ($entry.PSIsContainer) {
                $pending.Push($entry.FullName)
            }
        }
    }
    $relativeName = {
        param([string]$FullName)
        if (-not $FullName.StartsWith($root + [System.IO.Path]::DirectorySeparatorChar,
                                      [System.StringComparison]::OrdinalIgnoreCase)) {
            throw "Publish inventory entry escapes its root: $FullName"
        }
        return $FullName.Substring($root.Length + 1).Replace('\', '/')
    }
    $files = @($entries | Where-Object { -not $_.PSIsContainer } |
        Sort-Object { & $relativeName $_.FullName })
    if ($files.Count -eq 0) {
        throw "Publish staging directory contains no files: $root"
    }

    $builder = [System.Text.StringBuilder]::new()
    foreach ($file in $files) {
        $relative = & $relativeName $file.FullName
        $hash = Get-TransactionalFileSha256 -Path $file.FullName
        [void]$builder.Append($relative.Length).Append(':').Append($relative).Append('|')
        [void]$builder.Append($file.Length).Append('|').Append($hash).Append("`n")
    }
    $bytes = [System.Text.Encoding]::UTF8.GetBytes($builder.ToString())
    $sha = [System.Security.Cryptography.SHA256]::Create()
    try {
        return (($sha.ComputeHash($bytes) | ForEach-Object { $_.ToString('x2') }) -join '')
    }
    finally {
        $sha.Dispose()
    }
}

function Publish-TransactionalDirectory {
    [CmdletBinding()]
    param(
        [Parameter(Mandatory = $true)][string]$StagedDirectory,
        [Parameter(Mandatory = $true)][string]$DestinationDirectory,
        [Parameter(Mandatory = $true)][scriptblock]$Validate
    )

    $stage = [System.IO.Path]::GetFullPath($StagedDirectory)
    $destination = [System.IO.Path]::GetFullPath($DestinationDirectory)
    $parent = [System.IO.Path]::GetFullPath((Split-Path -Parent $destination))
    if (-not (Test-Path -LiteralPath $stage -PathType Container)) {
        throw "Publish staging directory does not exist: $stage"
    }
    if (-not (Test-Path -LiteralPath $parent -PathType Container)) {
        throw "Publish destination parent does not exist: $parent"
    }
    if ([System.StringComparer]::OrdinalIgnoreCase.Equals($stage, $destination) -or
        [System.StringComparer]::OrdinalIgnoreCase.Equals($destination, $parent) -or
        $stage.StartsWith($destination + [System.IO.Path]::DirectorySeparatorChar,
                          [System.StringComparison]::OrdinalIgnoreCase) -or
        $destination.StartsWith($stage + [System.IO.Path]::DirectorySeparatorChar,
                                [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Publish staging and destination directories are not safely separated."
    }
    if (Test-Path -LiteralPath $destination) {
        $destinationItem = Get-Item -LiteralPath $destination -Force
        if (-not $destinationItem.PSIsContainer -or
            ($destinationItem.Attributes -band [System.IO.FileAttributes]::ReparsePoint) -ne 0) {
            throw "Existing publish destination must be an ordinary directory: $destination"
        }
    }

    $identity = Get-TransactionalDirectoryIdentity -Directory $stage
    & $Validate $stage

    $token = "{0}.{1}" -f $PID, ([Guid]::NewGuid().ToString('N'))
    $leaf = Split-Path -Leaf $destination
    $candidate = Join-Path $parent (".$leaf.publish.$token")
    $backup = Join-Path $parent (".$leaf.rollback.$token")
    $failed = Join-Path $parent (".$leaf.failed.$token")
    foreach ($temporaryPath in @($candidate, $backup, $failed)) {
        if (Test-Path -LiteralPath $temporaryPath) {
            throw "Refusing to overwrite an existing transactional publish path: $temporaryPath"
        }
    }

    $hadDestination = Test-Path -LiteralPath $destination
    try {
        Move-Item -LiteralPath $stage -Destination $candidate
        if ((Get-TransactionalDirectoryIdentity -Directory $candidate) -ne $identity) {
            throw "Staged publish content changed while moving into the swap location."
        }
        & $Validate $candidate

        if ($hadDestination) {
            Move-Item -LiteralPath $destination -Destination $backup
        }
        try {
            Move-Item -LiteralPath $candidate -Destination $destination
            if ((Get-TransactionalDirectoryIdentity -Directory $destination) -ne $identity) {
                throw "Published directory content identity differs from validated staging."
            }
            & $Validate $destination
        }
        catch {
            if (Test-Path -LiteralPath $destination) {
                Move-Item -LiteralPath $destination -Destination $failed
            }
            if ($hadDestination -and (Test-Path -LiteralPath $backup)) {
                Move-Item -LiteralPath $backup -Destination $destination
            }
            throw
        }

        if ($hadDestination -and (Test-Path -LiteralPath $backup)) {
            Remove-Item -LiteralPath $backup -Recurse -Force
        }
        if (Test-Path -LiteralPath $failed) {
            Remove-Item -LiteralPath $failed -Recurse -Force
        }
        return $identity
    }
    catch {
        if (-not (Test-Path -LiteralPath $destination) -and
            $hadDestination -and (Test-Path -LiteralPath $backup)) {
            Move-Item -LiteralPath $backup -Destination $destination
        }
        throw
    }
    finally {
        foreach ($temporaryPath in @($candidate, $failed)) {
            if (Test-Path -LiteralPath $temporaryPath) {
                Remove-Item -LiteralPath $temporaryPath -Recurse -Force
            }
        }
    }
}
