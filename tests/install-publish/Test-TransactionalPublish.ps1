[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Helper,
    [Parameter(Mandatory = $true)][string]$FixtureRoot
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest
. $Helper

$root = [System.IO.Path]::GetFullPath($FixtureRoot)
$testsMarker = [System.IO.Path]::DirectorySeparatorChar + '.tests' +
    [System.IO.Path]::DirectorySeparatorChar
if ($root.IndexOf($testsMarker, [System.StringComparison]::OrdinalIgnoreCase) -lt 0 -or
    (Split-Path -Leaf $root) -ne 'install-publish-windows') {
    throw "Refusing transactional publish test fixture outside .tests/.../install-publish-windows: $root"
}
if (Test-Path -LiteralPath $root) {
    Remove-Item -LiteralPath $root -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $root | Out-Null
$destination = Join-Path $root 'release'
$stageOne = Join-Path $root 'stage-one'
New-Item -ItemType Directory -Force -Path $destination, $stageOne | Out-Null
Set-Content -LiteralPath (Join-Path $destination 'payload.txt') -Value 'old' -NoNewline
Set-Content -LiteralPath (Join-Path $stageOne 'payload.txt') -Value 'new-one' -NoNewline

$validateOne = {
    param($directory)
    if ((Get-Content -LiteralPath (Join-Path $directory 'payload.txt') -Raw) -ne 'new-one') {
        throw 'Unexpected staged payload.'
    }
}
$null = Publish-TransactionalDirectory -StagedDirectory $stageOne `
    -DestinationDirectory $destination -Validate $validateOne
if ((Get-Content -LiteralPath (Join-Path $destination 'payload.txt') -Raw) -ne 'new-one') {
    throw 'Successful transactional publish did not install the staged payload.'
}

$mutationStage = Join-Path $root 'stage-mutation'
New-Item -ItemType Directory -Force -Path $mutationStage | Out-Null
Set-Content -LiteralPath (Join-Path $mutationStage 'payload.txt') -Value 'mutation' -NoNewline
$script:mutatedStage = $false
$mutateAfterIdentity = {
    param($directory)
    if (-not $script:mutatedStage) {
        Add-Content -LiteralPath (Join-Path $directory 'payload.txt') -Value '-changed' -NoNewline
        $script:mutatedStage = $true
    }
}
$mutationRejected = $false
try {
    $null = Publish-TransactionalDirectory -StagedDirectory $mutationStage `
        -DestinationDirectory $destination -Validate $mutateAfterIdentity
}
catch {
    $mutationRejected = $true
}
if (-not $mutationRejected -or
    (Get-Content -LiteralPath (Join-Path $destination 'payload.txt') -Raw) -ne 'new-one') {
    throw 'Staged content mutation was not rejected before swapping the destination.'
}

$stageTwo = Join-Path $root 'stage-two'
New-Item -ItemType Directory -Force -Path $stageTwo | Out-Null
Set-Content -LiteralPath (Join-Path $stageTwo 'payload.txt') -Value 'new-two' -NoNewline
$script:validationCalls = 0
$failAfterSwap = {
    param($directory)
    $script:validationCalls++
    if ((Get-Content -LiteralPath (Join-Path $directory 'payload.txt') -Raw) -ne 'new-two') {
        throw 'Unexpected rollback-test payload.'
    }
    if ($script:validationCalls -ge 3) {
        throw 'Injected post-swap validation failure.'
    }
}
$failed = $false
try {
    $null = Publish-TransactionalDirectory -StagedDirectory $stageTwo `
        -DestinationDirectory $destination -Validate $failAfterSwap
}
catch {
    $failed = $true
}
if (-not $failed) {
    throw 'Injected post-swap validation failure was not reported.'
}
if ((Get-Content -LiteralPath (Join-Path $destination 'payload.txt') -Raw) -ne 'new-one') {
    throw 'Transactional publish did not restore the previous destination after failure.'
}
Write-Host 'Transactional Windows publish and rollback test passed.'
