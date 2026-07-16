[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$cmake = Get-Command cmake -ErrorAction Stop
& $cmake.Source "-DSMOKE_ROOT=$PSScriptRoot" -P (Join-Path $PSScriptRoot "prepare-data.cmake")
if ($LASTEXITCODE -ne 0) {
    throw "SCP-examples smoke data preparation failed with exit code $LASTEXITCODE."
}
& (Join-Path $PSScriptRoot "validate.ps1")
