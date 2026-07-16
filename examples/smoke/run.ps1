[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"
$moduleRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)

& (Join-Path $moduleRoot "install.ps1")
if ($LASTEXITCODE -ne 0) {
    throw "SCP-ResamplerIO installation failed with exit code $LASTEXITCODE."
}
& (Join-Path $PSScriptRoot "prepare_data.ps1")
