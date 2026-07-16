[CmdletBinding()]
param()

$ErrorActionPreference = "Stop"

function Resolve-ComponentSuperprojectRoot {
    param([string]$ModuleRoot, [string]$ComponentName)
    $modulesRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $ModuleRoot))
    if ((Split-Path -Leaf $modulesRoot) -ne "modules") { return $null }
    $candidateRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $modulesRoot))
    $manifestPath = Join-Path $modulesRoot "components.json"
    if (-not (Test-Path -LiteralPath (Join-Path $candidateRoot ".git")) -or
        -not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) { return $null }
    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    $component = @($manifest.components | Where-Object { $_.name -eq $ComponentName }) | Select-Object -First 1
    if ($null -eq $component) { return $null }
    $expectedMount = [System.IO.Path]::GetFullPath((Join-Path $modulesRoot ([string]$component.path)))
    if ([System.StringComparer]::OrdinalIgnoreCase.Equals(
            [System.IO.Path]::GetFullPath($ModuleRoot), $expectedMount)) { return $candidateRoot }
    return $null
}

$moduleRoot = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$tool = Join-Path $moduleRoot "release\tools\scp_nek5000_smoke.exe"
$superprojectRoot = Resolve-ComponentSuperprojectRoot -ModuleRoot $moduleRoot -ComponentName "SCP-ResamplerIO"
$testRoot = if ($null -ne $superprojectRoot) { Join-Path $superprojectRoot ".tests" } else { Join-Path $moduleRoot ".tests" }
$datasetDirectory = Join-Path $testRoot "smoke\SCP-ResamplerIO\data"

if (-not (Test-Path -LiteralPath $tool -PathType Leaf)) {
    throw "The smoke tool is not installed at its fixed location: $tool. Run ../../install.ps1 first."
}
& $tool validate $datasetDirectory
if ($LASTEXITCODE -ne 0) {
    throw "Nek5000 smoke validation failed with exit code $LASTEXITCODE."
}
