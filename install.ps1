[CmdletBinding()]
param(
    [switch]$Clean
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

function Resolve-ComponentSuperprojectRoot {
    param([string]$ModuleRoot, [string]$ComponentName)

    $modulesRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $ModuleRoot))
    if ((Split-Path -Leaf $modulesRoot) -ne "modules") {
        return $null
    }
    $candidateRoot = [System.IO.Path]::GetFullPath((Split-Path -Parent $modulesRoot))
    $manifestPath = Join-Path $modulesRoot "components.json"
    if (-not (Test-Path -LiteralPath (Join-Path $candidateRoot ".git")) -or
        -not (Test-Path -LiteralPath $manifestPath -PathType Leaf)) {
        return $null
    }

    $manifest = Get-Content -LiteralPath $manifestPath -Raw | ConvertFrom-Json
    $component = @($manifest.components | Where-Object { $_.name -eq $ComponentName }) | Select-Object -First 1
    if ($null -eq $component) {
        return $null
    }
    $expectedMount = [System.IO.Path]::GetFullPath((Join-Path $modulesRoot ([string]$component.path)))
    if ([System.StringComparer]::OrdinalIgnoreCase.Equals(
            [System.IO.Path]::GetFullPath($ModuleRoot), $expectedMount)) {
        return $candidateRoot
    }
    return $null
}

$moduleRoot = [System.IO.Path]::GetFullPath($PSScriptRoot)
. (Join-Path $moduleRoot "cmake\TransactionalPublish.ps1")
$vcpkgCommit = "42f37e46d37baf7b218966434b91116bddf73461"
$vcpkgRepository = "https://github.com/microsoft/vcpkg.git"
$sdkRepository = "https://github.com/hweifluids/SCP-SolverSDK.git"
$sdkTag = "779232d3ad4eaca5d359db5917b56a3c47289cde"
$sdkVersion = "0.2.0"
$BuildDir = Join-Path $moduleRoot ".build\windows-x64"
$InstallDir = Join-Path $moduleRoot "release"
$VcpkgInstalledDir = Join-Path $moduleRoot ".deps\vcpkg_installed"
$VcpkgTriplet = "x64-windows-release"
$Generator = "Visual Studio 17 2022"
$Architecture = "x64"
$Toolset = "v143"

$superprojectRoot = Resolve-ComponentSuperprojectRoot -ModuleRoot $moduleRoot -ComponentName "SCP-ResamplerIO"
$VcpkgRoot = if ($null -ne $superprojectRoot) {
    Join-Path $superprojectRoot ".deps\vcpkg"
} else {
    Join-Path $moduleRoot ".deps\vcpkg"
}

function Invoke-NativeCommand {
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments,
        [Parameter(Mandatory = $true)][string]$FailureMessage
    )

    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "$FailureMessage Exit code: $LASTEXITCODE"
    }
}

$BuildDir = [System.IO.Path]::GetFullPath($BuildDir)
$InstallDir = [System.IO.Path]::GetFullPath($InstallDir)
$VcpkgRoot = [System.IO.Path]::GetFullPath($VcpkgRoot)
$VcpkgInstalledDir = [System.IO.Path]::GetFullPath($VcpkgInstalledDir)

if ($Clean) {
    $moduleBuildRoot = [System.IO.Path]::GetFullPath((Join-Path $moduleRoot ".build"))
    if (-not $moduleBuildRoot.StartsWith($moduleRoot + [System.IO.Path]::DirectorySeparatorChar,
                                         [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to clean a build directory outside SCP-ResamplerIO: $moduleBuildRoot"
    }
    if (Test-Path -LiteralPath $moduleBuildRoot) {
        Remove-Item -LiteralPath $moduleBuildRoot -Recurse -Force
    }
}

$cmake = Get-Command cmake -ErrorAction SilentlyContinue
if ($null -eq $cmake) {
    throw "CMake 3.24 or newer is required and was not found on PATH."
}
$git = Get-Command git -ErrorAction SilentlyContinue
if ($null -eq $git) {
    throw "Git is required for the pinned vcpkg bootstrap and was not found on PATH."
}
$windowsPowerShell = Get-Command powershell.exe -ErrorAction SilentlyContinue
if ($null -eq $windowsPowerShell) {
    throw "Windows PowerShell is required by vcpkg's MSBuild deployment step."
}

$vcpkgExe = Join-Path $VcpkgRoot "vcpkg.exe"
$vcpkgToolchain = Join-Path $VcpkgRoot "scripts\buildsystems\vcpkg.cmake"
if (-not (Test-Path -LiteralPath $vcpkgExe) -or -not (Test-Path -LiteralPath $vcpkgToolchain)) {
    if (Test-Path -LiteralPath $VcpkgRoot) {
        throw "The selected vcpkg root is incomplete: $VcpkgRoot. Remove or repair it, then retry."
    }

    New-Item -ItemType Directory -Force -Path $VcpkgRoot | Out-Null
    Invoke-NativeCommand -Executable $git.Source `
        -Arguments @("-C", $VcpkgRoot, "init") `
        -FailureMessage "Could not initialize the module-local vcpkg repository."
    Invoke-NativeCommand -Executable $git.Source `
        -Arguments @("-C", $VcpkgRoot, "remote", "add", "origin", $vcpkgRepository) `
        -FailureMessage "Could not configure the vcpkg origin."
    Invoke-NativeCommand -Executable $git.Source `
        -Arguments @("-C", $VcpkgRoot, "fetch", "--depth", "1", "origin", $vcpkgCommit) `
        -FailureMessage "Could not fetch pinned vcpkg commit $vcpkgCommit."
    Invoke-NativeCommand -Executable $git.Source `
        -Arguments @("-C", $VcpkgRoot, "checkout", "--detach", "FETCH_HEAD") `
        -FailureMessage "Could not check out pinned vcpkg commit $vcpkgCommit."
    Invoke-NativeCommand -Executable (Join-Path $VcpkgRoot "bootstrap-vcpkg.bat") `
        -Arguments @("-disableMetrics") `
        -FailureMessage "Could not bootstrap vcpkg."
}

$actualVcpkgOrigin = ([string](& $git.Source -C $VcpkgRoot remote get-url origin 2>$null)).Trim()
if ($LASTEXITCODE -ne 0 -or $actualVcpkgOrigin -ne $vcpkgRepository) {
    throw "The selected vcpkg checkout must use the pinned origin ${vcpkgRepository}: $VcpkgRoot"
}
$actualVcpkgCommit = ([string](& $git.Source -C $VcpkgRoot rev-parse HEAD 2>$null)).Trim().ToLowerInvariant()
if ($LASTEXITCODE -ne 0 -or $actualVcpkgCommit -ne $vcpkgCommit) {
    throw "The selected vcpkg checkout must be exactly commit ${vcpkgCommit}: $VcpkgRoot"
}
$vcpkgStatus = (& $git.Source -C $VcpkgRoot status --porcelain 2>$null) -join "`n"
if ($LASTEXITCODE -ne 0 -or -not [string]::IsNullOrWhiteSpace($vcpkgStatus)) {
    throw "The selected vcpkg checkout must be clean: $VcpkgRoot"
}

$sdkSiblingRoot = Join-Path (Split-Path -Parent $moduleRoot) "SCP-SolverSDK"
$sdkSource = Join-Path $sdkSiblingRoot "CMakeLists.txt"
$sdkLocalRoot = Join-Path $moduleRoot ".deps\SCP-SolverSDK"
$sdkLocalSource = Join-Path $sdkLocalRoot "CMakeLists.txt"
$sdkInstalledConfigs = @(
    Get-ChildItem -Path (Join-Path $sdkSiblingRoot "release\*\lib\cmake\SCPSolverSDK\SCPSolverSDKConfig.cmake") -File -ErrorAction SilentlyContinue
    Get-ChildItem -Path (Join-Path $sdkLocalRoot "release\*\lib\cmake\SCPSolverSDK\SCPSolverSDKConfig.cmake") -File -ErrorAction SilentlyContinue
)
if (-not (Test-Path -LiteralPath $sdkSource) -and
    -not (Test-Path -LiteralPath $sdkLocalSource) -and
    $sdkInstalledConfigs.Count -eq 0) {
    if (Test-Path -LiteralPath $sdkLocalRoot) {
        throw "The fixed module-local SDK checkout is incomplete: $sdkLocalRoot. Remove or repair it, then retry."
    }
    New-Item -ItemType Directory -Force -Path (Split-Path -Parent $sdkLocalRoot) | Out-Null
    try {
        Invoke-NativeCommand -Executable $git.Source `
            -Arguments @("init", $sdkLocalRoot) -FailureMessage "Could not initialize SDK checkout."
        Invoke-NativeCommand -Executable $git.Source `
            -Arguments @("-C", $sdkLocalRoot, "remote", "add", "origin", $sdkRepository) `
            -FailureMessage "Could not configure SDK origin."
        Invoke-NativeCommand -Executable $git.Source `
            -Arguments @("-C", $sdkLocalRoot, "fetch", "--depth", "1", "origin", $sdkTag) `
            -FailureMessage "Could not fetch SCP-SolverSDK $sdkTag."
        Invoke-NativeCommand -Executable $git.Source `
            -Arguments @("-C", $sdkLocalRoot, "checkout", "--detach", $sdkTag) `
            -FailureMessage "Could not check out SCP-SolverSDK $sdkTag."
    }
    catch {
        throw "Could not bootstrap SCP-SolverSDK from $sdkRepository revision $sdkTag. The revision must be published before detached installation can succeed. Details: $($_.Exception.Message)"
    }
}

$sdkSourceRoot = if (Test-Path -LiteralPath $sdkSource) {
    $sdkSiblingRoot
} elseif (Test-Path -LiteralPath $sdkLocalSource) {
    $sdkLocalRoot
} else {
    $null
}
if ($null -ne $sdkSourceRoot) {
    $sdkVersionFile = Join-Path $sdkSourceRoot "VERSION"
    if (-not (Test-Path -LiteralPath $sdkVersionFile -PathType Leaf) -or
        (Get-Content -LiteralPath $sdkVersionFile -Raw).Trim() -ne $sdkVersion) {
        throw "SCP-ResamplerIO requires SCP-SolverSDK $sdkVersion at $sdkSourceRoot."
    }
}
if ([System.StringComparer]::OrdinalIgnoreCase.Equals($sdkSourceRoot, $sdkLocalRoot)) {
    $actualSdkOrigin = ([string](& $git.Source -C $sdkLocalRoot remote get-url origin 2>$null)).Trim()
    if ($LASTEXITCODE -ne 0 -or $actualSdkOrigin -ne $sdkRepository) {
        throw "The module-local SCP-SolverSDK checkout must use origin ${sdkRepository}: $sdkLocalRoot"
    }
    $checkedOutTag = (& $git.Source -C $sdkLocalRoot rev-parse HEAD 2>$null).Trim()
    if ($LASTEXITCODE -ne 0 -or $checkedOutTag -ne $sdkTag) {
        throw "The module-local SCP-SolverSDK checkout must be exactly revision ${sdkTag}: $sdkLocalRoot"
    }
}

$configureArgs = @(
    "-S", $moduleRoot,
    "-B", $BuildDir,
    "-G", $Generator,
    "-DCMAKE_BUILD_TYPE=Release",
    "-DCMAKE_TRY_COMPILE_CONFIGURATION=Release",
    "-DCMAKE_INSTALL_PREFIX=$InstallDir",
    "-DCMAKE_TOOLCHAIN_FILE=$vcpkgToolchain",
    "-DVCPKG_TARGET_TRIPLET=$VcpkgTriplet",
    "-DVCPKG_HOST_TRIPLET=$VcpkgTriplet",
    "-DVCPKG_OVERLAY_TRIPLETS=$(Join-Path $moduleRoot 'cmake\triplets')",
    "-DVCPKG_INSTALLED_DIR=$VcpkgInstalledDir",
    "-DZ_VCPKG_POWERSHELL_PATH=$($windowsPowerShell.Source)",
    "-DSCP_RESAMPLER_IO_BUILD_TOOL=ON"
)
if ($Generator -like "Visual Studio*") {
    $configureArgs += @("-A", $Architecture, "-T", $Toolset)
}
Invoke-NativeCommand -Executable $cmake.Source `
    -Arguments $configureArgs `
    -FailureMessage "SCP-ResamplerIO CMake configure failed."
Invoke-NativeCommand -Executable $cmake.Source `
    -Arguments @("--build", $BuildDir, "--config", "Release", "--parallel") `
    -FailureMessage "SCP-ResamplerIO Release build failed."

function Assert-ResamplerReleaseIdentity {
    param([Parameter(Mandatory = $true)][string]$ReleaseRoot)

    $toolsRoot = Join-Path $ReleaseRoot "tools"
    $tools = @(Get-ChildItem -LiteralPath $toolsRoot -Filter "Resample_solver_*.exe" `
        -File -ErrorAction SilentlyContinue)
    if ($tools.Count -ne 1 -or $tools[0].Length -le 0) {
        throw "Validated release must contain exactly one non-empty dated Resample tool under $toolsRoot."
    }
    $identity = (& $tools[0].FullName --solver-info | Out-String)
    if ($LASTEXITCODE -ne 0) {
        throw "The staged Resample tool failed its identity command. Exit code: $LASTEXITCODE"
    }
    $requiredLines = @(
        "streamcenterplus_solver_identity=1",
        "product=Streamcenter+",
        "type=tools",
        "solver=Resample",
        "version=0.3.0",
        "cpu=1",
        "cuda=0"
    )
    $identityLines = @($identity -split "`r?`n" | Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    foreach ($required in $requiredLines) {
        if ($identityLines -cnotcontains $required) {
            throw "The staged Resample identity is missing exact line '$required'."
        }
    }
}

$stageParent = [System.IO.Path]::GetFullPath((Join-Path $BuildDir "publish-stage"))
$stagedRelease = Join-Path $stageParent "release"
if (-not $stageParent.StartsWith($BuildDir + [System.IO.Path]::DirectorySeparatorChar,
                                 [System.StringComparison]::OrdinalIgnoreCase)) {
    throw "Refusing to stage installation outside the module build tree: $stageParent"
}
if (Test-Path -LiteralPath $stageParent) {
    Remove-Item -LiteralPath $stageParent -Recurse -Force
}
New-Item -ItemType Directory -Force -Path $stageParent | Out-Null
try {
    Invoke-NativeCommand -Executable $cmake.Source `
        -Arguments @("--install", $BuildDir, "--config", "Release", "--prefix", $stagedRelease) `
        -FailureMessage "SCP-ResamplerIO staged Release install failed."
    Assert-ResamplerReleaseIdentity -ReleaseRoot $stagedRelease
    $validatedIdentity = Publish-TransactionalDirectory `
        -StagedDirectory $stagedRelease `
        -DestinationDirectory $InstallDir `
        -Validate ${function:Assert-ResamplerReleaseIdentity}
    Write-Host "Validated release content SHA-256: $validatedIdentity"
}
finally {
    if (Test-Path -LiteralPath $stageParent) {
        Remove-Item -LiteralPath $stageParent -Recurse -Force
    }
}

$tool = @(Get-ChildItem -LiteralPath (Join-Path $InstallDir "tools") `
    -Filter "Resample_solver_*.exe" -File)[0]

Write-Host "SCP-ResamplerIO installed to $InstallDir"
Write-Host "Resample tool ready at $($tool.FullName)"
