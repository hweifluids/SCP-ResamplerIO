# SCP-ResamplerIO

`SCP-ResamplerIO` owns Streamcenter+'s mesh inspection, format-adapter, and unstructured-to-structured resampling boundary. It can be built outside the Streamcenter+ GUI and exports the CMake target `SCP::ResamplerIO`.

## Components

- `SCP::ResamplerIO`: a C++17 static library containing `FlowDataInspector` and the CPU resampling core.
- `Resample_solver_<YYYYMMDD>`: an optional command-line tool with `--solver-info` identity type `tools`.
- `scp_nek5000_smoke`: a deterministic generator/validator for the repository's 20-frame Nek5000 smoke asset.
- `FlowDataTypes.h`: the shared `FlowDataKind`, `FlowDataSource`, and `FlowDataEntry` compatibility types.
- `FlowDataInspector.h`: scanning and conversion APIs used by the GUI and standalone tools.

The command-line release is always installed below this repository at `release/tools`. Build trees, downloaded dependencies, and generated releases stay under `.build`, `.deps`, and `release`; all are ignored by Git.

## Format boundary

- **VTK:** scans legacy/XML VTK datasets, PVD collections, and VTKHDF time steps. The resampler accepts VTK datasets and writes uniform `vtkImageData` frames as VTI plus a PVD collection.
- **Nek5000:** the same generic Resample solver/API accepts a `.nek5000` descriptor directly. It enumerates the `TIME_STEPS` advertised by VTK's `vtkNek5000Reader`, applies the deck's zero-based `frame_range` to that advertised sequence, reads each selected time as a `vtkDataSet`, and sends it through the selected existing interpolation method (including `vtk_resample_to_image`). Output time values are preserved in the generated PVD. `tools/src/NekFtleResampler.cpp` is the migrated legacy MPI/gslib adapter source; it is not a default release target and no longer leaves gslib paths wired into the FTLE repository.
- **OpenFOAM:** scans serial or parallel OpenFOAM cases through VTK readers and can convert selected fields to VTU/PVD output.
- **EnSight Gold:** scans one or several case files through VTK's EnSight reader and can convert selected/common fields to VTU/PVD output.

The module defines the adapter boundary; it does not vendor separate Nek5000, OpenFOAM, or EnSight parser implementations. Reader behavior and supported format revisions follow the pinned VTK dependency.

## Requirements

- CMake 3.24 or newer and a C++17 compiler.
- Git for first-time vcpkg bootstrap.
- Internet access for the first installation so the pinned vcpkg checkout and manifest dependencies can be downloaded; verified caches are reused on later runs.
- `SCP-SolverSDK`, either checked out as sibling `../SCP-SolverSDK` or installed as a CMake package.
- Qt 6.5 or newer (`Core`) and VTK, obtained from the pinned vcpkg manifest by the one-command installers.
- CUDA and MPI are treated as target-system foundation components when required by a selected VTK build or a future accelerator backend. The current Resample tool reports `cpu=1` and `cuda=0` and does not build a CUDA kernel.

## One-command installation

Windows PowerShell:

```powershell
.\install.ps1
```

Linux:

```bash
bash ./install.sh
```

The installers bootstrap vcpkg at commit `42f37e46d37baf7b218966434b91116bddf73461` when the selected fixed checkout is absent. Package versions are anchored by the same manifest baseline. No dependency path, triplet, SDK path, or environment variable is accepted.

Windows target and host dependency builds use the module overlay triplet with the Visual Studio 2022 `v143` platform toolset, dynamic CRT/libraries, and Release-only packages. Linux uses the matching `x64-linux-release` overlay triplet for target and host packages. Pass `-Clean` on Windows or `--clean` on Linux to remove only this module's `.build` tree before rebuilding. Every installation replaces `release/tools` after a successful build so stale dated executables cannot survive.

Inside the Streamcenter+ superproject, the installer automatically uses the superproject checkout at `.deps/vcpkg`; detached use automatically uses this module's `.deps/vcpkg`. In both cases installed packages remain private at `SCP-ResamplerIO/.deps/vcpkg_installed`.

`SCP-SolverSDK` is resolved from the fixed sibling `../SCP-SolverSDK` checkout/installed release when present. A detached `SCP-ResamplerIO` checkout instead bootstraps `https://github.com/hweifluids/SCP-SolverSDK.git` at the exact tag `v0.1.0` into `.deps/SCP-SolverSDK`; the `VERSION` contract and exact checked-out tag are validated before configure. If the repository or tag has not yet been published, installation stops with an explicit publication error. The installers expose no external component-path override.

## CMake consumption

As a sibling module:

```cmake
add_subdirectory(path/to/SCP-ResamplerIO)
target_link_libraries(my_consumer PRIVATE SCP::ResamplerIO)
```

As an installed package:

```cmake
find_package(SCPResamplerIO CONFIG REQUIRED)
target_link_libraries(my_consumer PRIVATE SCP::ResamplerIO)
```

For installed-package use, make the `SCP-ResamplerIO` and `SCP-SolverSDK` installation prefixes available through the CMake `CMAKE_PREFIX_PATH` configure option.

The public format-adapter boundary is declared in `streamcenterplus_resample/Core.h`:

```cpp
auto frames = streamcenterplus::resample::EnumerateDataSetFrames("case.nek5000");
auto data = streamcenterplus::resample::ReadDataSet(frames.at(0));
auto bounds = streamcenterplus::resample::ReadDataSetBounds(frames.at(0));
```

`FrameInfo` records both the physical output time and, for multi-time sources such as Nek5000, the exact reader time needed to reproduce that logical frame. File-per-frame VTK callers can continue using the path overloads.

## Smoke example

`examples/smoke` provides no-argument Windows and Linux entrypoints. It consumes the published private `SCP-examples/hemi_resampled` dataset and independently reads all 20 frames to check geometry, `Velocity(float3)`, the valid-point mask, time ordering, PVD references, and JSON syntax. The integrated superproject uses its `examples/SCP-examples` submodule; detached use resolves a sibling checkout or an authenticated private `smoke-data-v1` Release download.

The smoke runner writes its prepared data only to the fixed ignored `.tests/smoke/SCP-ResamplerIO/data` path. The original 1,000-frame Nek5000 sequence is maintainer-only provenance and is not required on a target machine.

The no-argument CTest integration target `SCPResamplerIO.Nek5000Integration` separately exercises the generic adapter and solver on advertised source frames 490 and 491 with a small `16 x 8 x 6` grid. Its output is fixed at `.tests/integration/SCP-ResamplerIO/nek5000-general-adapter`. The test is reported as skipped when the fixed source dataset is not available; it never writes into the source or release trees.
