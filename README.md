# SCP-ResamplerIO

`SCP-ResamplerIO` owns Streamcenter+'s mesh inspection, native dataset-reader, and unstructured-to-structured resampling boundary. It can be built outside the Streamcenter+ GUI and exports both a lightweight reader target and the full resampling target.

## Components

- `SCP::ResamplerReaders`: the Qt-free C++17 reader/numerics library used by standalone unstructured solvers.
- `SCP::ResamplerIO`: the full C++17 library containing `FlowDataInspector`, `SCP::ResamplerReaders`, and the CPU resampling core.
- `Resample_solver_<YYYYMMDD>`: an optional command-line tool with `--solver-info` identity type `tools`.
- `scp_nek5000_smoke`: a deterministic generator/validator for the repository's 20-frame Nek5000 smoke asset.
- `FlowDataTypes.h`: the shared `FlowDataKind`, `FlowDataSource`, and `FlowDataEntry` compatibility types.
- `FlowDataInspector.h`: scanning and conversion APIs used by the GUI and standalone tools.

The command-line release is always installed below this repository at `release/tools`. Build trees, downloaded dependencies, and generated releases stay under `.build`, `.deps`, and `release`; all are ignored by Git.

## Format boundary

- **VTK:** scans legacy/XML VTK datasets, PVD collections, and VTKHDF time steps. The resampler accepts VTK datasets and writes uniform `vtkImageData` frames as VTI plus a PVD collection.
- **Nek5000:** accepts `.nek5000` and `.nek` descriptors and reads the referenced `#std` field files directly, without an intermediate VTU. The reader supports 4/8-byte words, either endian, 2-D/3-D and anisotropic GLL dimensions, multi-file element maps, mesh-bearing or mesh-reusing frames, velocity, pressure, temperature, and passive scalars. Values remain element-local `double` GLL nodal values; coincident interface nodes are intentionally not merged, and elements are returned in canonical global-ID order across pieces and frames. The output can carry reference GLL weights, the signed mapping Jacobian, physical quadrature weights `w_r w_s w_t |J|`, global spectral-element IDs, and tensor dimensions. `ComputeNek5000SpectralGradient` differentiates in reference space with the GLL derivative matrices and transforms with the inverse physical Jacobian (or the 2-D tangent metric). `Nek5000SpectralPointEvaluator` inverts the curved element mapping with Newton iteration and evaluates the tensor-product GLL polynomial. Its BVH uses conservative bounds for the complete high-order mapping, including curved regions outside the GLL-node box; `bounds()` exposes the resulting global physical box. The Resample solver automatically uses this spectral interpolation path for Nek data instead of interpolating the visualization-only linear subcells.
- **OpenFOAM:** opens reconstructed or decomposed cases and their native time sequence directly through `vtkPOpenFOAMReader`; only the internal volume mesh is assembled for analysis.
- **EnSight Gold:** opens ASCII or binary Gold `.case` sources and their native time sequence directly through `vtkGenericEnSightReader`.

Nek5000 parsing and spectral-element numerics are implemented in this module. EnSight Gold and OpenFOAM parsing use the corresponding readers from the pinned VTK dependency. Composite parts are retained as separate point sets while being assembled, and `SCP_SourceBlockId` records their source block.

The implementation follows the official [Nek5000 field-file description](https://nek5000.github.io/NekDoc/problem_setup/case_files.html) for the 132-byte `#std` header, endian tag, element map and field blocks. Its numerical treatment follows Nek5000's documented [GLL/GL placement and element-local derivative semantics](https://nek5000.github.io/NekDoc/faq.html). Linear quad/hex subcells are emitted only as a VTK-compatible visualization/topology view, matching the role described by VTK's [`vtkNek5000Reader`](https://vtk.org/doc/nightly/html/classvtkNek5000Reader.html); integration, gradients and resampling use the retained spectral metadata instead.

Three Nek5000 interpretation boundaries are explicit:

- Automatic spectral metadata assumes a standard field written at Nek's GLL solution points. The `#std` header does not reliably identify `ifreguo`/uniform or another custom reference-node distribution; callers that know the output is not standard GLL must disable `attachNek5000SpectralMetadata` and must not use the GLL gradient, quadrature or point-evaluator APIs.
- For 2-D data the generated Jacobian/quadrature is planar or embedded-surface area measure. An axisymmetric `2 pi r` volume factor cannot be inferred from a field file and must be supplied by the analysis layer through an explicit weight array, which takes precedence over automatic GLL weights.
- A pressure block is read exactly as the output point values stored in the field file, whose documented shape is `(nelt, nxyz)`. The reader does not reinterpret or reconstruct it on the solver's internal lower-order Pn-Pn-2 pressure grid.

## Requirements

- CMake 3.24 or newer and a C++17 compiler.
- Git for first-time vcpkg bootstrap.
- Internet access for the first installation so the pinned vcpkg checkout and manifest dependencies can be downloaded; verified caches are reused on later runs.
- `SCP-SolverSDK`, either checked out as sibling `../SCP-SolverSDK` or installed as a CMake package.
- VTK, obtained from the pinned vcpkg manifest by the one-command installers. Qt 6.5 or newer (`Core`) is required only for the full `SCP::ResamplerIO` target, not `SCP::ResamplerReaders`.
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

For a standalone solver that only needs native readers and Nek5000 numerics:

```cmake
set(SCP_RESAMPLER_IO_BUILD_FULL_LIBRARY OFF CACHE BOOL "" FORCE)
set(SCP_RESAMPLER_IO_BUILD_TOOL OFF CACHE BOOL "" FORCE)
set(SCP_RESAMPLER_IO_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(path/to/SCP-ResamplerIO)
target_link_libraries(my_solver PRIVATE SCP::ResamplerReaders)
```

The unstructured solver repositories use the sibling checkout when available and otherwise use CMake `FetchContent` to clone the public `https://github.com/hweifluids/SCP-ResamplerIO.git` repository. `SCP_RESAMPLER_IO_GIT_REPOSITORY` and `SCP_RESAMPLER_IO_GIT_TAG` are cache variables so release builds can pin a published tag or commit. Consequently a detached solver install automatically downloads this dependency when it is absent; that requires Git and network access on the first configure.

As an installed full package:

```cmake
find_package(SCPResamplerIO CONFIG REQUIRED)
target_link_libraries(my_consumer PRIVATE SCP::ResamplerIO)
```

Installed lightweight consumers can skip Qt discovery by requesting only the
reader component:

```cmake
find_package(SCPResamplerIO CONFIG REQUIRED COMPONENTS Readers)
target_link_libraries(my_solver PRIVATE SCP::ResamplerReaders)
```

For installed-package use, make the `SCP-ResamplerIO` and `SCP-SolverSDK` installation prefixes available through the CMake `CMAKE_PREFIX_PATH` configure option.

The shared public reader boundary is declared in `streamcenterplus/io/DataSetReader.h`:

```cpp
auto frames = streamcenterplus::io::EnumerateDataSetFrames("case.nek5000");
streamcenterplus::io::ReadOptions options;
options.buildNek5000LinearizedTopology = false;
options.attachNek5000QuadratureArrays = false;
auto data = streamcenterplus::io::ReadDataSet(frames.at(0), options);
auto grad = streamcenterplus::io::ComputeNek5000SpectralGradient(
    data, data->GetPointData()->GetArray("Velocity"));
```

`DataSetFrame` records physical time, native reader time, source kind and source index. File-per-frame VTK callers can continue using the path overloads. `ReadOptions::buildNek5000LinearizedTopology` can suppress visualization-only GLL subcells, while the master quadrature switch and the reference/Jacobian/physical selectors can suppress any unused point-sized arrays. Basic GLL dimensions, canonical element IDs and representation metadata remain attached when these optional representations are disabled. `HasSameNek5000SpectralGeometry` compares tensor dimensions, canonical IDs and coordinates before an evaluator is rebound to another time-frame field.

## Smoke example

`examples/smoke` provides no-argument Windows and Linux entrypoints. It consumes the published private `SCP-examples/hemi_resampled` dataset and independently reads all 20 frames to check geometry, `Velocity(float3)`, the valid-point mask, time ordering, PVD references, and JSON syntax. The integrated superproject uses its `examples/SCP-examples` submodule; detached use resolves a sibling checkout or an authenticated private `smoke-data-v1` Release download.

The smoke runner writes its prepared data only to the fixed ignored `.tests/smoke/SCP-ResamplerIO/data` path. The original 1,000-frame Nek5000 sequence is maintainer-only provenance and is not required on a target machine.

The no-argument CTest integration target `SCPResamplerIO.Nek5000Integration` separately exercises the generic adapter and solver on advertised source frames 490 and 491 with a small `16 x 8 x 6` grid. Its output is fixed at `.tests/integration/SCP-ResamplerIO/nek5000-general-adapter`. The test is reported as skipped when the fixed source dataset is not available; it never writes into the source or release trees.
