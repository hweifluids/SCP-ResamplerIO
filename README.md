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

Resample result publication is namespace- and manifest-bound. `output_basename`
must be one portable filename component (no absolute root, separator, `..`, or
reserved filename characters). Before changing an existing result, the solver
fully loads and fingerprints the selected input sequence. Its
`<basename>.scp-resample-manifest` binds ordered source content, frame times,
array schema, grid topology/coordinates, every interpolation and domain
transform setting, implementation version, and SHA-256 checksums for published
VTI/PVD/summary files. `solver_status=error` reuses only a contiguous checksum-
validated prefix with the same complete identity. Cleanup never sweeps unrelated
VTI/PVD files; without a valid ownership manifest, an occupied basename
namespace is preserved and the run fails closed. Frame files, PVD, summary, and
the manifest are published through same-directory temporary files and atomic
replacement.

Kriging uses a dense local system and therefore accepts at most 1024 neighbors.
PVD `file` attributes are XML-entity decoded, and explicit `timestep` values
must be complete finite numbers.

## Format boundary

- **VTK:** scans legacy/XML VTK datasets, PVD collections, and VTKHDF time steps. The resampler accepts VTK datasets and writes uniform `vtkImageData` frames as VTI plus a PVD collection.
- **Tecplot:** accepts ASCII `.dat` and classic binary `.plt` sources, including `ORDERED` zones (loaded as `vtkStructuredGrid`), fixed finite-element line, triangle, quadrilateral, tetrahedron, and brick zones, plus face-based `FEPOLYGON` and `FEPOLYHEDRON` zones (loaded as VTK polygon/polyhedron cells). Each file is exposed through the same single-file frame API used by the raw solvers and can be passed directly to the Resample tool or `FlowDataInspector::scanFolder`; folders are ordered as normal frame series. Classic `.plt` supports its documented v112/v191 layouts, byte order, nodal/cell variables, shared variables/connectivity, and shared face maps. Tecplot subzone `.szplt` and solution-only `.plt` files are intentionally not supported.
- **Nek5000:** accepts `.nek5000` and `.nek` descriptors and reads the referenced `#std` field files directly, without an intermediate VTU. Standard `visnek` file templates with one, two, or three integer conversions are supported, including multi-output-directory layouts such as `A%01d/case%01d.f%05d`. The reader supports 4/8-byte words, either endian, 2-D/3-D and anisotropic GLL dimensions, multi-file element maps, mesh-bearing or mesh-reusing frames, velocity, pressure, temperature, and passive scalars. Header dimensions, element/file counts, offsets, file extent, and all allocation arithmetic are checked before allocation; explicit resource ceilings reject unreasonable headers while retaining the fixed large integration corpus. Sparse positive global element IDs are supported and duplicates are rejected across pieces. Values remain element-local `double` GLL nodal values; coincident interface nodes are intentionally not merged, and elements are returned in canonical global-ID order across pieces and frames. The output can carry reference GLL weights, the signed mapping Jacobian, physical quadrature weights `w_r w_s w_t |J|`, global spectral-element IDs, and tensor dimensions. `ComputeNek5000SpectralGradient` differentiates in reference space with the GLL derivative matrices and transforms with the inverse physical Jacobian (or the 2-D tangent metric). `Nek5000SpectralPointEvaluator` inverts the curved element mapping with Newton iteration and evaluates the tensor-product GLL polynomial. Its BVH uses conservative bounds for the complete high-order mapping, including curved regions outside the GLL-node box; `bounds()` exposes the resulting global physical box. The Resample solver automatically uses this spectral interpolation path for Nek data instead of interpolating the visualization-only linear subcells.
- **OpenFOAM:** opens reconstructed or decomposed cases and their native time sequence directly through `vtkPOpenFOAMReader`; only the internal volume mesh is assembled for analysis.
- **EnSight Gold:** opens ASCII or binary Gold `.case` sources and their native time sequence directly through `vtkGenericEnSightReader`.

Nek5000 parsing and spectral-element numerics are implemented in this module. EnSight Gold and OpenFOAM parsing use the corresponding readers from the pinned VTK dependency. Composite parts are retained as separate point sets while being assembled, and `SCP_SourceBlockId` records their source block.

### Two-zone `.scpm` manifests

The reader and Resample tool also accept a lightweight `key=value` `.scpm`
mesh-set source. It is deliberately limited to exactly two zones and projects
them onto one fixed world-reference mesh: the larger reference zone is the
background (or select it explicitly with `background_zone=0|1`), while the
other zone supplies values at covered background points. This keeps the point
count, connectivity and world coordinates invariant for modal solvers while a
rigid foreground slides or rotates.

```ini
streamcenterplus_meshset=1
zone_count=2
mesh_kind=sliding                 # static | nested | multigrid | sliding
zone0_source=background.pvd
zone0_motion=static
zone0_velocity_frame=world
zone1_source=rotor.pvd
zone1_motion=rigid
zone1_translation=2,2,0
zone1_rotation_axis=0,0,1
zone1_rotation_origin=0,0,0        # optional, defaults to 0,0,0
zone1_angular_velocity=62.831853   # rad/s
zone1_phase=0                      # rad
zone1_linear_velocity=0,0,0        # optional world translation rate
zone1_velocity_frame=relative      # world | zone | relative
zone1_velocity_components=u,v,w    # or zone1_velocity_array=Velocity
overlap_policy=cut_background      # or blend_50_50
```

`zoneN_source` paths are relative to the manifest when not absolute. Dynamic
sources must have the same frame count and times; a one-frame source is held
static. `sliding` requires at least one `rigid` zone. The selected background
must have no angular or linear motion; otherwise the reader fails rather than
silently changing the modal reference geometry. The default background choice
is made once from the two reference frames, never per moving frame. Conversely,
a non-zero rigid angular or linear velocity requires `mesh_kind=sliding` so
solver capability checks cannot mistake moving data for a static mesh set.
Every later background frame is checked against that reference for point
coordinates, cell types, point counts and ordered cell connectivity.
For a dynamic foreground declared with `velocity_frame=zone` or `relative`,
its source-local coordinates and ordered connectivity must likewise remain
unchanged across frames; only the declared rigid transform may move it. This
prevents a changing local mesh from being silently interpreted in the wrong
velocity reference frame. The manifest schema is closed: only the documented
top-level and `zone0_`/`zone1_` keys are accepted, so misspelled settings fail
at load time instead of falling back to defaults.

Input velocity is normalized to one point-associated `Velocity(double3)` in
world components. `world` means the supplied velocity is already absolute in
world axes; `zone` rotates its components with the rigid transform; `relative`
also adds the rigid mesh velocity. Scalar `u/v/w` and `U/V/W` are recognized,
and arbitrary scalar names require `velocity_components`; a three-component
array can be named with `velocity_array`. Original velocity candidates are
removed from canonical output so consumers cannot accidentally select a
zone-frame array. `translation` and `linear_velocity` are world vectors;
`rotation_axis` is a world-axis direction and `rotation_origin` is expressed
in the zone's input coordinates before the world translation.

Velocity is optional for scalar-only mesh sets: when neither zone resolves a
velocity, geometry and shared scalar fields are still transformed/composed and
no canonical `Velocity` is emitted. A velocity is emitted only when both zones
resolve one; a one-zone velocity is removed rather than being presented as a
complete world field. Explicit `velocity_array` or `velocity_components`
declarations remain strict and fail when their arrays are unavailable or have
the wrong component count.

`cut_background` is a fixed-reference field ownership operation, not dynamic
VTK cell deletion: foreground values replace covered background point values.
`blend_50_50` averages covered floating-point point arrays exactly 50:50.
The diagnostic point arrays `SCP_MeshOwnership` (`0` background, `1`
foreground, `2` blend) and `SCP_MeshBlendWeight` (`0`, `1`, or `0.5`) must not
be applied a second time as integration weights. `SCP_MeshBackgroundCut` marks
covered background cell centres, and point/cell `SCP_MeshInterface` marks
ownership-crossing stencils. Derivative consumers must exclude
`SCP_MeshInterface` instead of differentiating across the interface;
`SCP_MeshBackgroundCut` is a provenance diagnostic, not a derivative mask, so
foreground-owned interior samples remain differentiable. The foreground must
be physically inside the background volume
(not merely inside its bounding box), and every foreground cell must contain a
fixed background point degree of freedom. A fine cell wholly embedded in a
coarse background cell therefore fails closed; refine or pre-resample the
background rather than silently losing foreground data. The canonical
reference mesh cannot represent a foreground that exits its bounds, so that
case also fails closed. `SCP_Nek5000*` quadrature/reference/Jacobian arrays
and recognized integration or spatial-weight arrays always retain the fixed
background values; they are never cut or blended from the foreground.
Ordinary analysis PointData must be present with the same component count and
exactly one tuple per zone point in both zones. A background-only,
component-mismatched, or truncated field is removed from the canonical result
rather than silently retaining stale background values in the
foreground-covered region.
When Resample promotes CellData to PointData, reader-composed point arrays win
over same-named CellData, `SCP_MeshInterface` remains its canonical point mask,
and the categorical cell cut is conservatively promoted to a point mask.
Cell-only data and
unidentified non-velocity vector/tensor fields are not reference-frame
transformed. A named/declared cell-associated velocity is first averaged onto
its zone points before world-frame conversion and its selected raw CellData
array plus standard velocity aliases are removed, so Resample cannot overwrite
canonical point `Velocity`.
Other cell-only analysis fields remain background diagnostics rather than
canonical remapped fields. For a reader-composed two-zone input, Resample
does not promote any remaining CellData to PointData: this prevents a
background-only CellData field from reappearing in the foreground-covered
region after the reader has correctly removed its unavailable PointData
counterpart.

The implementation follows the official [Nek5000 field-file description](https://nek5000.github.io/NekDoc/problem_setup/case_files.html) for the 132-byte `#std` header, endian tag, element map and field blocks. Its numerical treatment follows Nek5000's documented [GLL/GL placement and element-local derivative semantics](https://nek5000.github.io/NekDoc/faq.html). Native reads always return canonical element-local GLL point blocks with no linear VTK cells; equal coordinates on an element interface deliberately remain separate degrees of freedom. Integration, gradients, resampling, and `Nek5000SpectralPointEvaluator` use that canonical data. `CreateNek5000SpectralDisplayProxy()` evaluates the GLL polynomial onto an unmerged uniform-reference display/interoperability proxy. `CreateNek5000SpectralOperationPipeline()` is the only native Crop/Clip/Slice/Contour pathway: it samples the canonical spectral geometry and requested point field with the same evaluator, then emits an explicitly marked sampled-operation result with its maximum reference edge length. The canonical dataset is never altered or merged. Exterior PBR normals are evaluated from spectral mapping derivatives; normals on an operation-derived surface are explicitly labelled sampled-display normals.

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

Windows target and host dependency builds use the module overlay triplet with the Visual Studio 2022 `v143` platform toolset, dynamic CRT/libraries, and Release-only packages. Linux uses the matching `x64-linux-release` overlay triplet for target and host packages. Pass `-Clean` on Windows or `--clean` on Linux to remove only this module's `.build` tree before rebuilding. Every installation transactionally replaces the complete `release` prefix after a successful staged validation, so stale dated executables cannot survive.

Installers first install into a staging prefix under the fixed module build
tree. They reject links/reparse points, compute a complete staged-content
SHA-256 identity, and require exactly one non-empty Resample executable whose
`--solver-info` fields identify Streamcenter+ `type=tools`, `solver=Resample`,
version `0.3.0`, CPU enabled and CUDA disabled. Only then is the complete
`release` directory swapped into place; post-swap identity validation failure
restores the previous release.

Inside the Streamcenter+ superproject, the installer automatically uses the superproject checkout at `.deps/vcpkg`; detached use automatically uses this module's `.deps/vcpkg`. In both cases installed packages remain private at `SCP-ResamplerIO/.deps/vcpkg_installed`.

`SCP-SolverSDK` is resolved from the fixed sibling `../SCP-SolverSDK` checkout/installed release when present. A detached `SCP-ResamplerIO` checkout instead bootstraps `https://github.com/hweifluids/SCP-SolverSDK.git` at the exact tag `v0.2.0` into `.deps/SCP-SolverSDK`; the `VERSION` contract and exact checked-out tag are validated before configure. If the repository or tag has not yet been published, installation stops with an explicit publication error. The installers expose no external component-path override.

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
options.attachNek5000QuadratureArrays = false;
auto data = streamcenterplus::io::ReadDataSet(frames.at(0), options);
auto grad = streamcenterplus::io::ComputeNek5000SpectralGradient(
    data, data->GetPointData()->GetArray("Velocity"));
```

`DataSetFrame` records physical time, native reader time, source kind and source index. File-per-frame VTK callers can continue using the path overloads. Native reads are always cell-free canonical GLL data; the master quadrature switch and the reference/Jacobian/physical selectors can suppress unused point-sized arrays. Basic GLL dimensions, canonical element IDs and representation metadata remain attached when these optional arrays are disabled. `HasSameNek5000SpectralGeometry` compares tensor dimensions, canonical IDs and coordinates before an evaluator is rebound to another time-frame field.

## Smoke example

`examples/smoke` provides no-argument Windows and Linux entrypoints. It consumes the published private `SCP-examples/hemi_resampled` dataset and independently reads all 20 frames to check geometry, `Velocity(float3)`, the valid-point mask, time ordering, PVD references, and JSON syntax. The integrated superproject uses its `examples/SCP-examples` submodule; detached use resolves a sibling checkout or an authenticated private `smoke-data-v1` Release download.

The smoke runner writes its prepared data only to the fixed ignored `.tests/smoke/SCP-ResamplerIO/data` path. The original 1,000-frame Nek5000 sequence is maintainer-only provenance and is not required on a target machine.

The no-argument CTest integration target `SCPResamplerIO.Nek5000Integration` separately exercises the generic adapter and solver on advertised source frames 490 and 491 with a small `16 x 8 x 6` grid. Its output is fixed at `.tests/integration/SCP-ResamplerIO/nek5000-general-adapter`. The test is reported as skipped when the fixed source dataset is not available; it never writes into the source or release trees.
