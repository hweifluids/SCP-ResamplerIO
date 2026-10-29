#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <vtkSmartPointer.h>
#include <streamcenterplus/EnSightTimePolicy.h>

class vtkDataArray;
class vtkDataObject;
class vtkDataSet;
class vtkDoubleArray;
class vtkPolyData;

namespace streamcenterplus::io {

inline constexpr const char* kNek5000QuadratureWeightArray =
    "SCP_Nek5000QuadratureWeight";
inline constexpr const char* kNek5000ReferenceWeightArray =
    "SCP_Nek5000ReferenceWeight";
inline constexpr const char* kNek5000JacobianDeterminantArray =
    "SCP_Nek5000JacobianDeterminant";
inline constexpr const char* kNek5000SpectralElementIdArray =
    "SCP_Nek5000SpectralElementId";
inline constexpr const char* kNek5000GllDimensionsField =
    "SCP_Nek5000GLLDimensions";
inline constexpr const char* kNek5000ElementCountField =
    "SCP_Nek5000ElementCount";
inline constexpr const char* kNek5000ElementIdsField =
    "SCP_Nek5000ElementIds";
inline constexpr const char* kNek5000RepresentationField =
    "SCP_Nek5000Representation";
// Metadata and provenance on a display-only proxy created from canonical
// element-local GLL data.  These markers deliberately distinguish a sampled
// VTK representation from the canonical spectral source.
inline constexpr const char* kNek5000DisplayProxyField =
    "SCP_Nek5000DisplayProxy";
inline constexpr const char* kNek5000DisplayRepresentationField =
    "SCP_Nek5000DisplayRepresentation";
inline constexpr const char* kNek5000DisplayElementIdArray =
    "SCP_Nek5000DisplayElementId";
inline constexpr const char* kNek5000DisplayElementIndexArray =
    "SCP_Nek5000DisplayElementIndex";
inline constexpr const char* kNek5000DisplayLocalGllIndexArray =
    "SCP_Nek5000DisplayLocalGllIndex";
inline constexpr const char* kNek5000DisplayReferenceCoordinatesArray =
    "SCP_Nek5000DisplayReferenceCoordinates";
inline constexpr const char* kNek5000DisplayFaceIdArray =
    "SCP_Nek5000DisplayFaceId";
// Every display-only result produced from canonical Nek5000 data records the
// evaluator-driven operation and the largest uniform reference-space edge in
// the sampled representation.  The spacing is an approximation contract, not
// a claim of a global physical-space error bound for an arbitrary curved map.
inline constexpr const char* kNek5000DisplayOperationField =
    "SCP_Nek5000DisplayOperation";
inline constexpr const char* kNek5000DisplayReferenceSpacingField =
    "SCP_Nek5000DisplayReferenceSpacing";
inline constexpr const char* kNek5000DisplaySourceElementIdsField =
    "SCP_Nek5000DisplaySourceElementIds";
inline constexpr const char* kNek5000SpectralSurfaceNormalArray =
    "Normals";

// Point-data provenance emitted when a two-zone .scpm manifest is read.  The
// returned data set always uses the larger zone as a fixed world-reference
// mesh.  These arrays describe how each reference point was populated; they
// are diagnostics and must not be multiplied into integration weights.
inline constexpr const char* kMultiZoneCanonicalVelocityArray =
    "Velocity";
inline constexpr const char* kMultiZoneOwnershipArray =
    "SCP_MeshOwnership";
inline constexpr const char* kMultiZoneBlendWeightArray =
    "SCP_MeshBlendWeight";
inline constexpr const char* kMultiZoneBackgroundZoneArray =
    "SCP_MeshBackgroundZone";
inline constexpr const char* kMultiZoneForegroundZoneArray =
    "SCP_MeshForegroundZone";
// Cell-data diagnostic only.  A value of one says the cell centre is covered
// by the foreground zone.  It is deliberately not vtkGhostType: the fixed
// reference topology remains valid for modal analysis and resampling.
inline constexpr const char* kMultiZoneBackgroundCutArray =
    "SCP_MeshBackgroundCut";
// Point and cell arrays with this same name mark a canonical reference
// stencil that crosses ownership. Gradient-based consumers should exclude
// these locations rather than differentiating through a cut/sliding interface.
inline constexpr const char* kMultiZoneInterfaceArray =
    "SCP_MeshInterface";

enum class DataSetFrameKind {
    File,
    Nek5000TimeStep,
    EnSightGoldTimeStep,
    OpenFoamTimeStep,
    // A logical time frame assembled from exactly two data sources described
    // by a Streamcenter+ mesh-set manifest (.scpm).
    MultiZoneManifestTimeStep,
};

struct DataSetFrame {
    std::filesystem::path path;
    double time = 0.0;
    DataSetFrameKind kind = DataSetFrameKind::File;
    std::optional<double> readerTime;
    std::filesystem::path coordinatePath;
    int sourceIndex = 0;
    // For Nek5000 frames without an X block, this identifies the enumerated
    // mesh-bearing frame selected for coordinate reuse. A negative value
    // means that the source kind does not use a separate coordinate frame or
    // that a caller constructed the frame without enumeration metadata.
    int coordinateSourceIndex = -1;
    // Set only for MultiZoneManifestTimeStep frames.  It records the stable
    // source-zone choice made while enumerating the manifest, so every frame
    // in a modal sequence uses identical reference topology.
    int multiZoneBackgroundIndex = -1;
    // Present only after the complete EnSight source/PVD axis passed the
    // bounded decimal-quantization policy. Original timestamps stay in time;
    // this full-axis provenance survives contiguous or regular-stride subsets.
    std::optional<ensight_time::TimeFrameProvenance> ensightTime;
};

struct ReadOptions {
    bool readAllArrays = true;
    bool preserveCompositeBlockBoundaries = true;
    bool attachSourceBlockIds = true;
    bool attachNek5000SpectralMetadata = true;
    // Native Nek5000 output is always canonical element-local GLL data with
    // no VTK linear subcells.  Consumers that need renderable topology must
    // use CreateNek5000SpectralDisplayProxy or the operation API below.
    // Controls the three per-point reference/Jacobian/physical quadrature
    // arrays. Basic GLL dimensions, element ids and representation metadata
    // remain attached when spectral metadata is enabled.
    bool attachNek5000QuadratureArrays = true;
    // Granular selectors used when quadrature arrays are enabled. Analysis
    // modules that only need the physical integration measure can avoid two
    // additional point-sized arrays on production meshes.
    bool attachNek5000ReferenceWeightArray = true;
    bool attachNek5000JacobianDeterminantArray = true;
    bool attachNek5000QuadratureWeightArray = true;
    // The generic XML reader normally logs child-reader array failures and
    // still returns a partial data object. Mesh-set composition opts into a
    // fail-closed policy; ordinary single-grid reads retain their established
    // fast path and VTK error-handling behavior.
    bool failOnVtkXmlReadErrors = false;
    // EnSight commonly exports a volume together with separate boundary parts.
    // Analysis uses the 3-D parts when any exist; all parts remain available
    // for surface-only sources or when this option is explicitly enabled.
    bool includeEnSightBoundaryParts = false;
};

struct Nek5000SpectralMetadata {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::size_t elementCount = 0;
    std::vector<int> elementIds;

    bool is3D() const noexcept { return nz > 1; }
    std::size_t pointsPerElement() const noexcept {
        return static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) *
               static_cast<std::size_t>(nz);
    }
    std::size_t pointCount() const noexcept {
        return elementCount * pointsPerElement();
    }
};

// Controls a derived VTK representation for visualization/interoperability.
// It evaluates the original tensor-product GLL polynomial at uniformly spaced
// reference samples; it does not label non-uniform GLL nodes as VTK Lagrange
// or Bezier cells.  Every derived point keeps element/local provenance and
// interfaces are intentionally left unmerged.
struct Nek5000SpectralDisplayOptions {
    // Zero chooses a high-fidelity default: two sampled segments per GLL
    // interval, bounded by maxSamplesPerAxis.  One gives a lower-memory,
    // one-segment-per-interval proxy while still sampling at uniform reference
    // locations rather than treating raw GLL points as linear topology.
    int samplesPerGllInterval = 2;
    int maxSamplesPerAxis = 33;
    bool buildVolumeProxy = true;
    bool buildExteriorSurface = true;
    bool buildElementGrid = true;
    bool buildGllNodes = true;
};

struct Nek5000SpectralDisplayProxy {
    // Piecewise-linear, uniformly sampled visualization/interoperability
    // volume.  Generic VTK clip/cutter/contour operations performed on this
    // object are approximate display operations, never strict SEM operations.
    vtkSmartPointer<vtkDataSet> volume;
    // Only exterior spectral-element faces are emitted.  Shared faces are
    // discovered through face provenance and coordinate comparison but no
    // canonical point IDs are ever merged or modified.
    vtkSmartPointer<vtkPolyData> exteriorSurface;
    // Curved exterior-face GLL grid lines and all element-local GLL nodes for
    // optional visual overlays.
    vtkSmartPointer<vtkPolyData> elementGrid;
    vtkSmartPointer<vtkPolyData> gllNodes;
};

// Native Nek5000 visualization operations are always evaluated from the
// canonical element-local GLL source.  Their VTK output is an explicitly
// labelled, unmerged, uniform-reference sampled approximation; no operation
// uses raw GLL points as VTK Lagrange nodes or as linear GLL-to-GLL subcells.
enum class Nek5000SpectralOperationKind {
    Crop,
    Clip,
    Slice,
    Contour,
};

struct Nek5000SpectralOperationRequest {
    Nek5000SpectralOperationKind kind = Nek5000SpectralOperationKind::Crop;
    std::array<double, 6> cropBounds = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    std::array<double, 3> planeOrigin = {0.0, 0.0, 0.0};
    std::array<double, 3> planeNormal = {1.0, 0.0, 0.0};
    // The scalar is looked up in canonical point data, evaluated at every
    // display sample, and then contoured.  contourComponent >= 0 selects a
    // vector component; -1 requests the Euclidean magnitude.
    std::string contourArrayName;
    int contourComponent = 0;
    std::vector<double> contourValues;
};

struct Nek5000SpectralOperationOptions {
    // Reuses the same evaluator-backed reference sampling used by the display
    // proxy.  The actual maximum reference edge is returned in the result and
    // stored in field data for provenance.
    int samplesPerGllInterval = 2;
    int maxSamplesPerAxis = 33;
    // Builds a PBR-ready surface.  For original exterior faces its normals
    // are calculated from spectral mapping derivatives.  Normals on a
    // clipped/sliced/contoured result are explicitly marked sampled-display
    // normals because the operation surface is an approximation.
    bool preparePbrSurface = false;
};

struct Nek5000SpectralOperationResult {
    vtkSmartPointer<vtkDataSet> data;
    vtkSmartPointer<vtkPolyData> pbrSurface;
    double maximumReferenceEdgeLength = 0.0;
    bool isApproximation = true;
};

struct Nek5000SpectralProbeResult {
    int elementIndex = -1;
    int elementId = -1;
    std::array<double, 3> referenceCoordinates = {0.0, 0.0, 0.0};
    std::vector<double> values;
};

struct Nek5000PointEvaluationOptions {
    // Reference-space overshoot accepted at an element boundary. This is only
    // a containment tolerance; interpolation is still evaluated at the
    // converged reference coordinate.
    double referenceTolerance = 1.0e-8;
    // Physical convergence threshold is max(absoluteTolerance,
    // relativeTolerance * elementBoundingBoxDiagonal).
    double relativePhysicalTolerance = 1.0e-10;
    double absolutePhysicalTolerance = 1.0e-12;
    int maxNewtonIterations = 24;
};

// Evaluates an element-local tensor-product GLL polynomial without merging
// interface nodes or replacing the spectral element by linear subcells.  The
// optional elementHint is both input and output: a valid previous element is
// tried first, then the remaining elements are searched.  At a shared
// interface either adjacent element may be selected, which preserves Nek5000's
// C0 field semantics while avoiding derivative averaging.
class Nek5000SpectralPointEvaluator {
public:
    Nek5000SpectralPointEvaluator(
        vtkDataSet* dataSet,
        vtkDataArray* pointArray,
        const Nek5000PointEvaluationOptions& options = {});
    ~Nek5000SpectralPointEvaluator();

    Nek5000SpectralPointEvaluator(const Nek5000SpectralPointEvaluator&) noexcept;
    Nek5000SpectralPointEvaluator& operator=(
        const Nek5000SpectralPointEvaluator&) noexcept;
    Nek5000SpectralPointEvaluator(Nek5000SpectralPointEvaluator&&) noexcept;
    Nek5000SpectralPointEvaluator& operator=(
        Nek5000SpectralPointEvaluator&&) noexcept;

    int numberOfComponents() const noexcept;
    // Conservative physical-space bounds of every spectral element owned by
    // this evaluator, ordered as [xmin, xmax, ymin, ymax, zmin, zmax]. The
    // bounds enclose the complete high-order tensor-product mapping, not only
    // the supplied GLL coordinate nodes.
    std::array<double, 6> bounds() const noexcept;
    // Creates a field evaluator that reuses this evaluator's immutable
    // geometry, bounds and locator. The new array must use the same GLL point
    // ordering. This avoids duplicating large spectral meshes for time series.
    Nek5000SpectralPointEvaluator WithPointArray(
        vtkDataArray* pointArray) const;
    // Evaluates the element-local GLL tensor polynomial directly at a known
    // reference coordinate.  This is intended for SEM-aware derived display
    // adapters and diagnostics; it never performs a VTK-cell interpolation.
    // elementIndex is the canonical local element ordinal, not a VTK cell id.
    bool EvaluateReference(std::size_t elementIndex,
                           const double referencePoint[3],
                           double* tuple) const;
    bool Evaluate(const double physicalPoint[3],
                  double* tuple,
                  int* elementHint = nullptr) const;
    // Like Evaluate(), while also returning the element-local reference
    // position selected by the high-order inverse mapping.  This is the
    // native probe/picking path; it never interpolates a display proxy.
    bool EvaluateWithReference(const double physicalPoint[3],
                               double* tuple,
                               int* elementHint,
                               std::array<double, 3>* referencePoint) const;
    // Evaluates the spectral geometry mapping and its reference derivatives.
    // The derivative vectors are useful for evaluator-derived exterior PBR
    // normals and never rely on VTK's linear-cell derivatives.
    bool EvaluateGeometryReference(std::size_t elementIndex,
                                   const double referencePoint[3],
                                   std::array<double, 3>* physicalPoint,
                                   std::array<double, 3>* derivativeR,
                                   std::array<double, 3>* derivativeS,
                                   std::array<double, 3>* derivativeT) const;

private:
    struct Impl;
    explicit Nek5000SpectralPointEvaluator(std::shared_ptr<const Impl> impl);
    std::shared_ptr<const Impl> impl_;
};

std::string ToString(DataSetFrameKind kind);
std::optional<double> ValidatedEnSightTimeStep(const std::vector<DataSetFrame>& frames);
bool IsSupportedDataSource(const std::filesystem::path& sourcePath);
// Returns true only for a syntactically valid Streamcenter+ two-zone manifest.
// A manifest is intentionally a source, rather than a new solver input type:
// callers continue to use EnumerateDataSetFrames/ReadDataSet and receive a
// fixed-world-reference vtkDataSet.
bool IsMultiZoneMeshManifest(const std::filesystem::path& sourcePath);

// Enumerates logical frames without converting the source to an intermediate
// VTK file.  Direct sources include VTK data, Tecplot ASCII (.dat) and classic
// binary PLT (.plt) finite-element data, Nek5000 descriptors, and EnSight Gold
// descriptors (.case or Fluent .encas, with case-insensitive extensions).
// EnSight frames retain native times, including a static geometry time set
// shared by transient variable frames. A directory is treated as an OpenFOAM
// case when it contains a *.foam marker (or a recognizable
// system/constant/processor layout), and as a normal frame directory otherwise.
// Exactly one .pvd collection in a normal frame directory is authoritative;
// sibling files are not mixed into its ordered, timestamped frame list. A
// directory with multiple .pvd collections must be disambiguated by selecting
// a collection explicitly.
std::vector<DataSetFrame> EnumerateDataSetFrames(
    const std::filesystem::path& sourcePath,
    double fallbackPhysicalDt = 1.0);

vtkSmartPointer<vtkDataSet> ReadDataSet(
    const std::filesystem::path& sourcePath,
    const ReadOptions& options = {});
vtkSmartPointer<vtkDataSet> ReadDataSet(
    const DataSetFrame& frame,
    const ReadOptions& options = {});

// Shared EnSight assembly for native readers and conversion/import consumers.
// A part is volumetric if any of its cell types has dimension three. If such
// parts exist, separate lower-dimensional parts are excluded by default; cells
// of different dimensions within a retained part are kept together. Pure 2-D
// sources retain every part. Composite point boundaries and source-block IDs
// follow ReadOptions, with original nonempty part ordinals retained as IDs.
vtkSmartPointer<vtkDataSet> CollapseEnSightGoldDataObject(
    vtkDataObject* object,
    const ReadOptions& options = {});

// Returns true only when the complete Nek5000 metadata contract is present and
// self-consistent with the point count: tensor dimensions, element count and
// one unique positive global element ID per element. GetNek5000SpectralMetadata
// reports the corresponding reason when this validation fails.
bool HasNek5000SpectralMetadata(vtkDataSet* dataSet);
Nek5000SpectralMetadata GetNek5000SpectralMetadata(vtkDataSet* dataSet);

// Returns true for one of the explicit derived display representations below,
// rather than for the canonical element-local GLL source data.
bool IsNek5000SpectralDisplayProxy(vtkDataSet* dataSet);

// Creates a traceable display proxy from canonical Nek5000 data.  The
// canonical input is read only and remains the sole source for strict SEM
// operations.  The default proxy uses high-fidelity uniform-reference
// sampling; it never treats raw GLL nodes as VTK cells.
Nek5000SpectralDisplayProxy CreateNek5000SpectralDisplayProxy(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralDisplayOptions& options = {});

// Applies one evaluator-sampled Native Nek5000 operation.  For a chain of
// filters use CreateNek5000SpectralOperationPipeline so all stages share one
// unmerged SEM-derived sampling representation and provenance contract.
Nek5000SpectralOperationResult CreateNek5000SpectralOperation(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralOperationRequest& request,
    const Nek5000SpectralOperationOptions& options = {});

Nek5000SpectralOperationResult CreateNek5000SpectralOperationPipeline(
    vtkDataSet* canonicalDataSet,
    const std::vector<Nek5000SpectralOperationRequest>& requests,
    const Nek5000SpectralOperationOptions& options = {});

// Evaluates a Native Nek5000 point field directly at a physical point.  The
// output preserves the selected element-local degree-of-freedom provenance;
// on a shared interface either valid adjacent element may be reported.
bool ProbeNek5000SpectralField(
    vtkDataSet* canonicalDataSet,
    vtkDataArray* pointArray,
    const double physicalPoint[3],
    Nek5000SpectralProbeResult* result,
    const Nek5000PointEvaluationOptions& options = {});

// Returns true when two spectral datasets use the same tensor dimensions,
// canonical global-element ordering and point coordinates. tolerance is an
// absolute coordinate tolerance; the default compares coordinate values
// exactly.
bool HasSameNek5000SpectralGeometry(vtkDataSet* lhs,
                                    vtkDataSet* rhs,
                                    double tolerance = 0.0);

// Computes the element-local spectral derivative at every GLL point. The
// output contains 3 components per input component in component-major order:
// [df0/dx, df0/dy, df0/dz, df1/dx, ...]. Interface nodes remain duplicated so
// derivatives are never averaged across spectral elements.
vtkSmartPointer<vtkDoubleArray> ComputeNek5000SpectralGradient(
    vtkDataSet* dataSet,
    vtkDataArray* pointArray,
    const std::string& outputName = {});

}  // namespace streamcenterplus::io
