#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <vtkSmartPointer.h>

class vtkDataArray;
class vtkDataSet;
class vtkDoubleArray;

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

enum class DataSetFrameKind {
    File,
    Nek5000TimeStep,
    EnSightGoldTimeStep,
    OpenFoamTimeStep,
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
};

struct ReadOptions {
    bool readAllArrays = true;
    bool preserveCompositeBlockBoundaries = true;
    bool attachSourceBlockIds = true;
    bool attachNek5000SpectralMetadata = true;
    // The linear GLL-to-GLL subcells are useful for ordinary VTK consumers,
    // but are prohibitively large for high-order production meshes. Spectral
    // metadata and point arrays do not depend on these cells.
    bool buildNek5000LinearizedTopology = true;
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
    bool Evaluate(const double physicalPoint[3],
                  double* tuple,
                  int* elementHint = nullptr) const;

private:
    struct Impl;
    explicit Nek5000SpectralPointEvaluator(std::shared_ptr<const Impl> impl);
    std::shared_ptr<const Impl> impl_;
};

std::string ToString(DataSetFrameKind kind);
bool IsSupportedDataSource(const std::filesystem::path& sourcePath);

// Enumerates logical frames without converting the source to an intermediate
// VTK file. A directory is treated as an OpenFOAM case when it contains a
// *.foam marker (or a recognizable system/constant/processor layout), and as a
// normal VTK frame directory otherwise.
std::vector<DataSetFrame> EnumerateDataSetFrames(
    const std::filesystem::path& sourcePath,
    double fallbackPhysicalDt = 1.0);

vtkSmartPointer<vtkDataSet> ReadDataSet(
    const std::filesystem::path& sourcePath,
    const ReadOptions& options = {});
vtkSmartPointer<vtkDataSet> ReadDataSet(
    const DataSetFrame& frame,
    const ReadOptions& options = {});

bool HasNek5000SpectralMetadata(vtkDataSet* dataSet);
Nek5000SpectralMetadata GetNek5000SpectralMetadata(vtkDataSet* dataSet);

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
