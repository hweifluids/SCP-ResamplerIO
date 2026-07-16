#pragma once

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <vtkSmartPointer.h>

class vtkDataSet;
class vtkImageData;

namespace streamcenterplus::resample {

enum class Backend {
    Cpu,
    Gpu,
};

enum class InterpolationMethod {
    Linear,
    InverseDistanceWeighting,
    Kriging,
    VtkResampleToImage,
};

enum class NeighborMode {
    Radius,
    Nearest,
};

enum class NullStrategy {
    Mask,
    Closest,
    NullValue,
    Error,
};

enum class KrigingVariogram {
    Spherical,
    Exponential,
    Gaussian,
    Linear,
};

enum class KrigingFallback {
    Idw,
    Closest,
    Error,
};

enum class BoundsMode {
    Global,
    PerFrame,
};

// Describes how a logical frame is obtained from its source path. Most VTK
// inputs map one logical frame to one file. A Nek5000 descriptor instead maps
// every advertised reader time step to the same descriptor path.
enum class DataSetFrameKind {
    File,
    Nek5000TimeStep,
};

struct LinearParameters {
    bool computeTolerance = true;
    double tolerance = 1e-6;
    bool snapToCellWithClosestPoint = false;
};

struct IdwParameters {
    NeighborMode neighborMode = NeighborMode::Nearest;
    int neighbors = 32;
    double radius = 0.0;
    double power = 2.0;
    NullStrategy nullStrategy = NullStrategy::Mask;
    double nullValue = 0.0;
};

struct KrigingParameters {
    KrigingVariogram variogram = KrigingVariogram::Spherical;
    NeighborMode neighborMode = NeighborMode::Nearest;
    int neighbors = 16;
    int minNeighbors = 4;
    double radius = 0.0;
    double range = 0.0;
    double sill = 1.0;
    double nugget = 0.0;
    double regularization = 1e-10;
    double nullValue = 0.0;
    KrigingFallback fallback = KrigingFallback::Idw;
    double fallbackIdwPower = 2.0;
};

struct AxisRange {
    std::optional<double> min;
    std::optional<double> max;
};

struct FrameRange {
    std::optional<int> start;
    std::optional<int> end;
};

struct Parameters {
    std::filesystem::path deckFile;
    std::filesystem::path inputPath;
    std::filesystem::path outputDir = std::filesystem::path("./results/resample_results");
    std::string solverStatus = "unsolved";
    bool enabled = true;
    Backend backend = Backend::Cpu;
    InterpolationMethod method = InterpolationMethod::Linear;
    std::array<int, 3> dimensions = {64, 64, 64};
    BoundsMode boundsMode = BoundsMode::Global;
    std::string outputBasename;
    double physicalDt = 1.0;
    FrameRange frameRange;
    std::array<AxisRange, 3> cropRanges;
    std::array<double, 3> inputScale = {1.0, 1.0, 1.0};

    LinearParameters linear;
    IdwParameters idw;
    KrigingParameters kriging;
};

struct FrameInfo {
    std::filesystem::path path;
    double time = 0.0;
    DataSetFrameKind kind = DataSetFrameKind::File;
    std::optional<double> readerTime;
};

struct SequenceInfo {
    std::vector<FrameInfo> frames;
    std::array<double, 6> bounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    bool boundsInitialized = false;
    bool seriesBoundsPrecomputed = false;
};

struct FrameResult {
    int frameIndex = 0;
    double time = 0.0;
    std::filesystem::path relativeFile;
    std::array<double, 6> bounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    long long pointCount = 0;
    long long validPointCount = 0;
    long long hiddenCellCount = 0;
    double readSeconds = 0.0;
    double prepareSeconds = 0.0;
    double resampleSeconds = 0.0;
    double writeSeconds = 0.0;
    double statsSeconds = 0.0;
    double totalSeconds = 0.0;
};

std::string ToString(Backend value);
std::string ToString(InterpolationMethod value);
std::string ToString(NeighborMode value);
std::string ToString(NullStrategy value);
std::string ToString(KrigingVariogram value);
std::string ToString(KrigingFallback value);
std::string ToString(BoundsMode value);
std::string ToString(DataSetFrameKind value);

Parameters LoadParameters(const std::filesystem::path& deckFile);
std::string PrettyPrint(const Parameters& params);

SequenceInfo CollectInputSequence(const Parameters& params);
// Enumerates every logical frame advertised by a supported dataset source.
// For .nek5000 this returns one frame per VTK reader TIME_STEPS value.
std::vector<FrameInfo> EnumerateDataSetFrames(const std::filesystem::path& sourcePath,
                                              double physicalDt = 1.0);
std::array<double, 6> ReadDataSetBounds(const std::filesystem::path& filePath);
std::array<double, 6> ReadDataSetBounds(const FrameInfo& frame);
vtkSmartPointer<vtkDataSet> ReadDataSet(const std::filesystem::path& filePath);
vtkSmartPointer<vtkDataSet> ReadDataSet(const FrameInfo& frame);
vtkSmartPointer<vtkDataSet> PreparePointData(vtkDataSet* input);
std::array<double, 6> ResolveCropBounds(const Parameters& params,
                                        const std::array<double, 6>& inputBounds);
std::array<double, 6> ResolveConfiguredBounds(const Parameters& params,
                                              const std::array<double, 6>& inputBounds);
vtkSmartPointer<vtkImageData> ResampleFrame(vtkDataSet* preparedInput,
                                            const Parameters& params,
                                            const std::array<double, 6>& bounds);
void WriteVti(const std::filesystem::path& outputPath, vtkImageData* image);
void WritePvd(const std::filesystem::path& outputPath,
              const std::vector<std::pair<double, std::filesystem::path>>& frames);
void WriteTextFile(const std::filesystem::path& outputPath, const std::string& text);

std::string VtkSmpRuntimeSummary();
int RunVtkSmpSmokeTest();
int RunSolver(const Parameters& params);
int RunCli(int argc, char** argv);

}  // namespace streamcenterplus::resample
