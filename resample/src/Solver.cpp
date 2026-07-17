#include "streamcenterplus_resample/Core.h"
#include "streamcenterplus/io/DataSetReader.h"

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkBox.h>
#include <vtkCharArray.h>
#include <vtkClipDataSet.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdList.h>
#include <vtkImageData.h>
#include <vtkMultiThreader.h>
#include <vtkNew.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkPointInterpolator.h>
#include <vtkPoints.h>
#include <vtkProbeFilter.h>
#include <vtkResampleToImage.h>
#include <vtkSMP.h>
#include <vtkSMPTools.h>
#include <vtkShepardKernel.h>
#include <vtkSmartPointer.h>
#include <vtkStaticPointLocator.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace streamcenterplus::resample {

namespace {

constexpr const char* kValidPointMaskName = "vtkValidPointMask";

struct ImageStats {
    vtkIdType pointCount = 0;
    vtkIdType validPointCount = 0;
    vtkIdType hiddenCellCount = 0;
};

struct SourceArray {
    vtkDataArray* input = nullptr;
    vtkSmartPointer<vtkDoubleArray> output;
};

struct VtkSmpConfiguration {
    int logicalProcessorCount = 1;
    int processorGroupCount = 1;
    int requestedThreadCount = 1;
    int estimatedThreadCount = 1;
    std::string backend = "Sequential";
};

int DetectLogicalProcessorCount() {
#ifdef _WIN32
    const DWORD active = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
    if (active > 0) {
        return static_cast<int>(active);
    }
#endif
    const unsigned int fallback = std::thread::hardware_concurrency();
    return std::max(1, static_cast<int>(fallback));
}

int DetectProcessorGroupCount() {
#ifdef _WIN32
    const WORD groups = GetActiveProcessorGroupCount();
    return std::max(1, static_cast<int>(groups));
#else
    return 1;
#endif
}

unsigned int CurrentProcessorGroup() {
#ifdef _WIN32
    PROCESSOR_NUMBER processorNumber = {};
    GetCurrentProcessorNumberEx(&processorNumber);
    return static_cast<unsigned int>(processorNumber.Group);
#else
    return 0;
#endif
}

std::string LowerBackendName(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool IsSequentialBackend(const std::string& backend) {
    return LowerBackendName(backend) == "sequential";
}

std::vector<const char*> AvailableParallelBackends() {
    std::vector<const char*> backends;
#if defined(VTK_SMP_ENABLE_TBB) && VTK_SMP_ENABLE_TBB
    backends.push_back("TBB");
#endif
#if defined(VTK_SMP_ENABLE_STDTHREAD) && VTK_SMP_ENABLE_STDTHREAD
    backends.push_back("STDThread");
#endif
#if defined(VTK_SMP_ENABLE_OPENMP) && VTK_SMP_ENABLE_OPENMP
    backends.push_back("OpenMP");
#endif
    return backends;
}

VtkSmpConfiguration ConfigureVtkSmpForMaxConcurrency() {
    VtkSmpConfiguration config;
    config.logicalProcessorCount = DetectLogicalProcessorCount();
    config.processorGroupCount = DetectProcessorGroupCount();
    config.requestedThreadCount = std::max(1, config.logicalProcessorCount);

    std::string backend = vtkSMPTools::GetBackend() != nullptr ? vtkSMPTools::GetBackend() : "Sequential";
    if (IsSequentialBackend(backend) && config.requestedThreadCount > 1) {
        bool selected = false;
        for (const char* candidate : AvailableParallelBackends()) {
            if (vtkSMPTools::SetBackend(candidate)) {
                selected = true;
                break;
            }
        }
        if (!selected) {
            throw std::runtime_error(
                "VTK was built without an available parallel SMP backend. Rebuild VTK with "
                "VTK_SMP_IMPLEMENTATION_TYPE=STDThread, TBB, or OpenMP.");
        }
    }

    vtkSMPTools::Initialize(config.requestedThreadCount);
    vtkMultiThreader::SetGlobalMaximumNumberOfThreads(config.requestedThreadCount);
    vtkMultiThreader::SetGlobalDefaultNumberOfThreads(config.requestedThreadCount);

    config.backend = vtkSMPTools::GetBackend() != nullptr ? vtkSMPTools::GetBackend() : "unknown";
    config.estimatedThreadCount = std::max(1, vtkSMPTools::GetEstimatedNumberOfThreads());
    if (config.requestedThreadCount > 1 && IsSequentialBackend(config.backend)) {
        throw std::runtime_error(
            "VTK SMP is still running serially after initialization. Rebuild VTK with a parallel SMP backend.");
    }
    if (config.requestedThreadCount > 1 && config.estimatedThreadCount < config.requestedThreadCount) {
        std::ostringstream message;
        message << "VTK SMP backend " << config.backend << " exposes only "
                << config.estimatedThreadCount << " worker threads after requesting "
                << config.requestedThreadCount << ". Rebuild VTK with a backend that can use all active "
                << "logical processors, such as TBB.";
        throw std::runtime_error(message.str());
    }
    return config;
}

const VtkSmpConfiguration& EnsureVtkSmpConfigured() {
    static std::once_flag once;
    static VtkSmpConfiguration config;
    std::call_once(once, []() {
        config = ConfigureVtkSmpForMaxConcurrency();
    });
    return config;
}

std::string FormatSmpSummary(const VtkSmpConfiguration& config) {
    std::ostringstream out;
    out << "backend=" << config.backend
        << ", logical_processors=" << config.logicalProcessorCount
        << ", processor_groups=" << config.processorGroupCount
        << ", requested_threads=" << config.requestedThreadCount
        << ", estimated_threads=" << config.estimatedThreadCount;
    return out.str();
}

std::string FormatSeconds(double seconds) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(2) << seconds;
    return stream.str();
}

std::string FormatBounds(const std::array<double, 6>& bounds) {
    std::ostringstream out;
    out << std::setprecision(8)
        << "[" << bounds[0] << ", " << bounds[1] << "] x "
        << "[" << bounds[2] << ", " << bounds[3] << "] x "
        << "[" << bounds[4] << ", " << bounds[5] << "]";
    return out.str();
}

std::chrono::steady_clock::time_point LogStageStarted(const std::string& label,
                                                      const std::string& detail = std::string()) {
    const auto start = std::chrono::steady_clock::now();
    std::cout << "[Stage] " << label << " started";
    if (!detail.empty()) {
        std::cout << ": " << detail;
    }
    std::cout << "\n";
    return start;
}

void LogStageFinished(const std::string& label,
                      const std::chrono::steady_clock::time_point& start,
                      const std::string& detail = std::string()) {
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::cout << "[Stage] " << label << " finished in " << FormatSeconds(elapsed) << "s";
    if (!detail.empty()) {
        std::cout << ": " << detail;
    }
    std::cout << "\n";
}

bool ShouldLogProgress(vtkIdType stepIndex, vtkIdType totalSteps) {
    if (totalSteps <= 0) {
        return stepIndex == 0;
    }
    const vtkIdType current = stepIndex + 1;
    if (current == 1 || current == totalSteps) {
        return true;
    }
    const vtkIdType interval = std::max<vtkIdType>(1, totalSteps / 100);
    return current % interval == 0;
}

void LogProgress(const std::string& label,
                 vtkIdType stepIndex,
                 vtkIdType totalSteps,
                 const std::chrono::steady_clock::time_point& startTime,
                 const std::string& detail = std::string()) {
    const vtkIdType current = stepIndex + 1;
    const double percent = totalSteps <= 0
                               ? 100.0
                               : (100.0 * static_cast<double>(current) / static_cast<double>(totalSteps));
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    std::ostringstream stream;
    stream << "[Progress] " << label << ": step " << current << "/" << totalSteps
           << " (" << std::fixed << std::setprecision(1) << percent << "%)";
    if (!detail.empty()) {
        stream << ", " << detail;
    }
    stream << ", elapsed=" << std::setprecision(2) << elapsed << "s";
    std::cout << stream.str() << "\n";
}

void LogCompletedProgress(const std::string& label,
                          vtkIdType completedSteps,
                          vtkIdType totalSteps,
                          const std::chrono::steady_clock::time_point& startTime,
                          const std::string& detail = std::string()) {
    const double percent = totalSteps <= 0
                               ? 100.0
                               : (100.0 * static_cast<double>(completedSteps) / static_cast<double>(totalSteps));
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
    std::ostringstream stream;
    stream << "[Progress] " << label << ": step " << completedSteps << "/" << totalSteps
           << " (" << std::fixed << std::setprecision(1) << percent << "%)";
    if (!detail.empty()) {
        stream << ", " << detail;
    }
    stream << ", elapsed=" << std::setprecision(2) << elapsed << "s";
    std::cout << stream.str() << "\n";
}

void LogSnapshotProgress(const std::string& label,
                         vtkIdType solvedSnapshots,
                         vtkIdType totalSnapshots,
                         const std::chrono::steady_clock::time_point& startTime) {
    std::ostringstream stream;
    stream << "solved snapshot " << solvedSnapshots << "/" << totalSnapshots << " snapshot";
    LogCompletedProgress(label, solvedSnapshots, totalSnapshots, startTime, stream.str());
}

std::string LowerAscii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

int RemoveStaleVtiAndPvdFiles(const std::filesystem::path& outputDir) {
    if (!std::filesystem::exists(outputDir)) {
        return 0;
    }
    if (!std::filesystem::is_directory(outputDir)) {
        throw std::runtime_error("Output path exists but is not a directory: " + outputDir.string());
    }

    int removed = 0;
    for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(outputDir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string extension = LowerAscii(entry.path().extension().string());
        if (extension != ".vti" && extension != ".pvd") {
            continue;
        }
        std::filesystem::remove(entry.path());
        ++removed;
    }
    return removed;
}

std::array<double, 6> BoundsOf(vtkDataSet* dataSet) {
    if (dataSet == nullptr) {
        throw std::runtime_error("Internal error: null dataset when reading bounds.");
    }
    std::array<double, 6> bounds = {};
    dataSet->GetBounds(bounds.data());
    for (double value : bounds) {
        if (!std::isfinite(value)) {
            throw std::runtime_error("Input dataset has non-finite bounds.");
        }
    }
    for (int axis = 0; axis < 3; ++axis) {
        const int lo = axis * 2;
        const int hi = lo + 1;
        if (bounds[hi] > bounds[lo]) {
            continue;
        }
        const double center = bounds[lo];
        const double pad = std::max(1.0, std::abs(center) * 1e-6);
        bounds[lo] = center - 0.5 * pad;
        bounds[hi] = center + 0.5 * pad;
    }
    return bounds;
}

void MergeBounds(const std::array<double, 6>& frameBounds,
                 std::array<double, 6>* bounds,
                 bool* initialized) {
    if (bounds == nullptr || initialized == nullptr) {
        return;
    }
    if (!*initialized) {
        *bounds = frameBounds;
        *initialized = true;
        return;
    }
    (*bounds)[0] = std::min((*bounds)[0], frameBounds[0]);
    (*bounds)[1] = std::max((*bounds)[1], frameBounds[1]);
    (*bounds)[2] = std::min((*bounds)[2], frameBounds[2]);
    (*bounds)[3] = std::max((*bounds)[3], frameBounds[3]);
    (*bounds)[4] = std::min((*bounds)[4], frameBounds[4]);
    (*bounds)[5] = std::max((*bounds)[5], frameBounds[5]);
}

bool HasCropRange(const Parameters& params) {
    for (const AxisRange& range : params.cropRanges) {
        if (range.min.has_value() || range.max.has_value()) {
            return true;
        }
    }
    return false;
}

bool HasInputScale(const Parameters& params) {
    constexpr double kTolerance = 1e-14;
    return std::abs(params.inputScale[0] - 1.0) > kTolerance ||
           std::abs(params.inputScale[1] - 1.0) > kTolerance ||
           std::abs(params.inputScale[2] - 1.0) > kTolerance;
}

void ExpandDegenerateBounds(std::array<double, 6>* bounds) {
    if (bounds == nullptr) {
        return;
    }
    for (int axis = 0; axis < 3; ++axis) {
        const int lo = axis * 2;
        const int hi = lo + 1;
        if ((*bounds)[hi] > (*bounds)[lo]) {
            continue;
        }
        const double center = (*bounds)[lo];
        const double pad = std::max(1.0, std::abs(center) * 1e-6);
        (*bounds)[lo] = center - 0.5 * pad;
        (*bounds)[hi] = center + 0.5 * pad;
    }
}

double EstimateSourcePointSpacing(vtkDataSet* input) {
    if (input == nullptr || input->GetNumberOfPoints() <= 0) {
        return 0.0;
    }
    double dataBounds[6] = {};
    input->GetBounds(dataBounds);
    std::array<double, 3> extents = {
        std::max(0.0, dataBounds[1] - dataBounds[0]),
        std::max(0.0, dataBounds[3] - dataBounds[2]),
        std::max(0.0, dataBounds[5] - dataBounds[4]),
    };
    const double maxExtent = std::max({extents[0], extents[1], extents[2]});
    if (!(maxExtent > 0.0)) {
        return 0.0;
    }

    int activeDimensions = 0;
    double measure = 1.0;
    const double activeThreshold = maxExtent * 1.0e-2;
    for (double extent : extents) {
        if (extent > activeThreshold) {
            measure *= extent;
            ++activeDimensions;
        }
    }
    if (activeDimensions == 0 || !(measure > 0.0)) {
        return 0.0;
    }
    return std::pow(measure / static_cast<double>(input->GetNumberOfPoints()),
                    1.0 / static_cast<double>(activeDimensions));
}

std::array<double, 3> CropInterpolationHalo(vtkDataSet* input,
                                            const Parameters& params,
                                            const std::array<double, 6>& bounds) {
    const double sourceSpacing = EstimateSourcePointSpacing(input);
    std::array<double, 3> halo = {0.0, 0.0, 0.0};
    for (int axis = 0; axis < 3; ++axis) {
        const AxisRange& range = params.cropRanges[static_cast<std::size_t>(axis)];
        if (!range.min.has_value() && !range.max.has_value()) {
            continue;
        }
        const int dim = params.dimensions[static_cast<std::size_t>(axis)];
        if (dim <= 1) {
            continue;
        }
        const int lo = axis * 2;
        const int hi = lo + 1;
        const double extent = std::max(0.0, bounds[hi] - bounds[lo]);
        const double outputSpacing = extent > 0.0
                                         ? extent / static_cast<double>(std::max(1, dim - 1))
                                         : 0.0;
        const double numericPad = std::max({1.0, std::abs(bounds[lo]), std::abs(bounds[hi]), extent}) * 1.0e-9;
        double axisHalo = std::max({numericPad, outputSpacing, 4.0 * sourceSpacing});
        if (extent > 0.0) {
            axisHalo = std::min(axisHalo, extent * 0.1);
        }
        halo[static_cast<std::size_t>(axis)] = axisHalo;
    }
    return halo;
}

void ExpandClipBoundsForCropInterpolation(vtkDataSet* input,
                                           const Parameters& params,
                                           std::array<double, 6>* bounds) {
    if (bounds == nullptr) {
        return;
    }
    ExpandDegenerateBounds(bounds);
    const std::array<double, 3> halo = CropInterpolationHalo(input, params, *bounds);
    for (int axis = 0; axis < 3; ++axis) {
        const int lo = axis * 2;
        const int hi = lo + 1;
        const double scale = std::max({1.0,
                                       std::abs((*bounds)[lo]),
                                       std::abs((*bounds)[hi]),
                                       std::abs((*bounds)[hi] - (*bounds)[lo])});
        const double pad = std::max(scale * 1.0e-9, halo[static_cast<std::size_t>(axis)]);
        (*bounds)[lo] -= pad;
        (*bounds)[hi] += pad;
    }
}

vtkSmartPointer<vtkDataSet> ClipDataSetToBounds(vtkDataSet* input,
                                                const Parameters& params,
                                                const std::array<double, 6>& bounds) {
    if (input == nullptr) {
        throw std::runtime_error("Internal error: null dataset before resample crop.");
    }
    std::array<double, 6> clipBounds = bounds;
    // Keep a source-data halo around crop seams. The output grid still uses the
    // exact crop bounds, but interpolation must not sample from a freshly cut
    // boundary as if it were outside the domain.
    ExpandClipBoundsForCropInterpolation(input, params, &clipBounds);

    auto box = vtkSmartPointer<vtkBox>::New();
    box->SetBounds(clipBounds.data());

    auto clipper = vtkSmartPointer<vtkClipDataSet>::New();
    clipper->SetInputData(input);
    clipper->SetClipFunction(box);
    clipper->SetValue(0.0);
    clipper->InsideOutOn();
    clipper->Update();

    vtkUnstructuredGrid* rawOutput = clipper->GetOutput();
    if (rawOutput == nullptr) {
        throw std::runtime_error("Resample crop did not produce a dataset.");
    }
    if (rawOutput->GetNumberOfPoints() <= 0) {
        throw std::runtime_error("Resample crop produced no points; check x/y/z crop ranges.");
    }

    auto output = vtkSmartPointer<vtkUnstructuredGrid>::New();
    output->DeepCopy(rawOutput);
    return output;
}

vtkSmartPointer<vtkDataSet> ScalePointCoordinates(vtkDataSet* input,
                                                  const std::array<double, 3>& scale) {
    if (input == nullptr) {
        throw std::runtime_error("Internal error: null dataset before resample scale.");
    }
    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(input);
    if (pointSet == nullptr) {
        throw std::runtime_error("Internal error: resample input must be a vtkPointSet before coordinate scaling.");
    }

    vtkDataSet* rawCopy = vtkDataSet::SafeDownCast(input->NewInstance());
    if (rawCopy == nullptr) {
        throw std::runtime_error("Internal error: cannot allocate scaled resample dataset.");
    }
    vtkSmartPointer<vtkDataSet> output;
    output.TakeReference(rawCopy);
    output->DeepCopy(input);

    vtkPointSet* outputPointSet = vtkPointSet::SafeDownCast(output);
    if (outputPointSet == nullptr) {
        throw std::runtime_error("Internal error: scaled resample dataset is not a vtkPointSet.");
    }

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    const vtkIdType pointCount = input->GetNumberOfPoints();
    points->SetNumberOfPoints(pointCount);
    for (vtkIdType pointId = 0; pointId < pointCount; ++pointId) {
        double point[3] = {};
        input->GetPoint(pointId, point);
        point[0] *= scale[0];
        point[1] *= scale[1];
        point[2] *= scale[2];
        points->SetPoint(pointId, point);
    }
    outputPointSet->SetPoints(points);
    return output;
}

vtkSmartPointer<vtkDataSet> ApplyInputDomainTransform(vtkDataSet* preparedInput,
                                                      const Parameters& params,
                                                      const std::array<double, 6>& inputBounds) {
    if (preparedInput == nullptr) {
        throw std::runtime_error("Internal error: null prepared input.");
    }

    const bool cropInput = HasCropRange(params);
    const bool scaleInput = HasInputScale(params);
    if (!cropInput && !scaleInput) {
        return preparedInput;
    }

    vtkSmartPointer<vtkDataSet> working = preparedInput;
    const bool needsPointSetForScale = scaleInput && vtkPointSet::SafeDownCast(working) == nullptr;
    const bool spectralNek =
        streamcenterplus::io::HasNek5000SpectralMetadata(working);
    const bool shouldClipInput =
        (!spectralNek && params.method != InterpolationMethod::VtkResampleToImage &&
         (cropInput || vtkPointSet::SafeDownCast(working) == nullptr)) ||
        (!spectralNek && needsPointSetForScale);
    if (shouldClipInput) {
        const std::array<double, 6> cropBounds = ResolveCropBounds(params, inputBounds);
        working = ClipDataSetToBounds(working, params, cropBounds);
    }
    if (scaleInput) {
        working = ScalePointCoordinates(working, params.inputScale);
    }
    return working;
}

vtkSmartPointer<vtkImageData> CreateTargetImage(const Parameters& params,
                                                const std::array<double, 6>& bounds) {
    auto image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions(params.dimensions[0], params.dimensions[1], params.dimensions[2]);
    image->SetOrigin(bounds[0], bounds[2], bounds[4]);

    double spacing[3] = {1.0, 1.0, 1.0};
    for (int axis = 0; axis < 3; ++axis) {
        const int dim = params.dimensions[axis];
        const int lo = axis * 2;
        const int hi = lo + 1;
        spacing[axis] = dim > 1 ? (bounds[hi] - bounds[lo]) / static_cast<double>(dim - 1) : 1.0;
    }
    image->SetSpacing(spacing);
    image->AllocateScalars(VTK_DOUBLE, 1);
    image->GetPointData()->Initialize();
    return image;
}

std::vector<unsigned char> ExtractMask(vtkImageData* image) {
    std::vector<unsigned char> mask;
    if (image == nullptr) {
        return mask;
    }
    const vtkIdType pointCount = image->GetNumberOfPoints();
    mask.assign(static_cast<std::size_t>(pointCount), 1);
    vtkCharArray* charMask = vtkCharArray::SafeDownCast(image->GetPointData()->GetArray(kValidPointMaskName));
    if (charMask != nullptr && charMask->GetNumberOfTuples() == pointCount) {
        for (vtkIdType i = 0; i < pointCount; ++i) {
            mask[static_cast<std::size_t>(i)] = charMask->GetValue(i) != 0 ? 1 : 0;
        }
        return mask;
    }
    vtkDataArray* numericMask = vtkDataArray::SafeDownCast(image->GetPointData()->GetArray(kValidPointMaskName));
    if (numericMask != nullptr && numericMask->GetNumberOfTuples() == pointCount) {
        for (vtkIdType i = 0; i < pointCount; ++i) {
            mask[static_cast<std::size_t>(i)] = numericMask->GetComponent(i, 0) != 0.0 ? 1 : 0;
        }
    }
    return mask;
}

std::vector<unsigned char> ProbeDomainMask(vtkImageData* target, vtkDataSet* source) {
    auto prober = vtkSmartPointer<vtkProbeFilter>::New();
    prober->SetInputData(target);
    prober->SetSourceData(source);
    prober->PassCellArraysOff();
    prober->PassPointArraysOff();
    prober->PassFieldArraysOff();
    prober->SetValidPointMaskArrayName(kValidPointMaskName);
    prober->Update();
    vtkImageData* probed = vtkImageData::SafeDownCast(prober->GetOutput());
    if (probed == nullptr) {
        throw std::runtime_error("Internal error: vtkProbeFilter did not produce vtkImageData.");
    }
    return ExtractMask(probed);
}

vtkIdType PointIdFromIjk(int i, int j, int k, const int dims[3]) {
    return static_cast<vtkIdType>(i) +
           static_cast<vtkIdType>(dims[0]) * static_cast<vtkIdType>(j) +
           static_cast<vtkIdType>(dims[0]) * static_cast<vtkIdType>(dims[1]) * static_cast<vtkIdType>(k);
}

ImageStats MarkGhostArrays(vtkImageData* image, const std::vector<unsigned char>& validMask) {
    if (image == nullptr) {
        throw std::runtime_error("Internal error: null image when marking ghost arrays.");
    }
    const vtkIdType pointCount = image->GetNumberOfPoints();
    if (static_cast<vtkIdType>(validMask.size()) != pointCount) {
        throw std::runtime_error("Internal error: valid mask size does not match image point count.");
    }

    image->GetPointData()->RemoveArray(kValidPointMaskName);
    auto pointMask = vtkSmartPointer<vtkCharArray>::New();
    pointMask->SetName(kValidPointMaskName);
    pointMask->SetNumberOfComponents(1);
    pointMask->SetNumberOfTuples(pointCount);
    for (vtkIdType pointId = 0; pointId < pointCount; ++pointId) {
        pointMask->SetValue(
            pointId,
            validMask[static_cast<std::size_t>(pointId)] != 0 ? 1 : 0);
    }
    image->GetPointData()->AddArray(pointMask);

    image->GetPointData()->RemoveArray(vtkDataSetAttributes::GhostArrayName());
    auto pointGhost = vtkSmartPointer<vtkUnsignedCharArray>::New();
    pointGhost->SetName(vtkDataSetAttributes::GhostArrayName());
    pointGhost->SetNumberOfComponents(1);
    pointGhost->SetNumberOfTuples(pointCount);
    pointGhost->FillValue(0);

    ImageStats stats;
    stats.pointCount = pointCount;
    for (vtkIdType pointId = 0; pointId < pointCount; ++pointId) {
        if (validMask[static_cast<std::size_t>(pointId)] == 0) {
            pointGhost->SetValue(pointId, vtkDataSetAttributes::HIDDENPOINT);
        } else {
            ++stats.validPointCount;
        }
    }
    image->GetPointData()->AddArray(pointGhost);

    image->GetCellData()->RemoveArray(vtkDataSetAttributes::GhostArrayName());
    const vtkIdType cellCount = image->GetNumberOfCells();
    auto cellGhost = vtkSmartPointer<vtkUnsignedCharArray>::New();
    cellGhost->SetName(vtkDataSetAttributes::GhostArrayName());
    cellGhost->SetNumberOfComponents(1);
    cellGhost->SetNumberOfTuples(cellCount);
    cellGhost->FillValue(0);

    int pointDims[3] = {};
    image->GetDimensions(pointDims);
    const int cellDims[3] = {
        std::max(1, pointDims[0] - 1),
        std::max(1, pointDims[1] - 1),
        std::max(1, pointDims[2] - 1),
    };
    const int iterateDim[3] = {
        pointDims[0] > 1 ? 1 : 0,
        pointDims[1] > 1 ? 1 : 0,
        pointDims[2] > 1 ? 1 : 0,
    };

    vtkIdType cellId = 0;
    for (int k = 0; k < cellDims[2]; ++k) {
        for (int j = 0; j < cellDims[1]; ++j) {
            for (int i = 0; i < cellDims[0]; ++i, ++cellId) {
                bool validCell = true;
                for (int dz = 0; dz <= iterateDim[2]; ++dz) {
                    for (int dy = 0; dy <= iterateDim[1]; ++dy) {
                        for (int dx = 0; dx <= iterateDim[0]; ++dx) {
                            const vtkIdType pointId = PointIdFromIjk(i + dx, j + dy, k + dz, pointDims);
                            validCell = validCell && validMask[static_cast<std::size_t>(pointId)] != 0;
                        }
                    }
                }
                if (!validCell) {
                    cellGhost->SetValue(cellId, vtkDataSetAttributes::HIDDENCELL);
                    ++stats.hiddenCellCount;
                }
            }
        }
    }
    image->GetCellData()->AddArray(cellGhost);
    return stats;
}

std::vector<unsigned char> CombineMasks(const std::vector<unsigned char>& lhs,
                                        const std::vector<unsigned char>& rhs) {
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("Internal error: cannot combine masks with different sizes.");
    }
    std::vector<unsigned char> combined(lhs.size(), 0);
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        combined[i] = lhs[i] != 0 && rhs[i] != 0 ? 1 : 0;
    }
    return combined;
}

vtkSmartPointer<vtkImageData> DeepCopyImage(vtkImageData* image) {
    auto copy = vtkSmartPointer<vtkImageData>::New();
    copy->DeepCopy(image);
    return copy;
}

vtkSmartPointer<vtkImageData> ResampleVtkToImage(vtkDataSet* preparedInput,
                                                 const Parameters& params,
                                                 const std::array<double, 6>& bounds) {
    vtkSmartPointer<vtkDataSet> source = preparedInput;
    if (source == nullptr) {
        throw std::runtime_error("Internal error: null input before VTK ResampleToImage.");
    }

    auto resampler = vtkSmartPointer<vtkResampleToImage>::New();
    resampler->SetInputDataObject(source);
    resampler->UseInputBoundsOff();
    resampler->SetSamplingBounds(bounds.data());
    resampler->SetSamplingDimensions(params.dimensions[0], params.dimensions[1], params.dimensions[2]);
    resampler->Update();

    vtkImageData* rawOutput = resampler->GetOutput();
    if (rawOutput == nullptr) {
        throw std::runtime_error("vtkResampleToImage did not produce vtkImageData.");
    }

    auto output = vtkSmartPointer<vtkImageData>::New();
    output->ShallowCopy(rawOutput);
    output->GetFieldData()->PassData(source->GetFieldData());

    if (output->GetPointGhostArray() == nullptr || output->GetCellGhostArray() == nullptr) {
        std::vector<unsigned char> mask = ExtractMask(output);
        if (mask.empty()) {
            mask = ProbeDomainMask(output, source);
        }
        MarkGhostArrays(output, mask);
    }
    return output;
}

vtkSmartPointer<vtkImageData> ResampleLinear(vtkDataSet* preparedInput,
                                             const Parameters& params,
                                             const std::array<double, 6>& bounds) {
    vtkSmartPointer<vtkImageData> target = CreateTargetImage(params, bounds);
    auto prober = vtkSmartPointer<vtkProbeFilter>::New();
    prober->SetInputData(target);
    prober->SetSourceData(preparedInput);
    prober->SetValidPointMaskArrayName(kValidPointMaskName);
    prober->SetComputeTolerance(params.linear.computeTolerance);
    if (!params.linear.computeTolerance) {
        prober->SetTolerance(params.linear.tolerance);
    }
    prober->SetSnapToCellWithClosestPoint(params.linear.snapToCellWithClosestPoint);
    prober->PassCellArraysOff();
    prober->PassPointArraysOff();
    prober->PassFieldArraysOff();
    prober->Update();

    vtkImageData* rawOutput = vtkImageData::SafeDownCast(prober->GetOutput());
    if (rawOutput == nullptr) {
        throw std::runtime_error("vtkProbeFilter did not produce an image.");
    }
    vtkSmartPointer<vtkImageData> output = DeepCopyImage(rawOutput);
    output->GetFieldData()->PassData(preparedInput->GetFieldData());

    std::vector<unsigned char> mask = ExtractMask(output);
    if (mask.empty()) {
        mask = ProbeDomainMask(output, preparedInput);
    }
    MarkGhostArrays(output, mask);
    return output;
}

void ConfigureKernel(vtkShepardKernel* kernel, const IdwParameters& idw) {
    if (kernel == nullptr) {
        return;
    }
    kernel->SetPowerParameter(idw.power);
    if (idw.neighborMode == NeighborMode::Radius && idw.radius > 0.0) {
        kernel->SetKernelFootprintToRadius();
        kernel->SetRadius(idw.radius);
    } else {
        kernel->SetKernelFootprintToNClosest();
        kernel->SetNumberOfPoints(std::max(1, idw.neighbors));
    }
}

vtkSmartPointer<vtkImageData> ResampleIdw(vtkDataSet* preparedInput,
                                          const Parameters& params,
                                          const std::array<double, 6>& bounds) {
    vtkSmartPointer<vtkImageData> target = CreateTargetImage(params, bounds);
    const std::vector<unsigned char> domainMask = ProbeDomainMask(target, preparedInput);

    auto kernel = vtkSmartPointer<vtkShepardKernel>::New();
    ConfigureKernel(kernel, params.idw);

    auto interpolator = vtkSmartPointer<vtkPointInterpolator>::New();
    interpolator->SetInputData(target);
    interpolator->SetSourceData(preparedInput);
    interpolator->SetKernel(kernel);
    interpolator->SetValidPointsMaskArrayName(kValidPointMaskName);
    interpolator->SetNullValue(params.idw.nullValue);
    interpolator->PromoteOutputArraysOn();
    interpolator->PassCellArraysOff();
    interpolator->PassPointArraysOff();
    interpolator->PassFieldArraysOff();
    switch (params.idw.nullStrategy) {
    case NullStrategy::Mask:
    case NullStrategy::Error:
        interpolator->SetNullPointsStrategyToMaskPoints();
        break;
    case NullStrategy::Closest:
        interpolator->SetNullPointsStrategyToClosestPoint();
        break;
    case NullStrategy::NullValue:
        interpolator->SetNullPointsStrategyToNullValue();
        break;
    }
    interpolator->Update();

    vtkImageData* rawOutput = vtkImageData::SafeDownCast(interpolator->GetOutput());
    if (rawOutput == nullptr) {
        throw std::runtime_error("vtkPointInterpolator did not produce vtkImageData.");
    }
    vtkSmartPointer<vtkImageData> output = DeepCopyImage(rawOutput);
    output->GetFieldData()->PassData(preparedInput->GetFieldData());

    std::vector<unsigned char> finalMask = domainMask;
    if (params.idw.nullStrategy == NullStrategy::Mask || params.idw.nullStrategy == NullStrategy::Error) {
        const std::vector<unsigned char> interpolationMask = ExtractMask(output);
        finalMask = CombineMasks(domainMask, interpolationMask);
        if (params.idw.nullStrategy == NullStrategy::Error) {
            for (std::size_t i = 0; i < domainMask.size(); ++i) {
                if (domainMask[i] != 0 && interpolationMask[i] == 0) {
                    throw std::runtime_error("IDW interpolation has null points inside the source domain.");
                }
            }
        }
    }
    MarkGhostArrays(output, finalMask);
    return output;
}

bool IsInterpolatableArray(vtkDataArray* array) {
    if (array == nullptr || array->GetName() == nullptr) {
        return false;
    }
    const std::string name = array->GetName();
    return name != vtkDataSetAttributes::GhostArrayName() &&
           name != kValidPointMaskName &&
           name.rfind("SCP_Nek5000", 0) != 0;
}

std::vector<SourceArray> PrepareOutputArrays(vtkPointData* sourcePointData,
                                             vtkPointData* outputPointData,
                                             vtkIdType outputPointCount,
                                             double nullValue) {
    std::vector<SourceArray> arrays;
    if (sourcePointData == nullptr || outputPointData == nullptr) {
        return arrays;
    }
    for (int i = 0; i < sourcePointData->GetNumberOfArrays(); ++i) {
        vtkDataArray* inputArray = vtkDataArray::SafeDownCast(sourcePointData->GetArray(i));
        if (!IsInterpolatableArray(inputArray)) {
            continue;
        }
        auto outputArray = vtkSmartPointer<vtkDoubleArray>::New();
        outputArray->SetName(inputArray->GetName());
        outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
        outputArray->SetNumberOfTuples(outputPointCount);
        for (int component = 0; component < outputArray->GetNumberOfComponents(); ++component) {
            outputArray->FillComponent(component, nullValue);
        }
        outputPointData->AddArray(outputArray);
        if (inputArray == sourcePointData->GetScalars()) {
            outputPointData->SetScalars(outputArray);
        }
        arrays.push_back({inputArray, outputArray});
    }
    return arrays;
}

void CopyNonSpectralFieldData(vtkDataSet* source, vtkDataSet* target) {
    if (source == nullptr || target == nullptr || source->GetFieldData() == nullptr ||
        target->GetFieldData() == nullptr) {
        return;
    }
    target->GetFieldData()->Initialize();
    vtkFieldData* fieldData = source->GetFieldData();
    for (int i = 0; i < fieldData->GetNumberOfArrays(); ++i) {
        vtkAbstractArray* array = fieldData->GetAbstractArray(i);
        const std::string name =
            array == nullptr || array->GetName() == nullptr
                ? std::string()
                : std::string(array->GetName());
        if (array != nullptr && name.rfind("SCP_Nek5000", 0) != 0) {
            target->GetFieldData()->AddArray(array);
        }
    }
}

vtkSmartPointer<vtkImageData> ResampleNek5000Spectral(
    vtkDataSet* preparedInput,
    const Parameters& params,
    const std::array<double, 6>& bounds) {
    vtkSmartPointer<vtkImageData> output = CreateTargetImage(params, bounds);
    CopyNonSpectralFieldData(preparedInput, output);

    const vtkIdType outputPointCount = output->GetNumberOfPoints();
    std::vector<SourceArray> arrays =
        PrepareOutputArrays(preparedInput->GetPointData(), output->GetPointData(),
                            outputPointCount, 0.0);
    std::size_t totalComponents = 0;
    for (const SourceArray& array : arrays) {
        totalComponents +=
            static_cast<std::size_t>(array.input->GetNumberOfComponents());
    }
    if (totalComponents > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(
            "Nek5000 source has too many point-array components to interpolate.");
    }

    std::vector<unsigned char> validMask(
        static_cast<std::size_t>(outputPointCount), 0);
    if (totalComponents == 0) {
        MarkGhostArrays(output, validMask);
        return output;
    }

    vtkNew<vtkDoubleArray> packed;
    packed->SetName("SCP_Nek5000PackedInterpolationInput");
    packed->SetNumberOfComponents(static_cast<int>(totalComponents));
    packed->SetNumberOfTuples(preparedInput->GetNumberOfPoints());
    for (vtkIdType pointId = 0; pointId < preparedInput->GetNumberOfPoints();
         ++pointId) {
        int packedComponent = 0;
        for (const SourceArray& array : arrays) {
            for (int component = 0;
                 component < array.input->GetNumberOfComponents(); ++component) {
                packed->SetComponent(pointId, packedComponent++,
                                     array.input->GetComponent(pointId, component));
            }
        }
    }

    streamcenterplus::io::Nek5000SpectralPointEvaluator evaluator(
        preparedInput, packed);
    vtkSMPTools::For(
        static_cast<vtkIdType>(0), outputPointCount,
        [&](vtkIdType begin, vtkIdType end) {
            std::vector<double> tuple(totalComponents, 0.0);
            int elementHint = -1;
            double point[3] = {};
            for (vtkIdType pointId = begin; pointId < end; ++pointId) {
                output->GetPoint(pointId, point);
                if (!evaluator.Evaluate(point, tuple.data(), &elementHint)) {
                    continue;
                }
                validMask[static_cast<std::size_t>(pointId)] = 1;
                std::size_t packedComponent = 0;
                for (const SourceArray& array : arrays) {
                    const int components = array.input->GetNumberOfComponents();
                    for (int component = 0; component < components; ++component) {
                        array.output->SetComponent(
                            pointId, component, tuple[packedComponent++]);
                    }
                }
            }
        });
    MarkGhostArrays(output, validMask);
    return output;
}

double SquaredDistance(vtkDataSet* dataSet, vtkIdType pointId, const double point[3]) {
    double sourcePoint[3] = {};
    dataSet->GetPoint(pointId, sourcePoint);
    const double dx = sourcePoint[0] - point[0];
    const double dy = sourcePoint[1] - point[1];
    const double dz = sourcePoint[2] - point[2];
    return dx * dx + dy * dy + dz * dz;
}

std::vector<vtkIdType> CollectNeighborIds(vtkDataSet* source,
                                          vtkStaticPointLocator* locator,
                                          vtkUnsignedCharArray* sourceGhost,
                                          const double point[3],
                                          NeighborMode mode,
                                          int neighbors,
                                          double radius) {
    auto ids = vtkSmartPointer<vtkIdList>::New();
    if (mode == NeighborMode::Radius && radius > 0.0) {
        locator->FindPointsWithinRadius(radius, point, ids);
    } else {
        locator->FindClosestNPoints(std::max(1, neighbors), point, ids);
    }

    std::vector<std::pair<double, vtkIdType>> sorted;
    sorted.reserve(static_cast<std::size_t>(ids->GetNumberOfIds()));
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
        const vtkIdType id = ids->GetId(i);
        if (sourceGhost != nullptr &&
            (sourceGhost->GetValue(id) & vtkDataSetAttributes::HIDDENPOINT) != 0) {
            continue;
        }
        sorted.push_back({SquaredDistance(source, id, point), id});
    }
    std::sort(sorted.begin(), sorted.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.first != rhs.first) {
            return lhs.first < rhs.first;
        }
        return lhs.second < rhs.second;
    });
    if (static_cast<int>(sorted.size()) > neighbors) {
        sorted.resize(static_cast<std::size_t>(neighbors));
    }

    std::vector<vtkIdType> result;
    result.reserve(sorted.size());
    for (const auto& item : sorted) {
        result.push_back(item.second);
    }
    return result;
}

double VariogramValue(KrigingVariogram model,
                      double h,
                      double range,
                      double sill,
                      double nugget) {
    if (h <= 1e-14) {
        return 0.0;
    }
    const double effectiveRange = std::max(range, 1e-14);
    const double r = h / effectiveRange;
    switch (model) {
    case KrigingVariogram::Spherical:
        if (r >= 1.0) {
            return nugget + sill;
        }
        return nugget + sill * (1.5 * r - 0.5 * r * r * r);
    case KrigingVariogram::Exponential:
        return nugget + sill * (1.0 - std::exp(-3.0 * r));
    case KrigingVariogram::Gaussian:
        return nugget + sill * (1.0 - std::exp(-3.0 * r * r));
    case KrigingVariogram::Linear:
        return nugget + sill * std::min(1.0, r);
    }
    return nugget + sill;
}

bool SolveLinearSystem(std::vector<double> matrix, std::vector<double> rhs, std::vector<double>* solution) {
    const int n = static_cast<int>(rhs.size());
    if (static_cast<int>(matrix.size()) != n * n || solution == nullptr) {
        return false;
    }
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double pivotValue = std::abs(matrix[static_cast<std::size_t>(col * n + col)]);
        for (int row = col + 1; row < n; ++row) {
            const double value = std::abs(matrix[static_cast<std::size_t>(row * n + col)]);
            if (value > pivotValue) {
                pivot = row;
                pivotValue = value;
            }
        }
        if (pivotValue < 1e-14) {
            return false;
        }
        if (pivot != col) {
            for (int j = col; j < n; ++j) {
                std::swap(matrix[static_cast<std::size_t>(col * n + j)],
                          matrix[static_cast<std::size_t>(pivot * n + j)]);
            }
            std::swap(rhs[static_cast<std::size_t>(col)], rhs[static_cast<std::size_t>(pivot)]);
        }
        const double diag = matrix[static_cast<std::size_t>(col * n + col)];
        for (int row = col + 1; row < n; ++row) {
            const double factor = matrix[static_cast<std::size_t>(row * n + col)] / diag;
            if (factor == 0.0) {
                continue;
            }
            matrix[static_cast<std::size_t>(row * n + col)] = 0.0;
            for (int j = col + 1; j < n; ++j) {
                matrix[static_cast<std::size_t>(row * n + j)] -=
                    factor * matrix[static_cast<std::size_t>(col * n + j)];
            }
            rhs[static_cast<std::size_t>(row)] -= factor * rhs[static_cast<std::size_t>(col)];
        }
    }

    solution->assign(static_cast<std::size_t>(n), 0.0);
    for (int row = n - 1; row >= 0; --row) {
        double sum = rhs[static_cast<std::size_t>(row)];
        for (int col = row + 1; col < n; ++col) {
            sum -= matrix[static_cast<std::size_t>(row * n + col)] *
                   (*solution)[static_cast<std::size_t>(col)];
        }
        const double diag = matrix[static_cast<std::size_t>(row * n + row)];
        if (std::abs(diag) < 1e-14) {
            return false;
        }
        (*solution)[static_cast<std::size_t>(row)] = sum / diag;
    }
    return true;
}

std::vector<double> IdwWeights(vtkDataSet* source,
                               const std::vector<vtkIdType>& neighborIds,
                               const double point[3],
                               double power) {
    std::vector<double> weights(neighborIds.size(), 0.0);
    if (neighborIds.empty()) {
        return weights;
    }
    for (std::size_t i = 0; i < neighborIds.size(); ++i) {
        const double distanceSquared = SquaredDistance(source, neighborIds[i], point);
        if (distanceSquared <= 1e-24) {
            weights.assign(neighborIds.size(), 0.0);
            weights[i] = 1.0;
            return weights;
        }
        weights[i] = 1.0 / std::pow(std::sqrt(distanceSquared), power);
    }
    const double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sum > 0.0 && std::isfinite(sum)) {
        for (double& weight : weights) {
            weight /= sum;
        }
    }
    return weights;
}

std::vector<double> ClosestWeights(const std::vector<vtkIdType>& neighborIds) {
    std::vector<double> weights(neighborIds.size(), 0.0);
    if (!weights.empty()) {
        weights.front() = 1.0;
    }
    return weights;
}

bool KrigingWeights(vtkDataSet* source,
                    const KrigingParameters& params,
                    const std::vector<vtkIdType>& neighborIds,
                    const double point[3],
                    double effectiveRange,
                    std::vector<double>* weights) {
    const int n = static_cast<int>(neighborIds.size());
    if (n <= 0 || weights == nullptr) {
        return false;
    }
    for (int i = 0; i < n; ++i) {
        if (SquaredDistance(source, neighborIds[static_cast<std::size_t>(i)], point) <= 1e-24) {
            weights->assign(static_cast<std::size_t>(n), 0.0);
            (*weights)[static_cast<std::size_t>(i)] = 1.0;
            return true;
        }
    }

    const int systemSize = n + 1;
    std::vector<double> matrix(static_cast<std::size_t>(systemSize * systemSize), 0.0);
    std::vector<double> rhs(static_cast<std::size_t>(systemSize), 0.0);
    for (int i = 0; i < n; ++i) {
        double pi[3] = {};
        source->GetPoint(neighborIds[static_cast<std::size_t>(i)], pi);
        for (int j = 0; j < n; ++j) {
            double pj[3] = {};
            source->GetPoint(neighborIds[static_cast<std::size_t>(j)], pj);
            const double dx = pi[0] - pj[0];
            const double dy = pi[1] - pj[1];
            const double dz = pi[2] - pj[2];
            const double h = std::sqrt(dx * dx + dy * dy + dz * dz);
            matrix[static_cast<std::size_t>(i * systemSize + j)] =
                VariogramValue(params.variogram, h, effectiveRange, params.sill, params.nugget);
        }
        matrix[static_cast<std::size_t>(i * systemSize + i)] += params.regularization;
        matrix[static_cast<std::size_t>(i * systemSize + n)] = 1.0;
        matrix[static_cast<std::size_t>(n * systemSize + i)] = 1.0;

        const double dx = pi[0] - point[0];
        const double dy = pi[1] - point[1];
        const double dz = pi[2] - point[2];
        rhs[static_cast<std::size_t>(i)] =
            VariogramValue(params.variogram,
                           std::sqrt(dx * dx + dy * dy + dz * dz),
                           effectiveRange,
                           params.sill,
                           params.nugget);
    }
    rhs[static_cast<std::size_t>(n)] = 1.0;

    std::vector<double> solution;
    if (!SolveLinearSystem(std::move(matrix), std::move(rhs), &solution)) {
        return false;
    }
    weights->assign(solution.begin(), solution.begin() + n);
    return true;
}

void InterpolateArrays(const std::vector<SourceArray>& arrays,
                       const std::vector<vtkIdType>& neighborIds,
                       const std::vector<double>& weights,
                       vtkIdType outputPointId) {
    for (const SourceArray& array : arrays) {
        const int components = array.input->GetNumberOfComponents();
        for (int component = 0; component < components; ++component) {
            double value = 0.0;
            for (std::size_t i = 0; i < neighborIds.size(); ++i) {
                value += weights[i] * array.input->GetComponent(neighborIds[i], component);
            }
            array.output->SetComponent(outputPointId, component, value);
        }
    }
}

vtkSmartPointer<vtkImageData> ResampleKriging(vtkDataSet* preparedInput,
                                              const Parameters& params,
                                              const std::array<double, 6>& bounds) {
    vtkSmartPointer<vtkImageData> output = CreateTargetImage(params, bounds);
    output->GetFieldData()->PassData(preparedInput->GetFieldData());

    const std::vector<unsigned char> domainMask = ProbeDomainMask(output, preparedInput);
    std::vector<unsigned char> finalMask = domainMask;
    const vtkIdType outputPointCount = output->GetNumberOfPoints();

    std::vector<SourceArray> arrays = PrepareOutputArrays(preparedInput->GetPointData(),
                                                          output->GetPointData(),
                                                          outputPointCount,
                                                          params.kriging.nullValue);

    auto locator = vtkSmartPointer<vtkStaticPointLocator>::New();
    locator->SetDataSet(preparedInput);
    locator->BuildLocator();

    vtkUnsignedCharArray* sourceGhost = preparedInput->GetPointGhostArray();
    const double diagonal = std::sqrt((bounds[1] - bounds[0]) * (bounds[1] - bounds[0]) +
                                      (bounds[3] - bounds[2]) * (bounds[3] - bounds[2]) +
                                      (bounds[5] - bounds[4]) * (bounds[5] - bounds[4]));
    const double effectiveRange = params.kriging.range > 0.0
                                      ? params.kriging.range
                                      : std::max(diagonal / 3.0, 1e-12);

    const auto progressStart = std::chrono::steady_clock::now();
    for (vtkIdType pointId = 0; pointId < outputPointCount; ++pointId) {
        if (ShouldLogProgress(pointId, outputPointCount)) {
            LogProgress("Kriging interpolation", pointId, outputPointCount, progressStart);
        }
        if (domainMask[static_cast<std::size_t>(pointId)] == 0) {
            finalMask[static_cast<std::size_t>(pointId)] = 0;
            continue;
        }

        double point[3] = {};
        output->GetPoint(pointId, point);
        std::vector<vtkIdType> neighbors = CollectNeighborIds(preparedInput,
                                                              locator,
                                                              sourceGhost,
                                                              point,
                                                              params.kriging.neighborMode,
                                                              params.kriging.neighbors,
                                                              params.kriging.radius);
        if (static_cast<int>(neighbors.size()) < params.kriging.minNeighbors) {
            if (params.kriging.fallback == KrigingFallback::Error) {
                throw std::runtime_error("Kriging interpolation has fewer neighbors than min_neighbors.");
            }
            if (neighbors.empty()) {
                finalMask[static_cast<std::size_t>(pointId)] = 0;
                continue;
            }
        }

        std::vector<double> weights;
        const bool solved = static_cast<int>(neighbors.size()) >= params.kriging.minNeighbors &&
                            KrigingWeights(preparedInput,
                                           params.kriging,
                                           neighbors,
                                           point,
                                           effectiveRange,
                                           &weights);
        if (!solved) {
            if (params.kriging.fallback == KrigingFallback::Error) {
                throw std::runtime_error("Kriging linear system failed and fallback=error.");
            }
            weights = params.kriging.fallback == KrigingFallback::Closest
                          ? ClosestWeights(neighbors)
                          : IdwWeights(preparedInput, neighbors, point, params.kriging.fallbackIdwPower);
        }
        if (weights.empty()) {
            finalMask[static_cast<std::size_t>(pointId)] = 0;
            continue;
        }
        InterpolateArrays(arrays, neighbors, weights, pointId);
    }

    MarkGhostArrays(output, finalMask);
    return output;
}

ImageStats ReadStats(vtkImageData* image) {
    if (image == nullptr) {
        return {};
    }
    ImageStats stats;
    stats.pointCount = image->GetNumberOfPoints();
    vtkUnsignedCharArray* pointGhost = image->GetPointGhostArray();
    for (vtkIdType i = 0; i < stats.pointCount; ++i) {
        const bool hidden = pointGhost != nullptr &&
                            (pointGhost->GetValue(i) & vtkDataSetAttributes::HIDDENPOINT) != 0;
        if (!hidden) {
            ++stats.validPointCount;
        }
    }
    vtkUnsignedCharArray* cellGhost = image->GetCellGhostArray();
    if (cellGhost != nullptr) {
        for (vtkIdType i = 0; i < cellGhost->GetNumberOfTuples(); ++i) {
            if ((cellGhost->GetValue(i) & vtkDataSetAttributes::HIDDENCELL) != 0) {
                ++stats.hiddenCellCount;
            }
        }
    }
    return stats;
}

std::string FrameFilename(const std::string& basename, int frameIndex) {
    std::ostringstream out;
    out << basename << "_f" << std::setw(6) << std::setfill('0') << frameIndex << ".vti";
    return out.str();
}

bool ResumeExistingOutputRequested(const Parameters& params) {
    return LowerAscii(params.solverStatus) == "error";
}

double BoundsComparisonTolerance(double expected, double actual) {
    const double scale = std::max(1.0, std::max(std::abs(expected), std::abs(actual)));
    return std::max(1e-8, 1e-7 * scale);
}

bool BoundsNearlyEqual(const std::array<double, 6>& expected,
                       const std::array<double, 6>& actual,
                       std::string* detail) {
    for (int i = 0; i < 6; ++i) {
        const double delta = std::abs(expected[static_cast<std::size_t>(i)] -
                                      actual[static_cast<std::size_t>(i)]);
        const double tolerance = BoundsComparisonTolerance(expected[static_cast<std::size_t>(i)],
                                                           actual[static_cast<std::size_t>(i)]);
        if (delta <= tolerance) {
            continue;
        }
        if (detail != nullptr) {
            std::ostringstream message;
            message << "bounds component " << i << " differs by " << std::setprecision(8) << delta
                    << " (tolerance " << tolerance << ")";
            *detail = message.str();
        }
        return false;
    }
    return true;
}

std::optional<int> OutputFrameIndex(const std::filesystem::path& path, const std::string& basename) {
    const std::string fileName = path.filename().string();
    const std::string prefix = basename + "_f";
    constexpr const char* suffix = ".vti";
    constexpr std::size_t suffixLength = 4;
    if (fileName.size() <= prefix.size() + suffixLength ||
        fileName.compare(0, prefix.size(), prefix) != 0 ||
        fileName.compare(fileName.size() - suffixLength, suffixLength, suffix) != 0) {
        return std::nullopt;
    }
    const std::string digits = fileName.substr(prefix.size(),
                                               fileName.size() - prefix.size() - suffixLength);
    if (digits.empty() || !std::all_of(digits.begin(), digits.end(), [](unsigned char ch) {
            return std::isdigit(ch) != 0;
        })) {
        return std::nullopt;
    }
    try {
        return std::stoi(digits);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

int RemoveResumeOutputsFromFrame(const std::filesystem::path& outputDir,
                                 const std::string& basename,
                                 std::size_t startFrame) {
    if (!std::filesystem::exists(outputDir)) {
        return 0;
    }
    if (!std::filesystem::is_directory(outputDir)) {
        throw std::runtime_error("Output path exists but is not a directory: " + outputDir.string());
    }

    int removed = 0;
    const std::string pvdName = basename + ".pvd";
    const std::string summaryName = basename + "_run_summary.txt";
    for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(outputDir)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string fileName = entry.path().filename().string();
        bool removeFile = fileName == pvdName || fileName == summaryName;
        if (!removeFile) {
            if (std::optional<int> frameIndex = OutputFrameIndex(entry.path(), basename)) {
                removeFile = *frameIndex >= 0 &&
                             static_cast<std::size_t>(*frameIndex) >= startFrame;
            }
        }
        if (!removeFile) {
            continue;
        }
        std::filesystem::remove(entry.path());
        ++removed;
    }
    return removed;
}

std::array<double, 6> ExpectedOutputBoundsForFrame(const Parameters& params,
                                                   const SequenceInfo& sequence,
                                                   const FrameInfo& frame) {
    if (sequence.seriesBoundsPrecomputed) {
        return sequence.bounds;
    }
    return ResolveConfiguredBounds(params, ReadDataSetBounds(frame));
}

struct ValidatedResumeFrame {
    FrameResult result;
    std::array<double, 6> expectedBounds = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

ValidatedResumeFrame ValidateResumeFrame(const Parameters& params,
                                         const SequenceInfo& sequence,
                                         std::size_t frameIndex,
                                         const std::filesystem::path& outputPath) {
    vtkSmartPointer<vtkDataSet> dataSet = ReadDataSet(outputPath);
    vtkImageData* image = vtkImageData::SafeDownCast(dataSet);
    if (image == nullptr) {
        throw std::runtime_error("Existing resume output is not vtkImageData: " + outputPath.string());
    }
    int dims[3] = {0, 0, 0};
    image->GetDimensions(dims);
    for (int axis = 0; axis < 3; ++axis) {
        if (dims[axis] == params.dimensions[static_cast<std::size_t>(axis)]) {
            continue;
        }
        std::ostringstream message;
        message << "Existing resume output dimensions do not match target for frame " << frameIndex
                << ": got " << dims[0] << " x " << dims[1] << " x " << dims[2]
                << ", expected " << params.dimensions[0] << " x " << params.dimensions[1]
                << " x " << params.dimensions[2] << ".";
        throw std::runtime_error(message.str());
    }

    const std::array<double, 6> expectedBounds =
        ExpectedOutputBoundsForFrame(params, sequence, sequence.frames.at(frameIndex));
    const std::array<double, 6> actualBounds = BoundsOf(image);
    std::string boundsDetail;
    if (!BoundsNearlyEqual(expectedBounds, actualBounds, &boundsDetail)) {
        std::ostringstream message;
        message << "Existing resume output bounds do not match target for frame " << frameIndex
                << ": expected " << FormatBounds(expectedBounds)
                << ", got " << FormatBounds(actualBounds)
                << " (" << boundsDetail << ").";
        throw std::runtime_error(message.str());
    }

    const ImageStats stats = ReadStats(image);
    ValidatedResumeFrame validated;
    validated.expectedBounds = expectedBounds;
    validated.result.frameIndex = static_cast<int>(frameIndex);
    validated.result.time = sequence.frames.at(frameIndex).time;
    validated.result.relativeFile = FrameFilename(params.outputBasename, static_cast<int>(frameIndex));
    validated.result.bounds = expectedBounds;
    validated.result.pointCount = static_cast<long long>(stats.pointCount);
    validated.result.validPointCount = static_cast<long long>(stats.validPointCount);
    validated.result.hiddenCellCount = static_cast<long long>(stats.hiddenCellCount);
    return validated;
}

struct ResumePlan {
    std::size_t startFrame = 0;
    vtkIdType completedWork = 0;
    vtkIdType totalWork = 0;
    std::vector<FrameResult> reusedResults;
    std::vector<std::pair<double, std::filesystem::path>> reusedPvdEntries;
};

ResumePlan BuildResumePlan(const Parameters& params,
                           SequenceInfo* sequence,
                           const std::chrono::steady_clock::time_point& progressStart) {
    if (sequence == nullptr) {
        throw std::runtime_error("Internal error: null sequence for resume plan.");
    }

    ResumePlan plan;
    plan.totalWork = static_cast<vtkIdType>(sequence->frames.size());
    std::cout << "[Resume] Existing output resume requested; target folder will be inspected before cleanup.\n";

    std::size_t candidateCount = 0;
    for (; candidateCount < sequence->frames.size(); ++candidateCount) {
        const std::filesystem::path outputPath =
            params.outputDir / FrameFilename(params.outputBasename, static_cast<int>(candidateCount));
        if (!std::filesystem::is_regular_file(outputPath)) {
            break;
        }
    }

    const auto validationStage = LogStageStarted("Restart validation",
                                                 "candidate_frames=" + std::to_string(candidateCount));
    if (candidateCount == 0) {
        LogStageFinished("Restart validation", validationStage, "usable_prefix=0");
        const auto cleanupStage = LogStageStarted("Output cleanup", params.outputDir.string());
        const int removedStaleFiles = RemoveStaleVtiAndPvdFiles(params.outputDir);
        LogStageFinished("Output cleanup",
                         cleanupStage,
                         "removed_vti_pvd_files=" + std::to_string(removedStaleFiles));
        std::cout << "[Resume] No contiguous completed frame outputs were found; starting from frame 0.\n";
        return plan;
    }

    plan.totalWork = static_cast<vtkIdType>(sequence->frames.size() + 1);
    std::vector<ValidatedResumeFrame> validatedFrames;
    validatedFrames.reserve(candidateCount);
    for (std::size_t i = 0; i < candidateCount; ++i) {
        const std::filesystem::path outputPath =
            params.outputDir / FrameFilename(params.outputBasename, static_cast<int>(i));
        std::cout << "[Resume] Validating existing frame " << i
                  << " range: " << outputPath.filename().string() << "\n";
        try {
            ValidatedResumeFrame validated = ValidateResumeFrame(params, *sequence, i, outputPath);
            if (!sequence->seriesBoundsPrecomputed) {
                MergeBounds(validated.expectedBounds, &sequence->bounds, &sequence->boundsInitialized);
            }
            validatedFrames.push_back(validated);
        } catch (const std::exception& exception) {
            LogStageFinished("Restart validation",
                             validationStage,
                             "usable_prefix=" + std::to_string(validatedFrames.size()) +
                                 ", validation_failed=true");
            std::cout << "[Resume] Existing frame " << i << " failed validation: "
                      << exception.what() << "\n";
            const auto cleanupStage = LogStageStarted("Output cleanup", params.outputDir.string());
            const int removedStaleFiles = RemoveStaleVtiAndPvdFiles(params.outputDir);
            LogStageFinished("Output cleanup",
                             cleanupStage,
                             "removed_vti_pvd_files=" + std::to_string(removedStaleFiles));
            std::cout << "[Resume] Existing outputs cannot be reused; starting from frame 0.\n";
            plan.completedWork = 0;
            plan.totalWork = static_cast<vtkIdType>(sequence->frames.size());
            return plan;
        }

        plan.completedWork = static_cast<vtkIdType>(validatedFrames.size());
        std::ostringstream detail;
        detail << "validated snapshot " << validatedFrames.size() << "/" << sequence->frames.size()
               << " snapshot";
        LogCompletedProgress("Restart validation",
                             plan.completedWork,
                             plan.totalWork,
                             progressStart,
                             detail.str());
        std::cout << "[Resume] Existing frame " << i << " passed range validation.\n";
    }
    LogStageFinished("Restart validation",
                     validationStage,
                     "usable_prefix=" + std::to_string(validatedFrames.size()));

    plan.startFrame = validatedFrames.size() - 1;
    plan.reusedResults.reserve(plan.startFrame);
    plan.reusedPvdEntries.reserve(plan.startFrame);
    for (std::size_t i = 0; i < plan.startFrame; ++i) {
        plan.reusedResults.push_back(validatedFrames[i].result);
        plan.reusedPvdEntries.push_back({validatedFrames[i].result.time, validatedFrames[i].result.relativeFile});
    }

    const int removedRestartFiles = RemoveResumeOutputsFromFrame(params.outputDir,
                                                                 params.outputBasename,
                                                                 plan.startFrame);
    std::cout << "[Resume] Dropped validated frame " << plan.startFrame
              << " and newer outputs before continuing; removed_files=" << removedRestartFiles << ".\n";
    std::cout << "[Resume] Reusing " << plan.reusedResults.size()
              << " validated frame output(s); continuing at frame " << plan.startFrame << ".\n";
    return plan;
}

std::string BuildRunSummary(const Parameters& params,
                            const SequenceInfo& sequence,
                            const std::vector<FrameResult>& results,
                            double totalSeconds) {
    std::ostringstream out;
    out << "Streamcenter+ resample CPU solver summary\n";
    out << "Input: " << params.inputPath.string() << "\n";
    out << "Output: " << params.outputDir.string() << "\n";
    out << "PVD: " << (params.outputBasename + ".pvd") << "\n";
    out << "Frames: " << sequence.frames.size() << "\n";
    out << "Dimensions: " << params.dimensions[0] << " "
        << params.dimensions[1] << " " << params.dimensions[2] << "\n";
    out << "Bounds mode: " << ToString(params.boundsMode);
    if (!sequence.seriesBoundsPrecomputed) {
        out << " (streaming frame bounds)";
    }
    out << "\n";
    out << (sequence.seriesBoundsPrecomputed ? "Global bounds: " : "Bounds summary: ")
        << FormatBounds(sequence.bounds) << "\n";
    out << "Interpolation: " << ToString(params.method) << "\n";
    out << "VTK SMP: " << FormatSmpSummary(EnsureVtkSmpConfigured()) << "\n";
    out << "Elapsed seconds: " << std::setprecision(8) << totalSeconds << "\n\n";
    out << "frame,time,file,valid_points,total_points,hidden_cells,"
        << "read_seconds,prepare_seconds,resample_seconds,write_seconds,stats_seconds,total_seconds,bounds\n";
    for (const FrameResult& result : results) {
        out << result.frameIndex << ","
            << std::setprecision(17) << result.time << ","
            << result.relativeFile.generic_string() << ","
            << result.validPointCount << ","
            << result.pointCount << ","
            << result.hiddenCellCount << ","
            << std::setprecision(8) << result.readSeconds << ","
            << result.prepareSeconds << ","
            << result.resampleSeconds << ","
            << result.writeSeconds << ","
            << result.statsSeconds << ","
            << result.totalSeconds << ",\""
            << FormatBounds(result.bounds) << "\"\n";
    }
    return out.str();
}

}  // namespace

vtkSmartPointer<vtkImageData> ResampleFrame(vtkDataSet* preparedInput,
                                            const Parameters& params,
                                            const std::array<double, 6>& bounds) {
    if (preparedInput == nullptr) {
        throw std::runtime_error("Internal error: null prepared input.");
    }
    EnsureVtkSmpConfigured();
    if (streamcenterplus::io::HasNek5000SpectralMetadata(preparedInput)) {
        (void)streamcenterplus::io::GetNek5000SpectralMetadata(preparedInput);
        return ResampleNek5000Spectral(preparedInput, params, bounds);
    }
    switch (params.method) {
    case InterpolationMethod::Linear:
        return ResampleLinear(preparedInput, params, bounds);
    case InterpolationMethod::InverseDistanceWeighting:
        return ResampleIdw(preparedInput, params, bounds);
    case InterpolationMethod::Kriging:
        return ResampleKriging(preparedInput, params, bounds);
    case InterpolationMethod::VtkResampleToImage:
        return ResampleVtkToImage(preparedInput, params, bounds);
    }
    throw std::runtime_error("Unsupported interpolation method.");
}

std::string VtkSmpRuntimeSummary() {
    return FormatSmpSummary(EnsureVtkSmpConfigured());
}

int RunVtkSmpSmokeTest() {
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);

    const VtkSmpConfiguration& config = EnsureVtkSmpConfigured();
    std::mutex mutex;
    std::set<std::size_t> threadIds;
    std::set<unsigned int> processorGroups;
    std::atomic<unsigned long long> checksum{0};

    struct Worker {
        std::mutex* mutex = nullptr;
        std::set<std::size_t>* threadIds = nullptr;
        std::set<unsigned int>* processorGroups = nullptr;
        std::atomic<unsigned long long>* checksum = nullptr;

        void RecordThread() const {
            if (mutex == nullptr || threadIds == nullptr || processorGroups == nullptr) {
                return;
            }
            const std::size_t id = std::hash<std::thread::id>{}(std::this_thread::get_id());
            const unsigned int group = CurrentProcessorGroup();
            std::lock_guard<std::mutex> lock(*mutex);
            threadIds->insert(id);
            processorGroups->insert(group);
        }

        void Initialize() const {
            RecordThread();
        }

        void Reduce() const {
        }

        void operator()(vtkIdType begin, vtkIdType end) const {
            RecordThread();
            unsigned long long local = 1469598103934665603ull;
            for (vtkIdType i = begin; i < end; ++i) {
                unsigned long long value = static_cast<unsigned long long>(i + 1);
                for (int spin = 0; spin < 256; ++spin) {
                    value ^= value << 13;
                    value ^= value >> 7;
                    value ^= value << 17;
                    local ^= value + static_cast<unsigned long long>(spin);
                    local *= 1099511628211ull;
                }
            }
            if (checksum != nullptr) {
                checksum->fetch_add(local, std::memory_order_relaxed);
            }
        }
    };

    const vtkIdType threadCount = static_cast<vtkIdType>(std::max(1, config.requestedThreadCount));
    const vtkIdType iterations = std::max<vtkIdType>(threadCount * 8192, 8192);
    vtkSMPTools::For(static_cast<vtkIdType>(0),
                     iterations,
                     static_cast<vtkIdType>(1),
                     Worker{&mutex, &threadIds, &processorGroups, &checksum});

    const int observedThreads = static_cast<int>(threadIds.size());
    const int observedGroups = static_cast<int>(processorGroups.size());
    std::cout << "VTK SMP smoke test\n";
    std::cout << "  " << FormatSmpSummary(config) << "\n";
    std::cout << "  observed_worker_threads=" << observedThreads << "\n";
    std::cout << "  observed_processor_groups=" << observedGroups << "\n";
    std::cout << "  checksum=" << checksum.load(std::memory_order_relaxed) << "\n";

    const int expectedObservedThreads = std::min(config.requestedThreadCount, config.estimatedThreadCount);
    if (config.requestedThreadCount > 1 && observedThreads < expectedObservedThreads) {
        std::cerr << "VTK SMP smoke test failed: observed " << observedThreads
                  << " worker threads, expected at least " << expectedObservedThreads << ".\n";
        return 1;
    }
    if (config.processorGroupCount > 1 && observedGroups < 2) {
        std::cerr << "VTK SMP smoke test failed: worker threads did not run across processor groups.\n";
        return 1;
    }

    std::cout << "VTK SMP smoke test passed.\n";
    return 0;
}

int RunSolver(const Parameters& params) {
    std::cout.setf(std::ios::unitbuf);
    std::cerr.setf(std::ios::unitbuf);

    std::cout << "Resample solver started.\n\n";
    std::cout << "==============================================\n";
    std::cout << "=== Streamcenter+ VTK Resample CPU Solver ===\n";
    std::cout << "==============================================\n";
    std::cout << "VTK SMP runtime: " << VtkSmpRuntimeSummary() << "\n";
    std::cout << PrettyPrint(params);

    if (!params.enabled) {
        std::cout << "Resample solver is disabled by deck.\n";
        return 0;
    }
    if (params.backend == Backend::Gpu) {
        throw std::runtime_error("GPU backend was requested, but the resample GPU solver is not implemented in this release.");
    }

    const auto totalStart = std::chrono::steady_clock::now();
    const bool resumeExistingOutput = ResumeExistingOutputRequested(params);
    std::filesystem::create_directories(params.outputDir);
    if (!resumeExistingOutput) {
        const auto cleanupStage = LogStageStarted("Output cleanup", params.outputDir.string());
        const int removedStaleFiles = RemoveStaleVtiAndPvdFiles(params.outputDir);
        LogStageFinished("Output cleanup",
                         cleanupStage,
                         "removed_vti_pvd_files=" + std::to_string(removedStaleFiles));
    }

    const auto discoverStage = LogStageStarted("Data discovery", params.inputPath.string());
    SequenceInfo sequence = CollectInputSequence(params);
    LogStageFinished("Data discovery",
                     discoverStage,
                     "frames=" + std::to_string(sequence.frames.size()) +
                         (sequence.boundsInitialized
                              ? ", bounds=" + FormatBounds(sequence.bounds)
                              : ", bounds=streaming_per_frame"));

    std::vector<FrameResult> results;
    results.reserve(sequence.frames.size());
    std::vector<std::pair<double, std::filesystem::path>> pvdEntries;
    pvdEntries.reserve(sequence.frames.size());

    const auto progressStart = std::chrono::steady_clock::now();
    ResumePlan resumePlan;
    if (resumeExistingOutput) {
        resumePlan = BuildResumePlan(params, &sequence, progressStart);
        results = resumePlan.reusedResults;
        pvdEntries = resumePlan.reusedPvdEntries;
    } else {
        resumePlan.totalWork = static_cast<vtkIdType>(sequence.frames.size());
    }
    vtkIdType completedWork = resumePlan.completedWork;
    const vtkIdType totalWork = std::max<vtkIdType>(1, resumePlan.totalWork);

    const auto resampleStage = LogStageStarted("Frame resampling",
                                               "method=" + ToString(params.method) +
                                                   (resumeExistingOutput
                                                        ? ", start_frame=" + std::to_string(resumePlan.startFrame)
                                                        : std::string()));
    if (completedWork == 0) {
        LogCompletedProgress("Frame resampling",
                             0,
                             totalWork,
                             progressStart,
                             "computed snapshot 0/" + std::to_string(sequence.frames.size()) + " snapshot");
    } else {
        LogCompletedProgress("Frame resampling",
                             completedWork,
                             totalWork,
                             progressStart,
                             "resume validation accepted");
    }

    auto computeFrameOnce = [&](std::size_t i) -> FrameResult {
        const auto frameStart = std::chrono::steady_clock::now();
        const FrameInfo& frame = sequence.frames[i];
        const auto readStart = std::chrono::steady_clock::now();
        vtkSmartPointer<vtkDataSet> input = ReadDataSet(frame);
        const double readSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - readStart).count();

        const auto prepareStart = std::chrono::steady_clock::now();
        vtkSmartPointer<vtkDataSet> prepared = PreparePointData(input);
        const std::array<double, 6> inputBounds = BoundsOf(prepared);
        const std::array<double, 6> frameBounds =
            sequence.seriesBoundsPrecomputed ? sequence.bounds : ResolveConfiguredBounds(params, inputBounds);
        if (!sequence.seriesBoundsPrecomputed) {
            MergeBounds(frameBounds, &sequence.bounds, &sequence.boundsInitialized);
        }
        vtkSmartPointer<vtkDataSet> transformed = ApplyInputDomainTransform(prepared, params, inputBounds);
        const double prepareSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - prepareStart).count();

        const auto resampleStart = std::chrono::steady_clock::now();
        vtkSmartPointer<vtkImageData> output = ResampleFrame(transformed, params, frameBounds);
        const double resampleSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - resampleStart).count();

        const std::string fileName = FrameFilename(params.outputBasename, static_cast<int>(i));
        const auto writeStart = std::chrono::steady_clock::now();
        WriteVti(params.outputDir / fileName, output);
        const double writeSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - writeStart).count();

        const auto statsStart = std::chrono::steady_clock::now();
        const ImageStats stats = ReadStats(output);
        const double statsSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - statsStart).count();
        const double frameSeconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - frameStart).count();

        FrameResult result;
        result.frameIndex = static_cast<int>(i);
        result.time = frame.time;
        result.relativeFile = fileName;
        result.bounds = frameBounds;
        result.pointCount = static_cast<long long>(stats.pointCount);
        result.validPointCount = static_cast<long long>(stats.validPointCount);
        result.hiddenCellCount = static_cast<long long>(stats.hiddenCellCount);
        result.readSeconds = readSeconds;
        result.prepareSeconds = prepareSeconds;
        result.resampleSeconds = resampleSeconds;
        result.writeSeconds = writeSeconds;
        result.statsSeconds = statsSeconds;
        result.totalSeconds = frameSeconds;

        std::cout << "[Frame] index=" << i
                  << " read=" << FormatSeconds(readSeconds) << "s"
                  << " prepare=" << FormatSeconds(prepareSeconds) << "s"
                  << " resample=" << FormatSeconds(resampleSeconds) << "s"
                  << " write=" << FormatSeconds(writeSeconds) << "s"
                  << " stats=" << FormatSeconds(statsSeconds) << "s"
                  << " total=" << FormatSeconds(frameSeconds) << "s"
                  << ", output=" << fileName << "\n";
        return result;
    };

    auto computeFrameWithRetry = [&](std::size_t i) -> FrameResult {
        std::string lastError;
        for (int attempt = 0; attempt < 2; ++attempt) {
            try {
                return computeFrameOnce(i);
            } catch (const std::exception& exception) {
                lastError = exception.what();
            } catch (...) {
                lastError = "unknown exception";
            }

            const std::filesystem::path outputPath =
                params.outputDir / FrameFilename(params.outputBasename, static_cast<int>(i));
            std::error_code removeError;
            std::filesystem::remove(outputPath, removeError);
            if (attempt == 0) {
                std::cout << "[Frame] index=" << i << " failed: " << lastError << "\n";
                std::cout << "[Frame] Retrying index=" << i << " once after frame error.\n";
            }
        }
        throw std::runtime_error("Frame " + std::to_string(i) +
                                 " failed after one retry: " + lastError);
    };

    for (std::size_t i = resumePlan.startFrame; i < sequence.frames.size(); ++i) {
        FrameResult result = computeFrameWithRetry(i);
        pvdEntries.push_back({result.time, result.relativeFile});
        results.push_back(result);

        ++completedWork;
        if (ShouldLogProgress(completedWork - 1, totalWork)) {
            std::ostringstream detail;
            detail << "computed snapshot " << (i + 1) << "/" << sequence.frames.size()
                   << " snapshot";
            LogCompletedProgress("Frame resampling",
                                 completedWork,
                                 totalWork,
                                 progressStart,
                                 detail.str());
        }
    }
    LogStageFinished("Frame resampling", resampleStage, "frames=" + std::to_string(results.size()));

    const auto writeStage = LogStageStarted("Result writing", params.outputDir.string());
    const std::filesystem::path pvdPath = params.outputDir / (params.outputBasename + ".pvd");
    WritePvd(pvdPath, pvdEntries);
    const double totalSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - totalStart).count();
    WriteTextFile(params.outputDir / (params.outputBasename + "_run_summary.txt"),
                  BuildRunSummary(params, sequence, results, totalSeconds));
    LogStageFinished("Result writing", writeStage, pvdPath.filename().string());

    std::cout << "\n[Stage] Solver finished. total=" << FormatSeconds(totalSeconds)
              << "s, output=" << pvdPath.string() << "\n";
    return 0;
}

}  // namespace streamcenterplus::resample
