#include "streamcenterplus_resample/Core.h"
#include "Integrity.h"
#include "ParameterValidation.h"
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
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPointInterpolator.h>
#include <vtkPointSet.h>
#include <vtkProbeFilter.h>
#include <vtkResampleToImage.h>
#include <vtkSMP.h>
#include <vtkSMPTools.h>
#include <vtkShepardKernel.h>
#include <vtkSmartPointer.h>
#include <vtkStaticPointLocator.h>
#include <vtkType.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVariant.h>

#include <algorithm>
#include <atomic>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
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
#include <unordered_map>
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
    const std::size_t neighborCount = neighborIds.size();
    if (neighborCount == 0 || weights == nullptr ||
        neighborCount > static_cast<std::size_t>(kMaximumKrigingNeighbors)) {
        return false;
    }
    for (std::size_t i = 0; i < neighborCount; ++i) {
        if (SquaredDistance(source, neighborIds[i], point) <= 1e-24) {
            weights->assign(neighborCount, 0.0);
            (*weights)[i] = 1.0;
            return true;
        }
    }

    const std::size_t systemSize = neighborCount + 1;
    std::vector<double> matrix(systemSize * systemSize, 0.0);
    std::vector<double> rhs(systemSize, 0.0);
    for (std::size_t i = 0; i < neighborCount; ++i) {
        double pi[3] = {};
        source->GetPoint(neighborIds[i], pi);
        for (std::size_t j = 0; j < neighborCount; ++j) {
            double pj[3] = {};
            source->GetPoint(neighborIds[j], pj);
            const double dx = pi[0] - pj[0];
            const double dy = pi[1] - pj[1];
            const double dz = pi[2] - pj[2];
            const double h = std::sqrt(dx * dx + dy * dy + dz * dz);
            matrix[i * systemSize + j] =
                VariogramValue(params.variogram, h, effectiveRange, params.sill, params.nugget);
        }
        matrix[i * systemSize + i] += params.regularization;
        matrix[i * systemSize + neighborCount] = 1.0;
        matrix[neighborCount * systemSize + i] = 1.0;

        const double dx = pi[0] - point[0];
        const double dy = pi[1] - point[1];
        const double dz = pi[2] - point[2];
        rhs[i] =
            VariogramValue(params.variogram,
                           std::sqrt(dx * dx + dy * dy + dz * dz),
                           effectiveRange,
                           params.sill,
                           params.nugget);
    }
    rhs[neighborCount] = 1.0;

    std::vector<double> solution;
    if (!SolveLinearSystem(std::move(matrix), std::move(rhs), &solution)) {
        return false;
    }
    weights->assign(solution.begin(), solution.begin() + static_cast<std::ptrdiff_t>(neighborCount));
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

std::optional<int> OutputFrameIndex(const std::filesystem::path& path,
                                    const std::string& basename);

constexpr const char* kResumeManifestFormat = "streamcenterplus_resample_cache_v2";
constexpr const char* kResumeImplementationVersion = "SCP-ResamplerIO-0.3.0-resume-v2";
constexpr std::uintmax_t kMaximumResumeManifestBytes = 16u * 1024u * 1024u;

struct FrameInputIdentity {
    std::string source;
    double time = 0.0;
    std::string contentChecksum;
    std::string schemaChecksum;
    std::string gridChecksum;
};

struct RunIdentity {
    std::string parameterChecksum;
    std::string signature;
    std::vector<FrameInputIdentity> frames;
};

struct ResumeManifest {
    std::string outputBasename;
    std::string parameterChecksum;
    std::string runSignature;
    std::string state = "in_progress";
    std::vector<FrameInputIdentity> frames;
    std::vector<std::string> outputFiles;
    std::vector<std::string> outputChecksums;
    std::string pvdChecksum;
    std::string summaryChecksum;
};

struct OutputLayout {
    std::filesystem::path root;
    std::filesystem::path manifest;
    std::filesystem::path pvd;
    std::filesystem::path summary;
    std::vector<std::filesystem::path> frames;
};

void HashToken(integrity::Sha256* digest, const std::string& value) {
    if (digest == nullptr) {
        throw std::runtime_error("Internal error: null SHA-256 digest.");
    }
    const std::string length = std::to_string(value.size());
    digest->Update(length);
    digest->Update(":");
    digest->Update(value);
    digest->Update(";");
}

void HashInteger(integrity::Sha256* digest, long long value) {
    HashToken(digest, std::to_string(value));
}

void HashDouble(integrity::Sha256* digest, double value) {
    static_assert(sizeof(double) == sizeof(std::uint64_t),
                  "Dataset identity requires an IEEE-754-sized double.");
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    std::ostringstream text;
    text << std::hex << std::setw(16) << std::setfill('0') << bits;
    HashToken(digest, text.str());
}

void HashArraySchema(vtkFieldData* arrays,
                     const char* association,
                     integrity::Sha256* schemaDigest) {
    HashToken(schemaDigest, association != nullptr ? association : "unknown");
    const int arrayCount = arrays != nullptr ? arrays->GetNumberOfArrays() : 0;
    HashInteger(schemaDigest, arrayCount);
    for (int index = 0; index < arrayCount; ++index) {
        vtkAbstractArray* array = arrays->GetAbstractArray(index);
        if (array == nullptr) {
            throw std::runtime_error("Input dataset contains a null data-array entry.");
        }
        HashToken(schemaDigest, array->GetName() != nullptr ? array->GetName() : "");
        HashInteger(schemaDigest, array->GetDataType());
        HashInteger(schemaDigest, array->GetNumberOfComponents());
        HashInteger(schemaDigest, array->GetNumberOfTuples());
        for (int component = 0; component < array->GetNumberOfComponents(); ++component) {
            const char* componentName = array->GetComponentName(component);
            HashToken(schemaDigest, componentName != nullptr ? componentName : "");
        }
    }
    vtkDataSetAttributes* attributes = vtkDataSetAttributes::SafeDownCast(arrays);
    HashInteger(schemaDigest, attributes != nullptr ? 1 : 0);
    if (attributes != nullptr) {
        std::array<int, vtkDataSetAttributes::NUM_ATTRIBUTES> attributeIndices{};
        attributes->GetAttributeIndices(attributeIndices.data());
        for (int index : attributeIndices) {
            HashInteger(schemaDigest, index);
        }
    }
}

bool HashNumericArrayStorage(vtkDataArray* array,
                             integrity::Sha256* contentDigest) {
    if (array == nullptr || array->GetDataType() == VTK_BIT) {
        return false;
    }
    const vtkIdType valueCount = array->GetNumberOfValues();
    const int bytesPerValue = array->GetDataTypeSize();
    if (valueCount < 0 || bytesPerValue <= 0) {
        throw std::runtime_error("Input numeric array has an invalid storage extent.");
    }
    const std::uint64_t count = static_cast<std::uint64_t>(valueCount);
    const std::uint64_t width = static_cast<std::uint64_t>(bytesPerValue);
    if (count != 0 && width > std::numeric_limits<std::uint64_t>::max() / count) {
        throw std::runtime_error("Input numeric array byte count overflows uint64.");
    }
    const std::uint64_t byteCount = count * width;
    if (byteCount > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("Input numeric array byte count exceeds size_t capacity.");
    }
    HashToken(contentDigest, std::to_string(byteCount));
    if (byteCount == 0) {
        return true;
    }
    void* storage = array->GetVoidPointer(0);
    if (storage == nullptr) {
        return false;
    }
    contentDigest->Update(storage, static_cast<std::size_t>(byteCount));
    return true;
}

void HashArrayContents(vtkFieldData* arrays,
                       const char* association,
                       integrity::Sha256* contentDigest) {
    HashToken(contentDigest, association != nullptr ? association : "unknown");
    const int arrayCount = arrays != nullptr ? arrays->GetNumberOfArrays() : 0;
    HashInteger(contentDigest, arrayCount);
    for (int index = 0; index < arrayCount; ++index) {
        vtkAbstractArray* array = arrays->GetAbstractArray(index);
        if (array == nullptr) {
            throw std::runtime_error("Input dataset contains a null data-array entry.");
        }
        HashToken(contentDigest, array->GetName() != nullptr ? array->GetName() : "");
        vtkDataArray* numeric = vtkDataArray::SafeDownCast(array);
        if (numeric != nullptr) {
            if (HashNumericArrayStorage(numeric, contentDigest)) {
                continue;
            }
            HashToken(contentDigest, "component-fallback");
            const vtkIdType tuples = numeric->GetNumberOfTuples();
            const int components = numeric->GetNumberOfComponents();
            for (vtkIdType tuple = 0; tuple < tuples; ++tuple) {
                for (int component = 0; component < components; ++component) {
                    HashDouble(contentDigest, numeric->GetComponent(tuple, component));
                }
            }
            continue;
        }
        const vtkIdType values = array->GetNumberOfValues();
        for (vtkIdType value = 0; value < values; ++value) {
            HashToken(contentDigest, array->GetVariantValue(value).ToString());
        }
    }
}

FrameInputIdentity ComputeFrameInputIdentityFromDataSet(const Parameters& params,
                                                        const SequenceInfo& sequence,
                                                        std::size_t frameIndex,
                                                        vtkDataSet* input,
                                                        bool executablePreflight) {
    const FrameInfo& frame = sequence.frames.at(frameIndex);
    if (input == nullptr) {
        throw std::runtime_error("Cannot fingerprint a null input dataset.");
    }

    integrity::Sha256 gridDigest;
    HashToken(&gridDigest, input->GetClassName() != nullptr ? input->GetClassName() : "");
    HashInteger(&gridDigest, input->GetNumberOfPoints());
    HashInteger(&gridDigest, input->GetNumberOfCells());
    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(input);
    if (pointSet != nullptr && pointSet->GetPoints() != nullptr &&
        pointSet->GetPoints()->GetData() != nullptr) {
        vtkDataArray* coordinates = pointSet->GetPoints()->GetData();
        HashToken(&gridDigest, "explicit-point-storage");
        HashInteger(&gridDigest, coordinates->GetDataType());
        HashInteger(&gridDigest, coordinates->GetNumberOfComponents());
        HashInteger(&gridDigest, coordinates->GetNumberOfTuples());
        if (!HashNumericArrayStorage(coordinates, &gridDigest)) {
            throw std::runtime_error("Cannot access explicit point-coordinate storage for identity.");
        }
    } else {
        HashToken(&gridDigest, "evaluated-point-coordinates");
        for (vtkIdType pointId = 0; pointId < input->GetNumberOfPoints(); ++pointId) {
            double point[3] = {};
            input->GetPoint(pointId, point);
            HashDouble(&gridDigest, point[0]);
            HashDouble(&gridDigest, point[1]);
            HashDouble(&gridDigest, point[2]);
        }
    }
    vtkNew<vtkIdList> pointIds;
    for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
        HashInteger(&gridDigest, input->GetCellType(cellId));
        input->GetCellPoints(cellId, pointIds);
        HashInteger(&gridDigest, pointIds->GetNumberOfIds());
        for (vtkIdType local = 0; local < pointIds->GetNumberOfIds(); ++local) {
            HashInteger(&gridDigest, pointIds->GetId(local));
        }
    }

    integrity::Sha256 schemaDigest;
    HashArraySchema(input->GetPointData(), "point", &schemaDigest);
    HashArraySchema(input->GetCellData(), "cell", &schemaDigest);
    HashArraySchema(input->GetFieldData(), "field", &schemaDigest);

    integrity::Sha256 arrayDigest;
    HashArrayContents(input->GetPointData(), "point", &arrayDigest);
    HashArrayContents(input->GetCellData(), "cell", &arrayDigest);
    HashArrayContents(input->GetFieldData(), "field", &arrayDigest);

    FrameInputIdentity identity;
    std::error_code canonicalError;
    std::filesystem::path source = std::filesystem::weakly_canonical(
        std::filesystem::absolute(frame.path), canonicalError);
    if (canonicalError) {
        source = std::filesystem::absolute(frame.path).lexically_normal();
    }
    std::ostringstream sourceText;
    sourceText << source.generic_string() << "|kind=" << ToString(frame.kind)
               << "|source_index=" << frame.sourceIndex
               << "|coordinate_source_index=" << frame.coordinateSourceIndex
               << "|background_index=" << frame.multiZoneBackgroundIndex;
    if (frame.readerTime.has_value()) {
        sourceText << "|reader_time=" << std::setprecision(17) << *frame.readerTime;
    }
    if (!frame.coordinatePath.empty()) {
        sourceText << "|coordinate_path=" <<
            std::filesystem::absolute(frame.coordinatePath).lexically_normal().generic_string();
    }
    identity.source = sourceText.str();
    identity.time = frame.time;
    identity.gridChecksum = gridDigest.FinalHex();
    identity.schemaChecksum = schemaDigest.FinalHex();
    identity.contentChecksum = integrity::Sha256String(
        identity.gridChecksum + "\n" + identity.schemaChecksum + "\n" + arrayDigest.FinalHex());

    if (executablePreflight) {
        // Complete the non-destructive executable-input preflight before any old
        // run-owned result is removed. This catches unreadable arrays, invalid
        // crop bounds, and unsupported coordinate transforms up front.
        vtkSmartPointer<vtkDataSet> prepared = PreparePointData(input);
        const std::array<double, 6> inputBounds = BoundsOf(prepared);
        (void)(sequence.seriesBoundsPrecomputed
                   ? sequence.bounds
                   : ResolveConfiguredBounds(params, inputBounds));
        vtkSmartPointer<vtkDataSet> transformed =
            ApplyInputDomainTransform(prepared, params, inputBounds);
        if (transformed == nullptr || transformed->GetNumberOfPoints() <= 0) {
            throw std::runtime_error("Input preflight produced no usable points for frame " +
                                     std::to_string(frameIndex) + ".");
        }
    }
    return identity;
}

FrameInputIdentity ComputeFrameInputIdentity(const Parameters& params,
                                             const SequenceInfo& sequence,
                                             std::size_t frameIndex) {
    vtkSmartPointer<vtkDataSet> input = ReadDataSet(sequence.frames.at(frameIndex));
    return ComputeFrameInputIdentityFromDataSet(params, sequence, frameIndex, input, true);
}

std::string ParameterIdentityText(const Parameters& params) {
    const auto precise = [](double value) {
        std::ostringstream text;
        text << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
        return text.str();
    };
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "implementation=" << kResumeImplementationVersion << "\n";
    out << "enabled=" << params.enabled << "\n";
    out << "backend=" << ToString(params.backend) << "\n";
    out << "method=" << ToString(params.method) << "\n";
    out << "dimensions=" << params.dimensions[0] << ',' << params.dimensions[1] << ','
        << params.dimensions[2] << "\n";
    out << "bounds_mode=" << ToString(params.boundsMode) << "\n";
    out << "output_basename=" << params.outputBasename << "\n";
    out << "physical_dt=" << params.physicalDt << "\n";
    out << "frame_start=" << (params.frameRange.start.has_value()
                                  ? std::to_string(*params.frameRange.start) : "auto") << "\n";
    out << "frame_end=" << (params.frameRange.end.has_value()
                                ? std::to_string(*params.frameRange.end) : "auto") << "\n";
    for (int axis = 0; axis < 3; ++axis) {
        const AxisRange& range = params.cropRanges[static_cast<std::size_t>(axis)];
        out << "crop" << axis << "_min="
            << (range.min.has_value() ? precise(*range.min) : "auto") << "\n";
        out << "crop" << axis << "_max="
            << (range.max.has_value() ? precise(*range.max) : "auto") << "\n";
        out << "input_scale" << axis << '=' << params.inputScale[static_cast<std::size_t>(axis)] << "\n";
    }
    out << "linear.compute_tolerance=" << params.linear.computeTolerance << "\n";
    out << "linear.tolerance=" << params.linear.tolerance << "\n";
    out << "linear.snap=" << params.linear.snapToCellWithClosestPoint << "\n";
    out << "idw.neighbor_mode=" << ToString(params.idw.neighborMode) << "\n";
    out << "idw.neighbors=" << params.idw.neighbors << "\n";
    out << "idw.radius=" << params.idw.radius << "\n";
    out << "idw.power=" << params.idw.power << "\n";
    out << "idw.null_strategy=" << ToString(params.idw.nullStrategy) << "\n";
    out << "idw.null_value=" << params.idw.nullValue << "\n";
    out << "kriging.variogram=" << ToString(params.kriging.variogram) << "\n";
    out << "kriging.neighbor_mode=" << ToString(params.kriging.neighborMode) << "\n";
    out << "kriging.neighbors=" << params.kriging.neighbors << "\n";
    out << "kriging.min_neighbors=" << params.kriging.minNeighbors << "\n";
    out << "kriging.radius=" << params.kriging.radius << "\n";
    out << "kriging.range=" << params.kriging.range << "\n";
    out << "kriging.sill=" << params.kriging.sill << "\n";
    out << "kriging.nugget=" << params.kriging.nugget << "\n";
    out << "kriging.regularization=" << params.kriging.regularization << "\n";
    out << "kriging.null_value=" << params.kriging.nullValue << "\n";
    out << "kriging.fallback=" << ToString(params.kriging.fallback) << "\n";
    out << "kriging.fallback_idw_power=" << params.kriging.fallbackIdwPower << "\n";
    return out.str();
}

RunIdentity BuildRunIdentity(const Parameters& params, const SequenceInfo& sequence) {
    RunIdentity identity;
    identity.parameterChecksum = integrity::Sha256String(ParameterIdentityText(params));
    identity.frames.reserve(sequence.frames.size());
    integrity::Sha256 signature;
    HashToken(&signature, kResumeImplementationVersion);
    HashToken(&signature, identity.parameterChecksum);
    for (std::size_t index = 0; index < sequence.frames.size(); ++index) {
        FrameInputIdentity frame = ComputeFrameInputIdentity(params, sequence, index);
        HashInteger(&signature, static_cast<long long>(index));
        HashToken(&signature, frame.source);
        HashDouble(&signature, frame.time);
        HashToken(&signature, frame.contentChecksum);
        HashToken(&signature, frame.schemaChecksum);
        HashToken(&signature, frame.gridChecksum);
        identity.frames.push_back(std::move(frame));
    }
    identity.signature = signature.FinalHex();
    return identity;
}

bool SameFrameInputIdentity(const FrameInputIdentity& lhs,
                            const FrameInputIdentity& rhs) {
    return lhs.source == rhs.source && lhs.time == rhs.time &&
           lhs.contentChecksum == rhs.contentChecksum &&
           lhs.schemaChecksum == rhs.schemaChecksum &&
           lhs.gridChecksum == rhs.gridChecksum;
}

bool SameFrameInputIdentities(const std::vector<FrameInputIdentity>& lhs,
                              const std::vector<FrameInputIdentity>& rhs) {
    if (lhs.size() != rhs.size()) {
        return false;
    }
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        if (!SameFrameInputIdentity(lhs[index], rhs[index])) {
            return false;
        }
    }
    return true;
}

std::filesystem::path CanonicalContainedPath(const std::filesystem::path& root,
                                             const std::string& filename) {
    const std::filesystem::path name(filename);
    if (filename.empty() || name.is_absolute() || name.has_root_name() ||
        name.has_root_directory() || name.has_parent_path() ||
        filename.find('/') != std::string::npos || filename.find('\\') != std::string::npos ||
        filename.find("..") != std::string::npos) {
        throw std::runtime_error("Refusing output filename outside the validated run namespace: " +
                                 filename);
    }
    const std::filesystem::path directCandidate = (root / name).lexically_normal();
    if (directCandidate.parent_path() != root) {
        throw std::runtime_error("Output path is not a direct child of the selected directory: " +
                                 directCandidate.string());
    }
    std::error_code error;
    const std::filesystem::file_status linkStatus =
        std::filesystem::symlink_status(directCandidate, error);
    if (error && linkStatus.type() != std::filesystem::file_type::not_found) {
        throw std::runtime_error("Cannot inspect output path '" + directCandidate.string() +
                                 "': " + error.message());
    }
    if (std::filesystem::is_symlink(linkStatus)) {
        throw std::runtime_error("Output namespace entry must not be a symbolic link/reparse point: " +
                                 directCandidate.string());
    }
    error.clear();
    const std::filesystem::path candidate =
        std::filesystem::weakly_canonical(directCandidate, error);
    if (error) {
        throw std::runtime_error("Cannot canonicalize output path '" + directCandidate.string() +
                                 "': " + error.message());
    }
    const std::filesystem::path relative = candidate.lexically_relative(root);
    if (relative.empty() || relative.is_absolute() || relative.has_root_name() ||
        *relative.begin() == ".." || candidate.parent_path() != root) {
        throw std::runtime_error("Canonical output path escapes the selected output directory: " +
                                 candidate.string());
    }
    if (candidate != directCandidate) {
        throw std::runtime_error(
            "Requested output path aliases an existing filesystem entry with a different "
            "identity: " + candidate.string());
    }
    return directCandidate;
}

OutputLayout BuildOutputLayout(const Parameters& params, std::size_t frameCount) {
    std::error_code error;
    const std::filesystem::path absoluteRoot = std::filesystem::absolute(params.outputDir);
    if (std::filesystem::exists(absoluteRoot, error) && !error &&
        !std::filesystem::is_directory(absoluteRoot)) {
        throw std::runtime_error("Output path exists but is not a directory: " + absoluteRoot.string());
    }
    error.clear();
    const std::filesystem::path root = std::filesystem::weakly_canonical(absoluteRoot, error);
    if (error) {
        throw std::runtime_error("Cannot canonicalize output directory '" + absoluteRoot.string() +
                                 "': " + error.message());
    }
    OutputLayout layout;
    layout.root = root;
    layout.manifest = CanonicalContainedPath(root, params.outputBasename + ".scp-resample-manifest");
    layout.pvd = CanonicalContainedPath(root, params.outputBasename + ".pvd");
    layout.summary = CanonicalContainedPath(root, params.outputBasename + "_run_summary.txt");
    layout.frames.reserve(frameCount);
    for (std::size_t index = 0; index < frameCount; ++index) {
        if (index > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Frame count exceeds output filename index capacity.");
        }
        layout.frames.push_back(CanonicalContainedPath(
            root, FrameFilename(params.outputBasename, static_cast<int>(index))));
    }
    return layout;
}

bool IsChecksum(const std::string& value, bool allowEmpty) {
    if (allowEmpty && value.empty()) {
        return true;
    }
    return value.size() == 64 &&
           std::all_of(value.begin(), value.end(), [](unsigned char ch) {
               return std::isdigit(ch) != 0 || (ch >= 'a' && ch <= 'f');
           });
}

std::string SerializeResumeManifestBody(const ResumeManifest& manifest) {
    if (manifest.frames.size() != manifest.outputFiles.size() ||
        manifest.frames.size() != manifest.outputChecksums.size()) {
        throw std::runtime_error("Internal error: inconsistent resume manifest frame counts.");
    }
    std::ostringstream out;
    out << "format=" << kResumeManifestFormat << "\n";
    out << "implementation=" << integrity::HexEncode(kResumeImplementationVersion) << "\n";
    out << "namespace=" << integrity::HexEncode(manifest.outputBasename) << "\n";
    out << "state=" << manifest.state << "\n";
    out << "parameter_sha256=" << manifest.parameterChecksum << "\n";
    out << "run_sha256=" << manifest.runSignature << "\n";
    out << "frame_count=" << manifest.frames.size() << "\n";
    for (std::size_t index = 0; index < manifest.frames.size(); ++index) {
        const FrameInputIdentity& frame = manifest.frames[index];
        out << "frame." << index << ".source=" << integrity::HexEncode(frame.source) << "\n";
        out << "frame." << index << ".time="
            << std::setprecision(std::numeric_limits<double>::max_digits10) << frame.time << "\n";
        out << "frame." << index << ".content_sha256=" << frame.contentChecksum << "\n";
        out << "frame." << index << ".schema_sha256=" << frame.schemaChecksum << "\n";
        out << "frame." << index << ".grid_sha256=" << frame.gridChecksum << "\n";
        out << "output." << index << ".file=" << integrity::HexEncode(manifest.outputFiles[index]) << "\n";
        out << "output." << index << ".sha256="
            << (manifest.outputChecksums[index].empty() ? "-" : manifest.outputChecksums[index]) << "\n";
    }
    out << "pvd_file=" << integrity::HexEncode(manifest.outputBasename + ".pvd") << "\n";
    out << "pvd_sha256=" << (manifest.pvdChecksum.empty() ? "-" : manifest.pvdChecksum) << "\n";
    out << "summary_file=" << integrity::HexEncode(manifest.outputBasename + "_run_summary.txt") << "\n";
    out << "summary_sha256="
        << (manifest.summaryChecksum.empty() ? "-" : manifest.summaryChecksum) << "\n";
    return out.str();
}

std::string SerializeResumeManifest(const ResumeManifest& manifest) {
    const std::string body = SerializeResumeManifestBody(manifest);
    return body + "manifest_sha256=" + integrity::Sha256String(body) + "\n";
}

std::size_t ParseManifestSize(const std::string& value, const char* description) {
    try {
        std::size_t consumed = 0;
        const unsigned long long parsed = std::stoull(value, &consumed, 10);
        if (consumed != value.size() ||
            parsed > static_cast<unsigned long long>(std::numeric_limits<std::size_t>::max()) ||
            parsed > 1000000ull) {
            throw std::runtime_error("out of range");
        }
        return static_cast<std::size_t>(parsed);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid resume manifest ") + description + ".");
    }
}

double ParseManifestTime(const std::string& value) {
    try {
        std::size_t consumed = 0;
        const double parsed = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(parsed)) {
            throw std::runtime_error("invalid");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid resume manifest frame time.");
    }
}

ResumeManifest ParseResumeManifest(const std::string& text,
                                   const std::string& expectedBasename) {
    const std::string marker = "manifest_sha256=";
    const std::size_t markerPosition = text.rfind(marker);
    if (markerPosition == std::string::npos ||
        (markerPosition != 0 && text[markerPosition - 1] != '\n')) {
        throw std::runtime_error("Resume manifest is missing its integrity checksum.");
    }
    const std::string body = text.substr(0, markerPosition);
    std::string recordedChecksum = text.substr(markerPosition + marker.size());
    while (!recordedChecksum.empty() &&
           (recordedChecksum.back() == '\n' || recordedChecksum.back() == '\r')) {
        recordedChecksum.pop_back();
    }
    if (!IsChecksum(recordedChecksum, false) ||
        recordedChecksum != integrity::Sha256String(body)) {
        throw std::runtime_error("Resume manifest integrity checksum does not match its contents.");
    }

    std::unordered_map<std::string, std::string> values;
    std::istringstream lines(body);
    std::string line;
    while (std::getline(lines, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const std::size_t equal = line.find('=');
        if (equal == std::string::npos || equal == 0 ||
            !values.emplace(line.substr(0, equal), line.substr(equal + 1)).second) {
            throw std::runtime_error("Resume manifest contains a malformed or duplicate entry.");
        }
    }
    auto require = [&](const std::string& key) -> const std::string& {
        const auto found = values.find(key);
        if (found == values.end()) {
            throw std::runtime_error("Resume manifest is missing required key '" + key + "'.");
        }
        return found->second;
    };
    if (require("format") != kResumeManifestFormat ||
        integrity::HexDecode(require("implementation")) != kResumeImplementationVersion) {
        throw std::runtime_error("Resume manifest format/implementation version is unsupported.");
    }

    ResumeManifest manifest;
    manifest.outputBasename = integrity::HexDecode(require("namespace"));
    if (manifest.outputBasename != expectedBasename) {
        throw std::runtime_error("Resume manifest namespace does not match output_basename.");
    }
    manifest.state = require("state");
    if (manifest.state != "in_progress" && manifest.state != "complete") {
        throw std::runtime_error("Resume manifest has an invalid state.");
    }
    manifest.parameterChecksum = require("parameter_sha256");
    manifest.runSignature = require("run_sha256");
    if (!IsChecksum(manifest.parameterChecksum, false) ||
        !IsChecksum(manifest.runSignature, false)) {
        throw std::runtime_error("Resume manifest has an invalid parameter/run checksum.");
    }
    const std::size_t frameCount = ParseManifestSize(require("frame_count"), "frame_count");
    manifest.frames.resize(frameCount);
    manifest.outputFiles.resize(frameCount);
    manifest.outputChecksums.resize(frameCount);
    for (std::size_t index = 0; index < frameCount; ++index) {
        const std::string prefix = "frame." + std::to_string(index) + ".";
        FrameInputIdentity& frame = manifest.frames[index];
        frame.source = integrity::HexDecode(require(prefix + "source"));
        frame.time = ParseManifestTime(require(prefix + "time"));
        frame.contentChecksum = require(prefix + "content_sha256");
        frame.schemaChecksum = require(prefix + "schema_sha256");
        frame.gridChecksum = require(prefix + "grid_sha256");
        if (!IsChecksum(frame.contentChecksum, false) ||
            !IsChecksum(frame.schemaChecksum, false) ||
            !IsChecksum(frame.gridChecksum, false)) {
            throw std::runtime_error("Resume manifest has an invalid input identity checksum.");
        }
        const std::string outputPrefix = "output." + std::to_string(index) + ".";
        manifest.outputFiles[index] = integrity::HexDecode(require(outputPrefix + "file"));
        if (index > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
            manifest.outputFiles[index] !=
                FrameFilename(expectedBasename, static_cast<int>(index))) {
            throw std::runtime_error("Resume manifest contains an output outside its frame namespace.");
        }
        const std::string checksum = require(outputPrefix + "sha256");
        manifest.outputChecksums[index] = checksum == "-" ? std::string() : checksum;
        if (!IsChecksum(manifest.outputChecksums[index], true)) {
            throw std::runtime_error("Resume manifest has an invalid output checksum.");
        }
    }
    if (integrity::HexDecode(require("pvd_file")) != expectedBasename + ".pvd" ||
        integrity::HexDecode(require("summary_file")) !=
            expectedBasename + "_run_summary.txt") {
        throw std::runtime_error("Resume manifest metadata outputs escape their namespace.");
    }
    manifest.pvdChecksum = require("pvd_sha256") == "-" ? std::string() : require("pvd_sha256");
    manifest.summaryChecksum = require("summary_sha256") == "-"
                                   ? std::string() : require("summary_sha256");
    if (!IsChecksum(manifest.pvdChecksum, true) || !IsChecksum(manifest.summaryChecksum, true)) {
        throw std::runtime_error("Resume manifest has an invalid metadata-output checksum.");
    }
    return manifest;
}

std::optional<ResumeManifest> ReadResumeManifest(const OutputLayout& layout,
                                                 const std::string& basename) {
    if (!std::filesystem::exists(layout.manifest)) {
        return std::nullopt;
    }
    if (!std::filesystem::is_regular_file(layout.manifest) ||
        std::filesystem::file_size(layout.manifest) > kMaximumResumeManifestBytes) {
        throw std::runtime_error("Resume manifest is not a bounded regular file: " +
                                 layout.manifest.string());
    }
    std::ifstream input(layout.manifest, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open resume manifest: " + layout.manifest.string());
    }
    std::ostringstream text;
    text << input.rdbuf();
    if (!input.eof() && input.fail()) {
        throw std::runtime_error("Cannot read resume manifest: " + layout.manifest.string());
    }
    return ParseResumeManifest(text.str(), basename);
}

void WriteResumeManifest(const OutputLayout& layout, const ResumeManifest& manifest) {
    WriteTextFile(layout.manifest, SerializeResumeManifest(manifest));
}

ResumeManifest MakeResumeManifest(const Parameters& params,
                                  const RunIdentity& identity) {
    ResumeManifest manifest;
    manifest.outputBasename = params.outputBasename;
    manifest.parameterChecksum = identity.parameterChecksum;
    manifest.runSignature = identity.signature;
    manifest.frames = identity.frames;
    manifest.outputFiles.reserve(identity.frames.size());
    manifest.outputChecksums.resize(identity.frames.size());
    for (std::size_t index = 0; index < identity.frames.size(); ++index) {
        manifest.outputFiles.push_back(
            FrameFilename(params.outputBasename, static_cast<int>(index)));
    }
    return manifest;
}

int RemoveManifestOwnedOutputs(const ResumeManifest& manifest,
                               const std::filesystem::path& root,
                               std::size_t startFrame,
                               bool removeMetadata) {
    int removed = 0;
    for (std::size_t index = startFrame; index < manifest.outputFiles.size(); ++index) {
        const std::filesystem::path path =
            CanonicalContainedPath(root, manifest.outputFiles[index]);
        std::error_code error;
        const std::filesystem::file_status status = std::filesystem::symlink_status(path, error);
        if (error && status.type() != std::filesystem::file_type::not_found) {
            throw std::runtime_error("Cannot inspect manifest-owned output '" + path.string() +
                                     "': " + error.message());
        }
        if (status.type() == std::filesystem::file_type::not_found) {
            continue;
        }
        if (!std::filesystem::is_regular_file(status)) {
            throw std::runtime_error("Refusing to remove a non-regular manifest namespace entry: " +
                                     path.string());
        }
        const bool didRemove = std::filesystem::remove(path, error);
        if (error) {
            throw std::runtime_error("Cannot remove manifest-owned output '" + path.string() +
                                     "': " + error.message());
        }
        removed += didRemove ? 1 : 0;
    }
    if (removeMetadata) {
        for (const std::string& filename : {
                 manifest.outputBasename + ".pvd",
                 manifest.outputBasename + "_run_summary.txt"}) {
            const std::filesystem::path path = CanonicalContainedPath(root, filename);
            std::error_code error;
            const std::filesystem::file_status status = std::filesystem::symlink_status(path, error);
            if (error && status.type() != std::filesystem::file_type::not_found) {
                throw std::runtime_error("Cannot inspect manifest-owned metadata output '" +
                                         path.string() + "': " + error.message());
            }
            if (status.type() == std::filesystem::file_type::not_found) {
                continue;
            }
            if (!std::filesystem::is_regular_file(status)) {
                throw std::runtime_error(
                    "Refusing to remove a non-regular manifest metadata entry: " + path.string());
            }
            const bool didRemove = std::filesystem::remove(path, error);
            if (error) {
                throw std::runtime_error("Cannot remove manifest-owned metadata output '" +
                                         path.string() + "': " + error.message());
            }
            removed += didRemove ? 1 : 0;
        }
    }
    return removed;
}

void RequireNamespaceAvailableWithoutManifest(const Parameters& params,
                                              const OutputLayout& layout) {
    if (!std::filesystem::exists(layout.root)) {
        return;
    }
    for (const std::filesystem::directory_entry& entry :
         std::filesystem::directory_iterator(layout.root)) {
        const std::string filename = entry.path().filename().string();
        const std::string normalizedFilename = LowerAscii(filename);
        const std::string normalizedPvd = LowerAscii(params.outputBasename + ".pvd");
        const std::string normalizedSummary =
            LowerAscii(params.outputBasename + "_run_summary.txt");
        const bool metadataCollision =
            normalizedFilename == normalizedPvd || normalizedFilename == normalizedSummary;
        if (metadataCollision || OutputFrameIndex(entry.path(), params.outputBasename).has_value()) {
            throw std::runtime_error(
                "Output namespace already contains files without a valid ownership manifest; "
                "refusing to overwrite: " + entry.path().string());
        }
    }
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
    const std::string fileName = LowerAscii(path.filename().string());
    const std::string prefix = LowerAscii(basename) + "_f";
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
                           ResumeManifest* manifest,
                           const OutputLayout& layout,
                           const std::chrono::steady_clock::time_point& progressStart) {
    if (sequence == nullptr || manifest == nullptr) {
        throw std::runtime_error("Internal error: null sequence/manifest for resume plan.");
    }
    if (manifest->frames.size() != sequence->frames.size() ||
        manifest->outputChecksums.size() != sequence->frames.size() ||
        manifest->outputFiles.size() != sequence->frames.size()) {
        throw std::runtime_error("Resume manifest frame count does not match the validated input sequence.");
    }

    ResumePlan plan;
    plan.totalWork = static_cast<vtkIdType>(sequence->frames.size());
    std::cout << "[Resume] Manifest-bound output resume requested; checksums and VTK invariants will be validated.\n";

    std::size_t candidateCount = 0;
    for (; candidateCount < sequence->frames.size(); ++candidateCount) {
        if (manifest->outputChecksums[candidateCount].empty()) {
            break;
        }
    }

    const auto validationStage = LogStageStarted("Restart validation",
                                                 "candidate_frames=" + std::to_string(candidateCount));
    std::vector<ValidatedResumeFrame> validatedFrames;
    validatedFrames.reserve(candidateCount);
    for (std::size_t i = 0; i < candidateCount; ++i) {
        const std::filesystem::path& outputPath = layout.frames.at(i);
        std::cout << "[Resume] Validating existing frame " << i << ": "
                  << outputPath.filename().string() << "\n";
        try {
            if (!std::filesystem::is_regular_file(outputPath)) {
                throw std::runtime_error("Manifest-bound output file is missing or not regular: " +
                                         outputPath.string());
            }
            const std::string checksum = integrity::Sha256File(outputPath);
            if (checksum != manifest->outputChecksums[i]) {
                throw std::runtime_error("Existing resume output checksum does not match manifest: " +
                                         outputPath.string());
            }
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
            break;
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
        std::cout << "[Resume] Existing frame " << i
                  << " passed checksum, image type, dimension, and bounds validation.\n";
    }
    if (validatedFrames.size() == candidateCount) {
        LogStageFinished("Restart validation",
                         validationStage,
                         "usable_prefix=" + std::to_string(validatedFrames.size()));
    }

    plan.startFrame = validatedFrames.size();
    plan.completedWork = static_cast<vtkIdType>(plan.startFrame);
    plan.reusedResults.reserve(plan.startFrame);
    plan.reusedPvdEntries.reserve(plan.startFrame);
    for (std::size_t i = 0; i < plan.startFrame; ++i) {
        plan.reusedResults.push_back(validatedFrames[i].result);
        plan.reusedPvdEntries.push_back({validatedFrames[i].result.time, validatedFrames[i].result.relativeFile});
    }

    const int removedRestartFiles = RemoveManifestOwnedOutputs(
        *manifest, layout.root, plan.startFrame, true);
    for (std::size_t i = plan.startFrame; i < manifest->outputChecksums.size(); ++i) {
        manifest->outputChecksums[i].clear();
    }
    manifest->pvdChecksum.clear();
    manifest->summaryChecksum.clear();
    manifest->state = "in_progress";
    WriteResumeManifest(layout, *manifest);
    std::cout << "[Resume] Removed only manifest-owned outputs from frame " << plan.startFrame
              << " onward; removed_files=" << removedRestartFiles << ".\n";
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
    if (params.method == InterpolationMethod::Kriging &&
        (params.kriging.neighbors < 1 ||
         params.kriging.neighbors > kMaximumKrigingNeighbors)) {
        throw std::runtime_error(
            "Kriging neighbors must be between 1 and " +
            std::to_string(kMaximumKrigingNeighbors) + ".");
    }
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

int RunSolver(const Parameters& unvalidatedParams) {
    Parameters params = unvalidatedParams;
    ValidateParameters(&params);
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
    const auto discoverStage = LogStageStarted("Data discovery", params.inputPath.string());
    SequenceInfo sequence = CollectInputSequence(params);
    LogStageFinished("Data discovery",
                     discoverStage,
                     "frames=" + std::to_string(sequence.frames.size()) +
                         (sequence.boundsInitialized
                              ? ", bounds=" + FormatBounds(sequence.bounds)
                              : ", bounds=streaming_per_frame"));

    const auto identityStage = LogStageStarted(
        "Input executable preflight",
        "frames=" + std::to_string(sequence.frames.size()));
    const RunIdentity runIdentity = BuildRunIdentity(params, sequence);
    LogStageFinished("Input executable preflight",
                     identityStage,
                     "run_sha256=" + runIdentity.signature);

    const OutputLayout outputLayout = BuildOutputLayout(params, sequence.frames.size());
    const std::optional<ResumeManifest> existingManifest =
        ReadResumeManifest(outputLayout, params.outputBasename);
    std::filesystem::create_directories(outputLayout.root);
    ResumeManifest manifest = MakeResumeManifest(params, runIdentity);

    std::vector<FrameResult> results;
    results.reserve(sequence.frames.size());
    std::vector<std::pair<double, std::filesystem::path>> pvdEntries;
    pvdEntries.reserve(sequence.frames.size());

    const auto progressStart = std::chrono::steady_clock::now();
    ResumePlan resumePlan;
    const bool reusableManifest = resumeExistingOutput && existingManifest.has_value() &&
                                  existingManifest->runSignature == runIdentity.signature &&
                                  existingManifest->parameterChecksum == runIdentity.parameterChecksum &&
                                  SameFrameInputIdentities(existingManifest->frames,
                                                           runIdentity.frames);
    if (reusableManifest) {
        manifest = *existingManifest;
        resumePlan = BuildResumePlan(params, &sequence, &manifest, outputLayout, progressStart);
        results = resumePlan.reusedResults;
        pvdEntries = resumePlan.reusedPvdEntries;
    } else {
        const auto cleanupStage = LogStageStarted("Output ownership preparation",
                                                  outputLayout.root.string());
        int removedOwnedFiles = 0;
        if (existingManifest.has_value()) {
            removedOwnedFiles = RemoveManifestOwnedOutputs(
                *existingManifest, outputLayout.root, 0, true);
            if (resumeExistingOutput) {
                std::cout << "[Resume] Cached source/configuration identity changed; "
                             "reuse was rejected and only manifest-owned outputs were removed.\n";
            }
        }
        RequireNamespaceAvailableWithoutManifest(params, outputLayout);
        WriteResumeManifest(outputLayout, manifest);
        LogStageFinished("Output ownership preparation",
                         cleanupStage,
                         "removed_manifest_owned_files=" + std::to_string(removedOwnedFiles));
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
        const FrameInputIdentity executionIdentity =
            ComputeFrameInputIdentityFromDataSet(params, sequence, i, input, false);
        if (!SameFrameInputIdentity(executionIdentity, runIdentity.frames.at(i))) {
            throw std::runtime_error("Input frame changed after executable preflight; refusing "
                                     "to publish output for frame " + std::to_string(i) + ".");
        }
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

        const std::string& fileName = manifest.outputFiles.at(i);
        const auto writeStart = std::chrono::steady_clock::now();
        WriteVti(outputLayout.frames.at(i), output);
        manifest.outputChecksums.at(i) = integrity::Sha256File(outputLayout.frames.at(i));
        manifest.state = "in_progress";
        WriteResumeManifest(outputLayout, manifest);
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

            const std::filesystem::path& outputPath = outputLayout.frames.at(i);
            std::error_code removeError;
            std::filesystem::remove(outputPath, removeError);
            if (removeError) {
                throw std::runtime_error("Frame failure cleanup could not remove manifest-owned output '" +
                                         outputPath.string() + "': " + removeError.message() +
                                         "; original error: " + lastError);
            }
            manifest.outputChecksums.at(i).clear();
            manifest.state = "in_progress";
            try {
                WriteResumeManifest(outputLayout, manifest);
            } catch (const std::exception& manifestError) {
                throw std::runtime_error("Frame cleanup could not update its ownership manifest: " +
                                         std::string(manifestError.what()));
            }
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

    const auto writeStage = LogStageStarted("Result writing", outputLayout.root.string());
    WritePvd(outputLayout.pvd, pvdEntries);
    manifest.pvdChecksum = integrity::Sha256File(outputLayout.pvd);
    manifest.state = "in_progress";
    WriteResumeManifest(outputLayout, manifest);
    const double totalSeconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - totalStart).count();
    WriteTextFile(outputLayout.summary,
                  BuildRunSummary(params, sequence, results, totalSeconds));
    manifest.summaryChecksum = integrity::Sha256File(outputLayout.summary);
    manifest.state = "complete";
    WriteResumeManifest(outputLayout, manifest);
    LogStageFinished("Result writing", writeStage, outputLayout.pvd.filename().string());

    std::cout << "\n[Stage] Solver finished. total=" << FormatSeconds(totalSeconds)
              << "s, output=" << outputLayout.pvd.string() << "\n";
    return 0;
}

}  // namespace streamcenterplus::resample
