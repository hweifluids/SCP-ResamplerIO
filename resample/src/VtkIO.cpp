#include "streamcenterplus_resample/Core.h"
#include <streamcenterplus/io/DataSetReader.h>
#include <streamcenterplus/VtkHdfIO.h>

#include <vtkCellData.h>
#include <vtkCellDataToPointData.h>
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDataSetAttributes.h>
#include <vtkFieldData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkImageData.h>
#include <vtkIdList.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLGenericDataObjectReader.h>
#include <vtkXMLImageDataWriter.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace streamcenterplus::resample {

namespace {

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool IsSupportedDataFile(const std::filesystem::path& path) {
    const std::string ext = ToLower(streamcenterplus::vtkhdf::Extension(path));
    return ext == ".vtk" || ext == ".vti" || ext == ".vtr" ||
           ext == ".vts" || ext == ".vtu" || ext == ".vtp" ||
           ext == ".vtkhdf" || ext == ".nek5000" || ext == ".nek" ||
           ext == ".case" || ext == ".foam" || ext == ".scpm";
}

DataSetFrameKind FromIoFrameKind(streamcenterplus::io::DataSetFrameKind kind) {
    switch (kind) {
    case streamcenterplus::io::DataSetFrameKind::File:
        return DataSetFrameKind::File;
    case streamcenterplus::io::DataSetFrameKind::Nek5000TimeStep:
        return DataSetFrameKind::Nek5000TimeStep;
    case streamcenterplus::io::DataSetFrameKind::EnSightGoldTimeStep:
        return DataSetFrameKind::EnSightGoldTimeStep;
    case streamcenterplus::io::DataSetFrameKind::OpenFoamTimeStep:
        return DataSetFrameKind::OpenFoamTimeStep;
    case streamcenterplus::io::DataSetFrameKind::MultiZoneManifestTimeStep:
        return DataSetFrameKind::MultiZoneManifestTimeStep;
    }
    throw std::runtime_error("Unknown shared dataset frame kind.");
}

streamcenterplus::io::DataSetFrame ToIoFrame(const FrameInfo& frame) {
    streamcenterplus::io::DataSetFrame out;
    out.path = frame.path;
    out.time = frame.time;
    out.readerTime = frame.readerTime;
    out.coordinatePath = frame.coordinatePath;
    out.sourceIndex = frame.sourceIndex;
    out.coordinateSourceIndex = frame.coordinateSourceIndex;
    out.multiZoneBackgroundIndex = frame.multiZoneBackgroundIndex;
    switch (frame.kind) {
    case DataSetFrameKind::File:
        out.kind = streamcenterplus::io::DataSetFrameKind::File;
        break;
    case DataSetFrameKind::Nek5000TimeStep:
        out.kind = streamcenterplus::io::DataSetFrameKind::Nek5000TimeStep;
        break;
    case DataSetFrameKind::EnSightGoldTimeStep:
        out.kind = streamcenterplus::io::DataSetFrameKind::EnSightGoldTimeStep;
        break;
    case DataSetFrameKind::OpenFoamTimeStep:
        out.kind = streamcenterplus::io::DataSetFrameKind::OpenFoamTimeStep;
        break;
    case DataSetFrameKind::MultiZoneManifestTimeStep:
        out.kind = streamcenterplus::io::DataSetFrameKind::MultiZoneManifestTimeStep;
        break;
    }
    return out;
}

std::string XmlEscape(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size());
    for (const char ch : value) {
        switch (ch) {
        case '&':
            escaped += "&amp;";
            break;
        case '<':
            escaped += "&lt;";
            break;
        case '>':
            escaped += "&gt;";
            break;
        case '"':
            escaped += "&quot;";
            break;
        case '\'':
            escaped += "&apos;";
            break;
        default:
            escaped += ch;
            break;
        }
    }
    return escaped;
}

std::string XmlUnescape(std::string value) {
    const std::pair<const char*, const char*> replacements[] = {
        {"&quot;", "\""},
        {"&apos;", "'"},
        {"&lt;", "<"},
        {"&gt;", ">"},
        {"&amp;", "&"},
    };
    for (const auto& replacement : replacements) {
        std::string::size_type pos = 0;
        while ((pos = value.find(replacement.first, pos)) != std::string::npos) {
            value.replace(pos, std::strlen(replacement.first), replacement.second);
            pos += std::strlen(replacement.second);
        }
    }
    return value;
}

std::string AttributeValue(const std::string& tag, const std::string& name) {
    const std::regex pattern(name + R"(\s*=\s*(['"])(.*?)\1)", std::regex::icase);
    std::smatch match;
    if (!std::regex_search(tag, match, pattern)) {
        return {};
    }
    return XmlUnescape(match[2].str());
}

std::vector<FrameInfo> CollectFolderFiles(const std::filesystem::path& folder, double physicalDt) {
    if (!std::filesystem::is_directory(folder)) {
        throw std::runtime_error("Input folder does not exist: " + folder.string());
    }

    struct TaggedFrame {
        long long numericTag = -1;
        std::string fileName;
        std::size_t localIndex = 0;
        FrameInfo frame;
    };

    std::regex numberRegex(R"((\d+))");
    std::vector<TaggedFrame> tagged;
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
        if (!entry.is_regular_file() || !IsSupportedDataFile(entry.path())) {
            continue;
        }
        long long tag = -1;
        const std::string filename = entry.path().filename().string();
        for (std::sregex_iterator it(filename.begin(), filename.end(), numberRegex), end; it != end; ++it) {
            tag = std::stoll((*it)[1].str());
        }
        const std::vector<FrameInfo> fileFrames = EnumerateDataSetFrames(entry.path(), physicalDt);
        for (std::size_t localIndex = 0; localIndex < fileFrames.size(); ++localIndex) {
            tagged.push_back({tag, filename, localIndex, fileFrames[localIndex]});
        }
    }
    std::sort(tagged.begin(), tagged.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.numericTag != rhs.numericTag) {
            return lhs.numericTag < rhs.numericTag;
        }
        if (lhs.fileName != rhs.fileName) {
            return lhs.fileName < rhs.fileName;
        }
        return lhs.localIndex < rhs.localIndex;
    });

    std::vector<FrameInfo> frames;
    frames.reserve(tagged.size());
    for (std::size_t i = 0; i < tagged.size(); ++i) {
        FrameInfo frame = std::move(tagged[i].frame);
        if (!frame.readerTime.has_value() &&
            frame.kind != DataSetFrameKind::MultiZoneManifestTimeStep) {
            frame.time = static_cast<double>(i) * physicalDt;
        }
        frames.push_back(std::move(frame));
    }
    return frames;
}

std::vector<FrameInfo> CollectPvdFiles(const std::filesystem::path& pvdFile, double physicalDt) {
    std::ifstream input(pvdFile);
    if (!input) {
        throw std::runtime_error("Cannot open PVD file: " + pvdFile.string());
    }
    std::ostringstream text;
    text << input.rdbuf();

    const std::regex dataSetTag(R"(<\s*DataSet\b[^>]*>)", std::regex::icase);
    const std::string xml = text.str();
    std::vector<FrameInfo> frames;
    for (std::sregex_iterator it(xml.begin(), xml.end(), dataSetTag), end; it != end; ++it) {
        const std::string tag = it->str();
        std::string file = AttributeValue(tag, "file");
        if (file.empty()) {
            continue;
        }
        std::string timestep = AttributeValue(tag, "timestep");
        double time = static_cast<double>(frames.size()) * physicalDt;
        if (!timestep.empty()) {
            try {
                time = std::stod(timestep);
            } catch (...) {
                throw std::runtime_error("Invalid PVD timestep '" + timestep + "' in " + pvdFile.string());
            }
        }
        std::filesystem::path framePath(file);
        if (framePath.is_relative()) {
            framePath = pvdFile.parent_path() / framePath;
        }
        FrameInfo frame{framePath, time};
        if (ToLower(streamcenterplus::vtkhdf::Extension(framePath)) == ".nek5000") {
            frame.kind = DataSetFrameKind::Nek5000TimeStep;
            frame.readerTime = time;
        }
        frames.push_back(std::move(frame));
    }
    return frames;
}

bool IsFiniteBounds(const double bounds[6]) {
    for (int i = 0; i < 6; ++i) {
        if (!std::isfinite(bounds[i])) {
            return false;
        }
    }
    return bounds[0] <= bounds[1] && bounds[2] <= bounds[3] && bounds[4] <= bounds[5];
}

std::array<double, 6> BoundsFromDataObject(vtkDataObject* dataObject,
                                           const std::filesystem::path& filePath) {
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(dataObject);
    if (dataSet == nullptr) {
        throw std::runtime_error("The VTK file is not a vtkDataSet: " + filePath.string());
    }
    if (dataSet->GetNumberOfPoints() <= 0) {
        throw std::runtime_error("The VTK dataset has no points: " + filePath.string());
    }
    std::array<double, 6> bounds = {};
    dataSet->GetBounds(bounds.data());
    if (!IsFiniteBounds(bounds.data())) {
        throw std::runtime_error("Input data has invalid bounds: " + filePath.string());
    }
    return bounds;
}

vtkSmartPointer<vtkDataSet> RetainDataSet(vtkDataObject* dataObject,
                                          const std::filesystem::path& sourcePath) {
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(dataObject);
    if (dataSet == nullptr) {
        throw std::runtime_error("The input source did not produce a vtkDataSet: " +
                                 sourcePath.string());
    }
    if (dataSet->GetNumberOfPoints() <= 0) {
        throw std::runtime_error("The input dataset has no points: " + sourcePath.string());
    }

    vtkSmartPointer<vtkDataSet> retained;
    retained.TakeReference(vtkDataSet::SafeDownCast(dataSet->NewInstance()));
    if (retained == nullptr) {
        throw std::runtime_error("Cannot retain the VTK dataset produced by: " + sourcePath.string());
    }
    retained->ShallowCopy(dataSet);
    return retained;
}

void ExpandBounds(const double frameBounds[6], std::array<double, 6>* bounds, bool* initialized) {
    if (!IsFiniteBounds(frameBounds)) {
        throw std::runtime_error("Input data has invalid bounds.");
    }
    if (bounds == nullptr || initialized == nullptr) {
        return;
    }
    if (!*initialized) {
        for (int i = 0; i < 6; ++i) {
            (*bounds)[i] = frameBounds[i];
        }
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

void EnsureNonzeroBounds(std::array<double, 6>* bounds) {
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

bool ShouldPrecomputeSeriesBounds(const Parameters& params) {
    return params.boundsMode == BoundsMode::Global &&
           params.method != InterpolationMethod::VtkResampleToImage;
}

void ApplyFrameRange(const Parameters& params, std::vector<FrameInfo>* frames) {
    if (frames == nullptr || frames->empty()) {
        return;
    }
    const int totalFrames = static_cast<int>(frames->size());
    int startFrame = params.frameRange.start.value_or(0);
    int endFrame = params.frameRange.end.value_or(totalFrames - 1);
    startFrame = std::max(0, startFrame);
    endFrame = std::min(endFrame, totalFrames - 1);
    if (startFrame > endFrame) {
        throw std::runtime_error("Invalid resample frame range after clamping to the available input sequence.");
    }

    std::vector<FrameInfo> selected;
    selected.reserve(static_cast<std::size_t>(endFrame - startFrame + 1));
    for (int i = startFrame; i <= endFrame; ++i) {
        selected.push_back((*frames)[static_cast<std::size_t>(i)]);
    }
    *frames = std::move(selected);
}

}  // namespace

std::vector<FrameInfo> EnumerateDataSetFrames(const std::filesystem::path& sourcePath,
                                              double physicalDt) {
    const std::vector<streamcenterplus::io::DataSetFrame> sharedFrames =
        streamcenterplus::io::EnumerateDataSetFrames(sourcePath, physicalDt);
    std::vector<FrameInfo> frames;
    frames.reserve(sharedFrames.size());
    for (const auto& shared : sharedFrames) {
        FrameInfo frame;
        frame.path = shared.path;
        frame.time = shared.time;
        frame.kind = FromIoFrameKind(shared.kind);
        frame.readerTime = shared.readerTime;
        frame.coordinatePath = shared.coordinatePath;
        frame.sourceIndex = shared.sourceIndex;
        frame.coordinateSourceIndex = shared.coordinateSourceIndex;
        frame.multiZoneBackgroundIndex = shared.multiZoneBackgroundIndex;
        frames.push_back(std::move(frame));
    }
    return frames;
}

std::array<double, 6> ResolveCropBounds(const Parameters& params,
                                        const std::array<double, 6>& inputBounds) {
    std::array<double, 6> bounds = inputBounds;
    for (int axis = 0; axis < 3; ++axis) {
        const AxisRange& range = params.cropRanges[static_cast<std::size_t>(axis)];
        const int lo = axis * 2;
        const int hi = lo + 1;
        if (range.min.has_value()) {
            bounds[lo] = std::max(bounds[lo], *range.min);
        }
        if (range.max.has_value()) {
            bounds[hi] = std::min(bounds[hi], *range.max);
        }
        if (bounds[hi] < bounds[lo]) {
            throw std::runtime_error("Configured resample crop range does not overlap the input bounds.");
        }
    }
    return bounds;
}

std::array<double, 6> ResolveConfiguredBounds(const Parameters& params,
                                              const std::array<double, 6>& inputBounds) {
    std::array<double, 6> bounds = ResolveCropBounds(params, inputBounds);
    for (int axis = 0; axis < 3; ++axis) {
        const int lo = axis * 2;
        const int hi = lo + 1;
        const double scale = params.inputScale[static_cast<std::size_t>(axis)];
        bounds[lo] *= scale;
        bounds[hi] *= scale;
        if (bounds[hi] < bounds[lo]) {
            std::swap(bounds[lo], bounds[hi]);
        }
    }
    EnsureNonzeroBounds(&bounds);
    return bounds;
}

SequenceInfo CollectInputSequence(const Parameters& params) {
    std::vector<FrameInfo> frames =
        EnumerateDataSetFrames(params.inputPath, params.physicalDt);

    if (frames.empty()) {
        throw std::runtime_error("No supported dataset frames were found in input: " +
                                 params.inputPath.string());
    }
    ApplyFrameRange(params, &frames);
    for (const FrameInfo& frame : frames) {
        if (frame.kind != DataSetFrameKind::File) {
            continue;
        }
        const std::filesystem::path realPath = streamcenterplus::vtkhdf::RealPath(frame.path);
        if (!std::filesystem::exists(realPath)) {
            throw std::runtime_error("Referenced VTK frame does not exist: " + frame.path.string());
        }
        if (!IsSupportedDataFile(frame.path)) {
            throw std::runtime_error("Unsupported VTK frame extension: " + frame.path.string());
        }
    }

    SequenceInfo sequence;
    sequence.frames = std::move(frames);
    if (ShouldPrecomputeSeriesBounds(params)) {
        for (const FrameInfo& frame : sequence.frames) {
            const std::array<double, 6> frameBounds = ReadDataSetBounds(frame);
            ExpandBounds(frameBounds.data(), &sequence.bounds, &sequence.boundsInitialized);
        }
        EnsureNonzeroBounds(&sequence.bounds);
        sequence.bounds = ResolveConfiguredBounds(params, sequence.bounds);
        sequence.boundsInitialized = true;
        sequence.seriesBoundsPrecomputed = true;
    }
    return sequence;
}

std::array<double, 6> ReadDataSetBounds(const std::filesystem::path& filePath) {
    const std::string extension =
        ToLower(streamcenterplus::vtkhdf::Extension(filePath));
    if (extension == ".nek5000" || extension == ".nek" ||
        extension == ".case" || extension == ".foam" || extension == ".scpm") {
        const std::vector<FrameInfo> frames = EnumerateDataSetFrames(filePath);
        return ReadDataSetBounds(frames.front());
    }
    return ReadDataSetBounds(FrameInfo{filePath, 0.0});
}

std::array<double, 6> ReadDataSetBounds(const FrameInfo& frame) {
    const std::filesystem::path& filePath = frame.path;
    const std::string ext = ToLower(streamcenterplus::vtkhdf::Extension(filePath));
    if (frame.kind != DataSetFrameKind::File || ext == ".nek5000" ||
        ext == ".nek" ||
        ext == ".case" || ext == ".foam" || ext == ".scpm") {
        if (!frame.readerTime.has_value() &&
            frame.kind != DataSetFrameKind::MultiZoneManifestTimeStep) {
            const std::vector<FrameInfo> frames = EnumerateDataSetFrames(filePath);
            return ReadDataSetBounds(frames.front());
        }
        streamcenterplus::io::ReadOptions options;
        options.readAllArrays = false;
        options.attachSourceBlockIds = false;
        options.attachNek5000SpectralMetadata = false;
        options.buildNek5000LinearizedTopology = false;
        options.attachNek5000QuadratureArrays = false;
        return BoundsFromDataObject(
            streamcenterplus::io::ReadDataSet(ToIoFrame(frame), options),
            filePath);
    }
    if (ext == ".vtkhdf") {
        return BoundsFromDataObject(streamcenterplus::vtkhdf::ReadDataSet(filePath), filePath);
    }

    if (ext == ".vtk") {
        auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(filePath.string().c_str());
        reader->ReadAllScalarsOff();
        reader->ReadAllVectorsOff();
        reader->ReadAllNormalsOff();
        reader->ReadAllTensorsOff();
        reader->ReadAllColorScalarsOff();
        reader->ReadAllTCoordsOff();
        reader->ReadAllFieldsOff();
        reader->SetScalarsName("__streamcenterplus_bounds_only_no_scalars__");
        reader->SetVectorsName("__streamcenterplus_bounds_only_no_vectors__");
        reader->SetNormalsName("__streamcenterplus_bounds_only_no_normals__");
        reader->SetTensorsName("__streamcenterplus_bounds_only_no_tensors__");
        reader->SetTCoordsName("__streamcenterplus_bounds_only_no_tcoords__");
        reader->SetFieldDataName("__streamcenterplus_bounds_only_no_fields__");
        reader->Update();
        return BoundsFromDataObject(reader->GetOutput(), filePath);
    }

    auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(filePath.string().c_str());
    reader->UpdateInformation();
    if (vtkDataArraySelection* selection = reader->GetPointDataArraySelection()) {
        selection->DisableAllArrays();
    }
    if (vtkDataArraySelection* selection = reader->GetCellDataArraySelection()) {
        selection->DisableAllArrays();
    }
    if (vtkDataArraySelection* selection = reader->GetColumnArraySelection()) {
        selection->DisableAllArrays();
    }
    reader->Update();
    return BoundsFromDataObject(reader->GetOutputDataObject(0), filePath);
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const std::filesystem::path& filePath) {
    const std::string extension =
        ToLower(streamcenterplus::vtkhdf::Extension(filePath));
    if (extension == ".nek5000" || extension == ".nek" ||
        extension == ".case" || extension == ".foam" || extension == ".scpm") {
        const std::vector<FrameInfo> frames = EnumerateDataSetFrames(filePath);
        return ReadDataSet(frames.front());
    }
    return ReadDataSet(FrameInfo{filePath, 0.0});
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const FrameInfo& frame) {
    const std::filesystem::path& filePath = frame.path;
    const std::string ext = ToLower(streamcenterplus::vtkhdf::Extension(filePath));
    if (frame.kind != DataSetFrameKind::File || ext == ".nek5000" ||
        ext == ".nek" ||
        ext == ".case" || ext == ".foam" || ext == ".scpm") {
        if (!frame.readerTime.has_value() &&
            frame.kind != DataSetFrameKind::MultiZoneManifestTimeStep) {
            const std::vector<FrameInfo> frames = EnumerateDataSetFrames(filePath);
            return ReadDataSet(frames.front());
        }
        streamcenterplus::io::ReadOptions options;
        if (frame.kind == DataSetFrameKind::Nek5000TimeStep ||
            ext == ".nek5000" || ext == ".nek") {
            // Spectral resampling uses the element-local GLL coordinates and
            // its own AABB/Newton locator.  Tens of millions of visualization
            // subcells and three quadrature arrays are therefore unnecessary.
            options.buildNek5000LinearizedTopology = false;
            options.attachNek5000QuadratureArrays = false;
        }
        return streamcenterplus::io::ReadDataSet(ToIoFrame(frame), options);
    }

    if (ext == ".vtkhdf") {
        return streamcenterplus::vtkhdf::ReadDataSet(filePath);
    }

    if (ext == ".vtk") {
        auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
        reader->SetFileName(filePath.string().c_str());
        reader->ReadAllScalarsOn();
        reader->ReadAllVectorsOn();
        reader->ReadAllTensorsOn();
        reader->ReadAllFieldsOn();
        reader->Update();
        return RetainDataSet(reader->GetOutput(), filePath);
    }

    auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
    reader->SetFileName(filePath.string().c_str());
    reader->Update();
    return RetainDataSet(reader->GetOutputDataObject(0), filePath);
}

vtkSmartPointer<vtkDataSet> PreparePointData(vtkDataSet* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal error: null input dataset.");
    }

    vtkSmartPointer<vtkDataSet> copy;
    copy.TakeReference(input->NewInstance());
    copy->ShallowCopy(input);

    vtkPointData* pointData = copy->GetPointData();
    vtkCellData* cellData = copy->GetCellData();
    const bool multiZoneCanonical =
        pointData != nullptr &&
        pointData->GetArray(streamcenterplus::io::kMultiZoneOwnershipArray) != nullptr;
    if (multiZoneCanonical && pointData != nullptr && cellData != nullptr) {
        // The reader's cell cut is categorical.  Promote it conservatively
        // before the generic cell-to-point filter can average/truncate it.
        vtkDataArray* cellCut = cellData->GetArray(
            streamcenterplus::io::kMultiZoneBackgroundCutArray);
        if (cellCut != nullptr && cellCut->GetNumberOfComponents() == 1 &&
            cellCut->GetNumberOfTuples() >= copy->GetNumberOfCells()) {
            vtkNew<vtkUnsignedCharArray> pointCut;
            pointCut->SetName(streamcenterplus::io::kMultiZoneBackgroundCutArray);
            pointCut->SetNumberOfTuples(copy->GetNumberOfPoints());
            pointCut->FillValue(0);
            vtkNew<vtkIdList> pointIds;
            for (vtkIdType cellId = 0; cellId < copy->GetNumberOfCells(); ++cellId) {
                if (cellCut->GetComponent(cellId, 0) == 0.0) {
                    continue;
                }
                copy->GetCellPoints(cellId, pointIds);
                for (vtkIdType pointIndex = 0; pointIndex < pointIds->GetNumberOfIds();
                     ++pointIndex) {
                    pointCut->SetValue(pointIds->GetId(pointIndex), 1);
                }
            }
            pointData->RemoveArray(streamcenterplus::io::kMultiZoneBackgroundCutArray);
            pointData->AddArray(pointCut);
            cellData->RemoveArray(streamcenterplus::io::kMultiZoneBackgroundCutArray);
        }

        // The mesh-set reader deliberately composes only point-associated
        // physical fields.  Any remaining CellData is zone-local/background
        // diagnostic data: promoting it here would manufacture a PointData
        // field with stale background values in the foreground-covered
        // region (including when the same-named PointData was correctly
        // rejected by the reader as unavailable in the foreground).  The
        // canonical point diagnostics above are the only cell diagnostics
        // that Resample exposes for a two-zone input.
        cellData->Initialize();
        return copy;
    }
    if (pointData != nullptr && cellData != nullptr) {
        std::vector<std::string> conflictingNames;
        for (int i = 0; i < cellData->GetNumberOfArrays(); ++i) {
            vtkAbstractArray* cellArray = cellData->GetAbstractArray(i);
            const char* name = cellArray != nullptr ? cellArray->GetName() : nullptr;
            if (name != nullptr && std::string(name) != vtkDataSetAttributes::GhostArrayName() &&
                pointData->GetAbstractArray(name) != nullptr) {
                conflictingNames.emplace_back(name);
            }
        }
        for (const std::string& name : conflictingNames) {
            if (multiZoneCanonical) {
                // Reader-composed point data is already in the common world
                // reference frame.  A background-only cell array with the
                // same name must not replace it during point conversion.
                cellData->RemoveArray(name.c_str());
            } else {
                pointData->RemoveArray(name.c_str());
            }
        }
    }

    if (cellData == nullptr || cellData->GetNumberOfArrays() == 0) {
        return copy;
    }

    auto converter = vtkSmartPointer<vtkCellDataToPointData>::New();
    converter->SetInputData(copy);
    converter->PassCellDataOff();
    converter->ProcessAllArraysOn();
    converter->SetContributingCellOption(vtkCellDataToPointData::All);
    converter->Update();

    vtkSmartPointer<vtkDataSet> output = converter->GetOutput();
    output->GetCellData()->Initialize();
    return output;
}

void WriteVti(const std::filesystem::path& outputPath, vtkImageData* image) {
    if (image == nullptr) {
        throw std::runtime_error("Internal error: null output image.");
    }
    std::filesystem::create_directories(outputPath.parent_path());
    auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName(outputPath.string().c_str());
    writer->SetInputData(image);
    writer->SetDataModeToBinary();
    if (!writer->Write()) {
        throw std::runtime_error("Failed to write VTI file: " + outputPath.string());
    }
}

void WritePvd(const std::filesystem::path& outputPath,
              const std::vector<std::pair<double, std::filesystem::path>>& frames) {
    std::filesystem::create_directories(outputPath.parent_path());
    std::ofstream output(outputPath);
    if (!output) {
        throw std::runtime_error("Failed to write PVD file: " + outputPath.string());
    }
    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output << "  <Collection>\n";
    for (const auto& frame : frames) {
        output << "    <DataSet timestep=\"" << std::setprecision(17) << frame.first
               << "\" group=\"\" part=\"0\" file=\"" << XmlEscape(frame.second.generic_string()) << "\"/>\n";
    }
    output << "  </Collection>\n";
    output << "</VTKFile>\n";
}

void WriteTextFile(const std::filesystem::path& outputPath, const std::string& text) {
    std::filesystem::create_directories(outputPath.parent_path());
    std::ofstream output(outputPath);
    if (!output) {
        throw std::runtime_error("Failed to write text file: " + outputPath.string());
    }
    output << text;
}

}  // namespace streamcenterplus::resample
