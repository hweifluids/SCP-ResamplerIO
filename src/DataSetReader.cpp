#include "streamcenterplus/io/DataSetReader.h"

#include <streamcenterplus/VtkHdfIO.h>

#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericCell.h>
#include <vtkGenericEnSightReader.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkOpenFOAMReader.h>
#include <vtkPOpenFOAMReader.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkSMPTools.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkStaticCellLocator.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLGenericDataObjectReader.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace streamcenterplus::io {

namespace {

constexpr const char* kSourceBlockIdArray = "SCP_SourceBlockId";
constexpr const char* kTimeValueArray = "TimeValue";

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string LowerExtension(const std::filesystem::path& path) {
    return ToLower(streamcenterplus::vtkhdf::Extension(path));
}

std::string TrimAscii(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string StripComment(const std::string& line) {
    const std::size_t hash = line.find('#');
    return TrimAscii(hash == std::string::npos ? line : line.substr(0, hash));
}

bool IsFiniteTime(double value) {
    return std::isfinite(value);
}

bool IsRegularFileWithExtension(const std::filesystem::path& path,
                                const std::string& extension) {
    return std::filesystem::is_regular_file(path) && LowerExtension(path) == extension;
}

bool IsOrdinaryVtkExtension(const std::string& extension) {
    return extension == ".vtk" || extension == ".vti" || extension == ".vtr" ||
           extension == ".vts" || extension == ".vtu" || extension == ".vtp" ||
           extension == ".vtm" || extension == ".vtkhdf";
}

vtkSmartPointer<vtkDataSet> CloneDataSet(vtkDataSet* input) {
    if (input == nullptr) {
        return nullptr;
    }
    vtkSmartPointer<vtkDataSet> output;
    output.TakeReference(vtkDataSet::SafeDownCast(input->NewInstance()));
    if (output == nullptr) {
        throw std::runtime_error("Cannot clone a VTK dataset instance.");
    }
    output->DeepCopy(input);
    return output;
}

void CollectDataSets(vtkDataObject* object,
                     std::vector<vtkSmartPointer<vtkDataSet>>* dataSets) {
    if (object == nullptr || dataSets == nullptr) {
        return;
    }
    if (vtkDataSet* dataSet = vtkDataSet::SafeDownCast(object)) {
        if (dataSet->GetNumberOfPoints() > 0) {
            dataSets->push_back(CloneDataSet(dataSet));
        }
        return;
    }
    vtkCompositeDataSet* composite = vtkCompositeDataSet::SafeDownCast(object);
    if (composite == nullptr) {
        return;
    }
    vtkSmartPointer<vtkCompositeDataIterator> iterator;
    iterator.TakeReference(composite->NewIterator());
    if (iterator == nullptr) {
        return;
    }
    iterator->SkipEmptyNodesOn();
    for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal();
         iterator->GoToNextItem()) {
        CollectDataSets(iterator->GetCurrentDataObject(), dataSets);
    }
}

void AttachSourceBlockId(vtkDataSet* dataSet, int blockId) {
    if (dataSet == nullptr) {
        return;
    }
    vtkNew<vtkIntArray> cells;
    cells->SetName(kSourceBlockIdArray);
    cells->SetNumberOfTuples(dataSet->GetNumberOfCells());
    cells->FillValue(blockId);
    dataSet->GetCellData()->AddArray(cells);
}

vtkSmartPointer<vtkDataSet> CollapseDataObject(vtkDataObject* object,
                                               const ReadOptions& options,
                                               const std::filesystem::path& sourcePath) {
    std::vector<vtkSmartPointer<vtkDataSet>> dataSets;
    CollectDataSets(object, &dataSets);
    if (dataSets.empty()) {
        throw std::runtime_error("The source produced no non-empty VTK dataset: " +
                                 sourcePath.string());
    }
    if (dataSets.size() == 1) {
        if (options.attachSourceBlockIds) {
            AttachSourceBlockId(dataSets.front(), 0);
        }
        return dataSets.front();
    }

    vtkNew<vtkAppendFilter> append;
    if (options.preserveCompositeBlockBoundaries) {
        append->MergePointsOff();
    } else {
        append->MergePointsOn();
    }
    for (std::size_t i = 0; i < dataSets.size(); ++i) {
        if (options.attachSourceBlockIds) {
            AttachSourceBlockId(dataSets[i], static_cast<int>(i));
        }
        append->AddInputData(dataSets[i]);
    }
    append->Update();
    vtkSmartPointer<vtkDataSet> output = CloneDataSet(append->GetOutput());
    if (output == nullptr || output->GetNumberOfPoints() <= 0) {
        throw std::runtime_error("Failed to assemble the volume blocks from: " +
                                 sourcePath.string());
    }
    return output;
}

vtkSmartPointer<vtkDataSet> ReadOrdinaryVtk(const std::filesystem::path& path,
                                            const ReadOptions& options) {
    const std::string extension = LowerExtension(path);
    vtkSmartPointer<vtkDataObject> object;
    if (extension == ".vtkhdf") {
        object = streamcenterplus::vtkhdf::ReadDataObject(path);
    } else if (extension == ".vtk") {
        vtkNew<vtkGenericDataObjectReader> reader;
        reader->SetFileName(path.string().c_str());
        if (options.readAllArrays) {
            reader->ReadAllScalarsOn();
            reader->ReadAllVectorsOn();
            reader->ReadAllNormalsOn();
            reader->ReadAllTensorsOn();
            reader->ReadAllFieldsOn();
        }
        reader->Update();
        object = reader->GetOutputDataObject(0);
    } else {
        vtkNew<vtkXMLGenericDataObjectReader> reader;
        reader->SetFileName(path.string().c_str());
        reader->UpdateInformation();
        if (!options.readAllArrays) {
            if (vtkDataArraySelection* selection = reader->GetPointDataArraySelection()) {
                selection->DisableAllArrays();
            }
            if (vtkDataArraySelection* selection = reader->GetCellDataArraySelection()) {
                selection->DisableAllArrays();
            }
            if (vtkDataArraySelection* selection = reader->GetColumnArraySelection()) {
                selection->DisableAllArrays();
            }
        }
        reader->Update();
        object = reader->GetOutputDataObject(0);
    }
    return CollapseDataObject(object, options, path);
}

std::vector<double> PipelineTimeSteps(vtkInformation* information) {
    std::vector<double> values;
    if (information == nullptr) {
        return values;
    }
    auto* key = vtkStreamingDemandDrivenPipeline::TIME_STEPS();
    if (!information->Has(key)) {
        return values;
    }
    const int count = information->Length(key);
    values.reserve(static_cast<std::size_t>(std::max(0, count)));
    for (int i = 0; i < count; ++i) {
        const double value = information->Get(key, i);
        if (IsFiniteTime(value)) {
            values.push_back(value);
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

std::optional<std::filesystem::path> ResolveOpenFoamCase(
    const std::filesystem::path& sourcePath) {
    if (IsRegularFileWithExtension(sourcePath, ".foam")) {
        return std::filesystem::absolute(sourcePath);
    }
    if (!std::filesystem::is_directory(sourcePath)) {
        return std::nullopt;
    }
    std::vector<std::filesystem::path> markers;
    for (const auto& entry : std::filesystem::directory_iterator(sourcePath)) {
        if (entry.is_regular_file() && LowerExtension(entry.path()) == ".foam") {
            markers.push_back(entry.path());
        }
    }
    std::sort(markers.begin(), markers.end());
    if (!markers.empty()) {
        return std::filesystem::absolute(markers.front());
    }

    const bool recognizable =
        std::filesystem::is_regular_file(sourcePath / "system" / "controlDict") ||
        std::filesystem::is_directory(sourcePath / "constant" / "polyMesh") ||
        std::filesystem::is_directory(sourcePath / "processor0");
    if (!recognizable) {
        return std::nullopt;
    }
    // VTK treats the .foam path as a marker and reads the case from its parent;
    // no marker is created in the user's case directory.
    return std::filesystem::absolute(sourcePath /
                                     (sourcePath.filename().string() + ".foam"));
}

bool IsDecomposedOpenFoamCase(const std::filesystem::path& marker) {
    const std::filesystem::path root = marker.parent_path();
    return std::filesystem::is_directory(root / "processor0");
}

vtkSmartPointer<vtkPOpenFOAMReader> CreateOpenFoamReader(
    const std::filesystem::path& marker,
    bool readAllArrays) {
    vtkSmartPointer<vtkPOpenFOAMReader> reader =
        vtkSmartPointer<vtkPOpenFOAMReader>::New();
    reader->SetCaseType(IsDecomposedOpenFoamCase(marker)
                            ? vtkPOpenFOAMReader::DECOMPOSED_CASE
                            : vtkPOpenFOAMReader::RECONSTRUCTED_CASE);
    reader->SetFileName(marker.string().c_str());
    reader->SetSkipZeroTime(false);
    reader->SetListTimeStepsByControlDict(false);
    reader->SetCreateCellToPoint(false);
    reader->SetCacheMesh(true);
    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->DisableAllPatchArrays();
    if (readAllArrays) {
        reader->EnableAllCellArrays();
        reader->EnableAllPointArrays();
    } else {
        reader->DisableAllCellArrays();
        reader->DisableAllPointArrays();
    }
    reader->UpdateInformation();
    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->DisableAllPatchArrays();
    for (int i = 0; i < reader->GetNumberOfPatchArrays(); ++i) {
        const char* rawName = reader->GetPatchArrayName(i);
        const std::string name = rawName == nullptr ? std::string() : ToLower(rawName);
        if (name == "internalmesh" || name.find("internal") != std::string::npos) {
            reader->SetPatchArrayStatus(rawName, 1);
        }
    }
    return reader;
}

std::vector<double> OpenFoamTimeSteps(vtkOpenFOAMReader* reader) {
    std::vector<double> values;
    if (reader == nullptr) {
        return values;
    }
    if (vtkDoubleArray* times = reader->GetTimeValues()) {
        values.reserve(static_cast<std::size_t>(times->GetNumberOfTuples()));
        for (vtkIdType i = 0; i < times->GetNumberOfTuples(); ++i) {
            const double value = times->GetValue(i);
            if (IsFiniteTime(value)) {
                values.push_back(value);
            }
        }
    }
    if (values.empty()) {
        values = PipelineTimeSteps(reader->GetOutputInformation(0));
    }
    if (values.empty()) {
        values.push_back(0.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

vtkSmartPointer<vtkDataSet> ReadOpenFoam(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    vtkSmartPointer<vtkPOpenFOAMReader> reader =
        CreateOpenFoamReader(frame.path, options.readAllArrays);
    const double time = frame.readerTime.value_or(frame.time);
    reader->SetTimeValue(time);
    if (reader->UpdateTimeStep(time) == 0) {
        throw std::runtime_error("VTK failed to read OpenFOAM time " +
                                 std::to_string(time) + " from " + frame.path.string());
    }
    return CollapseDataObject(reader->GetOutputDataObject(0), options, frame.path);
}

vtkSmartPointer<vtkGenericEnSightReader> CreateEnSightGoldReader(
    const std::filesystem::path& casePath,
    bool readAllArrays) {
    vtkSmartPointer<vtkGenericEnSightReader> reader =
        vtkSmartPointer<vtkGenericEnSightReader>::New();
    reader->SetCaseFileName(casePath.string().c_str());
    if (readAllArrays) {
        reader->ReadAllVariablesOn();
    } else {
        reader->ReadAllVariablesOff();
    }
    reader->UpdateInformation();
    const int version = reader->GetEnSightVersion();
    if (version != vtkGenericEnSightReader::ENSIGHT_GOLD &&
        version != vtkGenericEnSightReader::ENSIGHT_GOLD_BINARY) {
        throw std::runtime_error("The .case source is not EnSight Gold: " +
                                 casePath.string());
    }
    if (readAllArrays) {
        reader->ReadAllVariablesOn();
        for (int i = 0; i < reader->GetNumberOfPointArrays(); ++i) {
            if (const char* name = reader->GetPointArrayName(i)) {
                reader->SetPointArrayStatus(name, 1);
            }
        }
        for (int i = 0; i < reader->GetNumberOfCellArrays(); ++i) {
            if (const char* name = reader->GetCellArrayName(i)) {
                reader->SetCellArrayStatus(name, 1);
            }
        }
    }
    return reader;
}

std::vector<double> EnSightGoldTimeSteps(vtkGenericEnSightReader* reader) {
    std::vector<double> values =
        reader == nullptr ? std::vector<double>()
                          : PipelineTimeSteps(reader->GetOutputInformation(0));
    if (values.empty() && reader != nullptr) {
        if (vtkDataArrayCollection* sets = reader->GetTimeSets()) {
            for (int setIndex = 0; setIndex < sets->GetNumberOfItems(); ++setIndex) {
                vtkDataArray* array = sets->GetItem(setIndex);
                if (array == nullptr) {
                    continue;
                }
                for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
                    const double value = array->GetComponent(i, 0);
                    if (IsFiniteTime(value)) {
                        values.push_back(value);
                    }
                }
            }
        }
    }
    if (values.empty()) {
        values.push_back(0.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

vtkSmartPointer<vtkDataSet> ReadEnSightGold(const DataSetFrame& frame,
                                           const ReadOptions& options) {
    vtkSmartPointer<vtkGenericEnSightReader> reader =
        CreateEnSightGoldReader(frame.path, options.readAllArrays);
    const double time = frame.readerTime.value_or(frame.time);
    reader->SetTimeValue(time);
    if (reader->UpdateTimeStep(time) == 0) {
        throw std::runtime_error("VTK failed to read EnSight Gold time " +
                                 std::to_string(time) + " from " + frame.path.string());
    }
    return CollapseDataObject(reader->GetOutputDataObject(0), options, frame.path);
}

bool HostIsLittleEndian() {
    const std::uint16_t value = 1;
    return *reinterpret_cast<const unsigned char*>(&value) == 1;
}

template <typename T>
T ByteSwapValue(T value) {
    std::array<unsigned char, sizeof(T)> bytes{};
    std::memcpy(bytes.data(), &value, sizeof(T));
    std::reverse(bytes.begin(), bytes.end());
    std::memcpy(&value, bytes.data(), sizeof(T));
    return value;
}

struct NekDescriptor {
    std::filesystem::path path;
    std::filesystem::path caseDirectory;
    std::string fileTemplate;
    int firstTimeStep = 0;
    int numberOfTimeSteps = 0;
    int integerConversions = 0;
};

struct NekFieldMeta {
    std::filesystem::path path;
    int wordSize = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int elementCount = 0;
    int globalElementCount = 0;
    int fileIndex = 0;
    int fileCount = 1;
    int nxyz = 0;
    int dimensions = 0;
    double time = 0.0;
    std::string code;
    bool littleEndian = true;
    std::uint64_t dataBaseOffset = 0;
    std::vector<std::string> labels;

    std::size_t pointCount() const {
        return static_cast<std::size_t>(elementCount) *
               static_cast<std::size_t>(nxyz);
    }
};

int ValidateNekPrintfTemplate(const std::string& format) {
    int conversions = 0;
    for (std::size_t i = 0; i < format.size(); ++i) {
        if (format[i] != '%') {
            continue;
        }
        if (i + 1 < format.size() && format[i + 1] == '%') {
            ++i;
            continue;
        }
        ++i;
        while (i < format.size() && std::string("-+ #0").find(format[i]) != std::string::npos) {
            ++i;
        }
        while (i < format.size() && std::isdigit(static_cast<unsigned char>(format[i]))) {
            ++i;
        }
        if (i < format.size() && format[i] == '.') {
            ++i;
            while (i < format.size() && std::isdigit(static_cast<unsigned char>(format[i]))) {
                ++i;
            }
        }
        if (i < format.size() &&
            std::string("hljztL").find(format[i]) != std::string::npos) {
            throw std::runtime_error(
                "Nek5000 filetemplate integer conversions must use the default int length.");
        }
        if (i >= format.size() || std::string("diuoxX").find(format[i]) == std::string::npos) {
            throw std::runtime_error(
                "Nek5000 filetemplate contains a non-integer or unsafe printf conversion.");
        }
        ++conversions;
    }
    if (conversions < 1 || conversions > 2) {
        throw std::runtime_error(
            "Nek5000 filetemplate must contain one timestep conversion or file/timestep conversions.");
    }
    return conversions;
}

NekDescriptor ParseNekDescriptor(const std::filesystem::path& descriptorPath) {
    std::ifstream input(descriptorPath);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 descriptor: " +
                                 descriptorPath.string());
    }
    NekDescriptor descriptor;
    descriptor.path = std::filesystem::absolute(descriptorPath);
    descriptor.caseDirectory = descriptor.path.parent_path();
    bool hasTemplate = false;
    bool hasFirst = false;
    bool hasCount = false;
    std::string line;
    while (std::getline(input, line)) {
        line = StripComment(line);
        if (line.empty()) {
            continue;
        }
        const std::size_t colon = line.find(':');
        if (colon == std::string::npos) {
            continue;
        }
        const std::string key = ToLower(TrimAscii(line.substr(0, colon)));
        const std::string value = TrimAscii(line.substr(colon + 1));
        if (key == "filetemplate") {
            descriptor.fileTemplate = value;
            hasTemplate = !value.empty();
        } else if (key == "firsttimestep") {
            descriptor.firstTimeStep = std::stoi(value);
            hasFirst = true;
        } else if (key == "numtimesteps") {
            descriptor.numberOfTimeSteps = std::stoi(value);
            hasCount = true;
        } else if (key == "nek5000" || key == "endian" || key == "version") {
            continue;
        } else {
            throw std::runtime_error("Unsupported Nek5000 descriptor tag: " + key);
        }
    }
    if (!hasTemplate || !hasFirst || !hasCount || descriptor.numberOfTimeSteps <= 0) {
        throw std::runtime_error(
            "Nek5000 descriptor must define filetemplate, firsttimestep, and numtimesteps.");
    }
    descriptor.integerConversions = ValidateNekPrintfTemplate(descriptor.fileTemplate);
    return descriptor;
}

std::filesystem::path FormatNekFieldPath(const NekDescriptor& descriptor,
                                         int relativeFrameIndex,
                                         int fileIndex) {
    const int timeStep = descriptor.firstTimeStep + relativeFrameIndex;
    std::array<char, 4096> buffer{};
    const int written = descriptor.integerConversions == 2
                            ? std::snprintf(buffer.data(), buffer.size(),
                                            descriptor.fileTemplate.c_str(), fileIndex, timeStep)
                            : std::snprintf(buffer.data(), buffer.size(),
                                            descriptor.fileTemplate.c_str(), timeStep);
    if (written < 0 || static_cast<std::size_t>(written) >= buffer.size()) {
        throw std::runtime_error("Nek5000 filetemplate produced an invalid path.");
    }
    std::filesystem::path path(buffer.data());
    if (path.is_relative()) {
        path = descriptor.caseDirectory / path;
    }
    return path;
}

float DecodeEndianTag(const std::array<char, 4>& bytes, bool littleEndian) {
    std::array<char, 4> ordered = bytes;
    if (littleEndian != HostIsLittleEndian()) {
        std::reverse(ordered.begin(), ordered.end());
    }
    float value = 0.0f;
    std::memcpy(&value, ordered.data(), sizeof(value));
    return value;
}

std::vector<std::string> SplitWhitespace(const std::string& value) {
    std::istringstream stream(value);
    std::vector<std::string> parts;
    std::string token;
    while (stream >> token) {
        parts.push_back(token);
    }
    return parts;
}

std::vector<std::string> NekLabelsForCode(std::string code, int dimensions) {
    std::transform(code.begin(), code.end(), code.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    std::vector<std::string> labels;
    if (code.find('X') != std::string::npos) {
        labels.push_back("x");
        labels.push_back("y");
        if (dimensions == 3) {
            labels.push_back("z");
        }
    }
    if (code.find('U') != std::string::npos) {
        labels.push_back("u");
        labels.push_back("v");
        if (dimensions == 3) {
            labels.push_back("w");
        }
    }
    if (code.find('P') != std::string::npos) {
        labels.push_back("p");
    }
    if (code.find('T') != std::string::npos) {
        labels.push_back("t");
    }
    const std::size_t scalar = code.find('S');
    if (scalar != std::string::npos) {
        std::size_t end = scalar + 1;
        while (end < code.size() && std::isdigit(static_cast<unsigned char>(code[end]))) {
            ++end;
        }
        const std::string countToken = code.substr(scalar + 1, end - scalar - 1);
        const int count = countToken.empty() ? 1 : std::stoi(countToken);
        for (int i = 0; i < count; ++i) {
            std::ostringstream name;
            name << 's' << std::setw(2) << std::setfill('0') << (i + 1);
            labels.push_back(name.str());
        }
    }
    return labels;
}

NekFieldMeta ReadNekHeader(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 field file: " + path.string());
    }
    std::array<char, 132> headerBytes{};
    std::array<char, 4> tagBytes{};
    input.read(headerBytes.data(), static_cast<std::streamsize>(headerBytes.size()));
    input.read(tagBytes.data(), static_cast<std::streamsize>(tagBytes.size()));
    if (!input) {
        throw std::runtime_error("Cannot read Nek5000 field header: " + path.string());
    }
    const std::vector<std::string> parts =
        SplitWhitespace(std::string(headerBytes.data(), headerBytes.size()));
    if (parts.size() < 12 || parts[0] != "#std") {
        throw std::runtime_error("Not a supported Nek5000 #std field file: " + path.string());
    }

    NekFieldMeta meta;
    meta.path = path;
    meta.wordSize = std::stoi(parts[1]);
    meta.nx = std::stoi(parts[2]);
    meta.ny = std::stoi(parts[3]);
    meta.nz = std::stoi(parts[4]);
    meta.elementCount = std::stoi(parts[5]);
    meta.globalElementCount = std::stoi(parts[6]);
    std::string timeToken = parts[7];
    std::replace(timeToken.begin(), timeToken.end(), 'D', 'E');
    std::replace(timeToken.begin(), timeToken.end(), 'd', 'E');
    meta.time = std::stod(timeToken);
    meta.fileIndex = std::stoi(parts[9]);
    meta.fileCount = std::max(1, std::stoi(parts[10]));
    meta.code = parts[11];
    meta.dimensions = meta.nz > 1 ? 3 : 2;
    meta.nxyz = meta.nx * meta.ny * meta.nz;
    const float littleTag = DecodeEndianTag(tagBytes, true);
    const float bigTag = DecodeEndianTag(tagBytes, false);
    if (std::abs(littleTag - 6.54321f) < 1.0e-4f) {
        meta.littleEndian = true;
    } else if (std::abs(bigTag - 6.54321f) < 1.0e-4f) {
        meta.littleEndian = false;
    } else {
        throw std::runtime_error("Cannot determine Nek5000 field endianness: " +
                                 path.string());
    }
    if (meta.wordSize != 4 && meta.wordSize != 8) {
        throw std::runtime_error("Nek5000 field word size must be 4 or 8 bytes: " +
                                 path.string());
    }
    if (meta.nx < 2 || meta.ny < 2 || meta.nz < 1 || meta.elementCount <= 0 ||
        meta.globalElementCount <= 0 || !IsFiniteTime(meta.time)) {
        throw std::runtime_error("Invalid Nek5000 field dimensions or metadata: " +
                                 path.string());
    }
    meta.dataBaseOffset = 136ull +
                          4ull * static_cast<std::uint64_t>(meta.elementCount);
    meta.labels = NekLabelsForCode(meta.code, meta.dimensions);
    return meta;
}

void ValidateNekPieceCompatible(const NekFieldMeta& reference,
                                const NekFieldMeta& current,
                                bool requireCode) {
    if ((requireCode && reference.wordSize != current.wordSize) ||
        reference.nx != current.nx ||
        reference.ny != current.ny || reference.nz != current.nz ||
        reference.globalElementCount != current.globalElementCount ||
        (requireCode && reference.fileCount != current.fileCount)) {
        throw std::runtime_error("Nek5000 field pieces have incompatible dimensions: " +
                                 current.path.string());
    }
    if (requireCode && reference.code != current.code) {
        throw std::runtime_error("Nek5000 field pieces have different rdcode values: " +
                                 current.path.string());
    }
}

std::vector<NekFieldMeta> ReadNekPieceSet(const NekDescriptor& descriptor,
                                          int relativeFrameIndex) {
    NekFieldMeta first =
        ReadNekHeader(FormatNekFieldPath(descriptor, relativeFrameIndex, 0));
    if (first.fileCount > 1 && descriptor.integerConversions != 2) {
        throw std::runtime_error(
            "Multi-file Nek5000 output requires filetemplate to contain file and timestep conversions.");
    }
    std::vector<NekFieldMeta> pieces;
    pieces.reserve(static_cast<std::size_t>(first.fileCount));
    pieces.push_back(first);
    for (int fileIndex = 1; fileIndex < first.fileCount; ++fileIndex) {
        NekFieldMeta piece =
            ReadNekHeader(FormatNekFieldPath(descriptor, relativeFrameIndex, fileIndex));
        ValidateNekPieceCompatible(first, piece, true);
        pieces.push_back(std::move(piece));
    }
    std::size_t totalElements = 0;
    for (const NekFieldMeta& piece : pieces) {
        totalElements += static_cast<std::size_t>(piece.elementCount);
    }
    if (totalElements != static_cast<std::size_t>(first.globalElementCount)) {
        throw std::runtime_error("Nek5000 field pieces do not cover the advertised global elements.");
    }
    return pieces;
}

std::vector<int> ReadNekElementMap(const NekFieldMeta& meta) {
    std::ifstream input(meta.path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 element map: " + meta.path.string());
    }
    input.seekg(136, std::ios::beg);
    std::vector<std::int32_t> raw(static_cast<std::size_t>(meta.elementCount));
    input.read(reinterpret_cast<char*>(raw.data()),
               static_cast<std::streamsize>(raw.size() * sizeof(std::int32_t)));
    if (!input) {
        throw std::runtime_error("Cannot read Nek5000 element map: " + meta.path.string());
    }
    std::vector<int> ids(raw.size(), 0);
    for (std::size_t i = 0; i < raw.size(); ++i) {
        std::int32_t value = raw[i];
        if (meta.littleEndian != HostIsLittleEndian()) {
            value = ByteSwapValue(value);
        }
        ids[i] = static_cast<int>(value);
    }
    return ids;
}

void ReadNekRawValues(std::ifstream& input,
                      const NekFieldMeta& meta,
                      std::uint64_t byteOffset,
                      std::size_t count,
                      double* destination) {
    input.seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);
    if (!input) {
        throw std::runtime_error("Cannot seek Nek5000 field data: " + meta.path.string());
    }
    constexpr std::size_t kChunk = 1u << 20;
    std::size_t written = 0;
    if (meta.wordSize == 4) {
        std::vector<float> buffer;
        while (written < count) {
            const std::size_t current = std::min(kChunk, count - written);
            buffer.resize(current);
            input.read(reinterpret_cast<char*>(buffer.data()),
                       static_cast<std::streamsize>(current * sizeof(float)));
            if (!input) {
                throw std::runtime_error("Cannot read Nek5000 float field data: " +
                                         meta.path.string());
            }
            for (std::size_t i = 0; i < current; ++i) {
                float value = buffer[i];
                if (meta.littleEndian != HostIsLittleEndian()) {
                    value = ByteSwapValue(value);
                }
                destination[written + i] = static_cast<double>(value);
            }
            written += current;
        }
        return;
    }
    std::vector<double> buffer;
    while (written < count) {
        const std::size_t current = std::min(kChunk, count - written);
        buffer.resize(current);
        input.read(reinterpret_cast<char*>(buffer.data()),
                   static_cast<std::streamsize>(current * sizeof(double)));
        if (!input) {
            throw std::runtime_error("Cannot read Nek5000 double field data: " +
                                     meta.path.string());
        }
        for (std::size_t i = 0; i < current; ++i) {
            double value = buffer[i];
            if (meta.littleEndian != HostIsLittleEndian()) {
                value = ByteSwapValue(value);
            }
            destination[written + i] = value;
        }
        written += current;
    }
}

struct NekVectorLayout {
    std::uint64_t blockOffset = 0;
    int component = 0;
};

std::optional<NekVectorLayout> ResolveNekVectorLayout(const NekFieldMeta& meta,
                                                       const std::string& label) {
    const std::string lower = ToLower(label);
    const int dimensions = meta.dimensions;
    if (lower == "x" || lower == "y" || lower == "z") {
        if (meta.code.find('X') == std::string::npos ||
            (lower == "z" && dimensions != 3)) {
            return std::nullopt;
        }
        return NekVectorLayout{meta.dataBaseOffset,
                               lower == "x" ? 0 : (lower == "y" ? 1 : 2)};
    }
    if (lower == "u" || lower == "v" || lower == "w") {
        if (meta.code.find('U') == std::string::npos ||
            (lower == "w" && dimensions != 3)) {
            return std::nullopt;
        }
        const std::uint64_t coordinateComponents =
            meta.code.find('X') == std::string::npos
                ? 0ull
                : static_cast<std::uint64_t>(dimensions);
        const std::uint64_t blockOffset =
            meta.dataBaseOffset + coordinateComponents *
                                      static_cast<std::uint64_t>(meta.elementCount) *
                                      static_cast<std::uint64_t>(meta.nxyz) *
                                      static_cast<std::uint64_t>(meta.wordSize);
        return NekVectorLayout{blockOffset,
                               lower == "u" ? 0 : (lower == "v" ? 1 : 2)};
    }
    return std::nullopt;
}

int NekLabelIndex(const NekFieldMeta& meta, const std::string& label) {
    const std::string lower = ToLower(label);
    for (std::size_t i = 0; i < meta.labels.size(); ++i) {
        if (meta.labels[i] == lower) {
            return static_cast<int>(i);
        }
    }
    throw std::runtime_error("Nek5000 field has no component '" + label +
                             "': " + meta.path.string());
}

void ReadNekComponent(const NekFieldMeta& meta,
                      const std::string& label,
                      std::vector<double>* values) {
    if (values == nullptr) {
        throw std::runtime_error("Internal error: null Nek5000 component destination.");
    }
    values->assign(meta.pointCount(), 0.0);
    std::ifstream input(meta.path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 field data: " + meta.path.string());
    }
    if (const std::optional<NekVectorLayout> layout =
            ResolveNekVectorLayout(meta, label)) {
        const std::size_t nxyz = static_cast<std::size_t>(meta.nxyz);
        const std::size_t valuesPerElement =
            static_cast<std::size_t>(meta.dimensions) * nxyz;
        std::vector<double> raw(valuesPerElement);
        for (int element = 0; element < meta.elementCount; ++element) {
            const std::uint64_t offset =
                layout->blockOffset +
                static_cast<std::uint64_t>(element) *
                    static_cast<std::uint64_t>(valuesPerElement) *
                    static_cast<std::uint64_t>(meta.wordSize);
            ReadNekRawValues(input, meta, offset, valuesPerElement, raw.data());
            const std::size_t source =
                static_cast<std::size_t>(layout->component) * nxyz;
            const std::size_t target = static_cast<std::size_t>(element) * nxyz;
            std::copy_n(raw.data() + source, nxyz, values->data() + target);
        }
        return;
    }
    const int componentIndex = NekLabelIndex(meta, label);
    const std::uint64_t offset =
        meta.dataBaseOffset + static_cast<std::uint64_t>(componentIndex) *
                                  static_cast<std::uint64_t>(meta.elementCount) *
                                  static_cast<std::uint64_t>(meta.nxyz) *
                                  static_cast<std::uint64_t>(meta.wordSize);
    ReadNekRawValues(input, meta, offset, meta.pointCount(), values->data());
}

struct NekPieceLookup {
    std::vector<std::vector<int>> ids;
    std::vector<std::pair<std::size_t, std::size_t>> outputToPieceElement;
    std::vector<int> outputElementIds;
};

NekPieceLookup BuildNekPieceLookup(const std::vector<NekFieldMeta>& dataPieces) {
    NekPieceLookup lookup;
    struct PieceElement {
        int id = 0;
        std::size_t piece = 0;
        std::size_t local = 0;
    };
    std::vector<PieceElement> elements;
    for (std::size_t piece = 0; piece < dataPieces.size(); ++piece) {
        lookup.ids.push_back(ReadNekElementMap(dataPieces[piece]));
        for (std::size_t local = 0; local < lookup.ids.back().size(); ++local) {
            elements.push_back({lookup.ids.back()[local], piece, local});
        }
    }
    std::sort(elements.begin(), elements.end(), [](const PieceElement& lhs,
                                                    const PieceElement& rhs) {
        return lhs.id < rhs.id;
    });
    lookup.outputToPieceElement.reserve(elements.size());
    lookup.outputElementIds.reserve(elements.size());
    for (std::size_t output = 0; output < elements.size(); ++output) {
        if (output > 0 && elements[output - 1].id == elements[output].id) {
            throw std::runtime_error("Duplicate Nek5000 global element id " +
                                     std::to_string(elements[output].id) + ".");
        }
        lookup.outputToPieceElement.push_back(
            {elements[output].piece, elements[output].local});
        lookup.outputElementIds.push_back(elements[output].id);
    }
    return lookup;
}

std::unordered_map<int, std::pair<std::size_t, std::size_t>> BuildNekElementMap(
    const std::vector<NekFieldMeta>& pieces) {
    std::unordered_map<int, std::pair<std::size_t, std::size_t>> result;
    for (std::size_t piece = 0; piece < pieces.size(); ++piece) {
        const std::vector<int> ids = ReadNekElementMap(pieces[piece]);
        for (std::size_t local = 0; local < ids.size(); ++local) {
            if (!result.emplace(ids[local], std::make_pair(piece, local)).second) {
                throw std::runtime_error("Duplicate Nek5000 coordinate element id " +
                                         std::to_string(ids[local]) + ".");
            }
        }
    }
    return result;
}

void ScatterNekComponent(const std::vector<NekFieldMeta>& pieces,
                         const NekPieceLookup& outputLookup,
                         const std::string& label,
                         std::vector<double>* output) {
    if (output == nullptr || pieces.empty()) {
        throw std::runtime_error("Internal error while assembling Nek5000 component.");
    }
    const std::size_t nxyz = static_cast<std::size_t>(pieces.front().nxyz);
    output->assign(outputLookup.outputToPieceElement.size() * nxyz, 0.0);
    std::vector<std::vector<double>> pieceValues(pieces.size());
    for (std::size_t piece = 0; piece < pieces.size(); ++piece) {
        ReadNekComponent(pieces[piece], label, &pieceValues[piece]);
    }
    for (std::size_t outputElement = 0;
         outputElement < outputLookup.outputToPieceElement.size(); ++outputElement) {
        const auto [piece, local] = outputLookup.outputToPieceElement[outputElement];
        std::copy_n(pieceValues[piece].data() + local * nxyz, nxyz,
                    output->data() + outputElement * nxyz);
    }
}

void ScatterNekCoordinateComponent(const std::vector<NekFieldMeta>& coordinatePieces,
                                   const NekPieceLookup& outputLookup,
                                   const std::string& label,
                                   std::vector<double>* output) {
    if (output == nullptr || coordinatePieces.empty()) {
        throw std::runtime_error("Internal error while assembling Nek5000 coordinates.");
    }
    const std::size_t nxyz = static_cast<std::size_t>(coordinatePieces.front().nxyz);
    const auto coordinateMap = BuildNekElementMap(coordinatePieces);
    std::vector<std::vector<double>> pieceValues(coordinatePieces.size());
    for (std::size_t piece = 0; piece < coordinatePieces.size(); ++piece) {
        ReadNekComponent(coordinatePieces[piece], label, &pieceValues[piece]);
    }
    output->assign(outputLookup.outputElementIds.size() * nxyz, 0.0);
    for (std::size_t outputElement = 0;
         outputElement < outputLookup.outputElementIds.size(); ++outputElement) {
        const int id = outputLookup.outputElementIds[outputElement];
        const auto found = coordinateMap.find(id);
        if (found == coordinateMap.end()) {
            throw std::runtime_error("Nek5000 coordinate field is missing global element id " +
                                     std::to_string(id) + ".");
        }
        const auto [piece, local] = found->second;
        std::copy_n(pieceValues[piece].data() + local * nxyz, nxyz,
                    output->data() + outputElement * nxyz);
    }
}

double LegendreP(int order, double x) {
    if (order == 0) {
        return 1.0;
    }
    if (order == 1) {
        return x;
    }
    double previous = 1.0;
    double current = x;
    for (int k = 2; k <= order; ++k) {
        const double next =
            ((2.0 * static_cast<double>(k) - 1.0) * x * current -
             (static_cast<double>(k) - 1.0) * previous) /
            static_cast<double>(k);
        previous = current;
        current = next;
    }
    return current;
}

double LegendreDerivative(int order, double x) {
    if (order == 0) {
        return 0.0;
    }
    return static_cast<double>(order) *
           (x * LegendreP(order, x) - LegendreP(order - 1, x)) /
           (x * x - 1.0);
}

std::vector<double> GllNodes(int count) {
    if (count < 2) {
        throw std::runtime_error("GLL node count must be at least two.");
    }
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const int order = count - 1;
    std::vector<double> nodes(static_cast<std::size_t>(count), 0.0);
    nodes.front() = -1.0;
    nodes.back() = 1.0;
    for (int i = 1; i < order; ++i) {
        double x = -std::cos(kPi * static_cast<double>(i) /
                             static_cast<double>(order));
        for (int iteration = 0; iteration < 60; ++iteration) {
            const double p = LegendreP(order, x);
            const double dp = LegendreDerivative(order, x);
            const double ddp =
                (2.0 * x * dp - static_cast<double>(order) *
                                      static_cast<double>(order + 1) * p) /
                (1.0 - x * x);
            const double delta = dp / ddp;
            x -= delta;
            if (std::abs(delta) < 2.0e-15) {
                break;
            }
        }
        nodes[static_cast<std::size_t>(i)] = x;
    }
    return nodes;
}

struct OutwardInterval {
    double lower = 0.0;
    double upper = 0.0;
};

double RoundDown(double value) {
    return std::nextafter(value, -std::numeric_limits<double>::infinity());
}

double RoundUp(double value) {
    return std::nextafter(value, std::numeric_limits<double>::infinity());
}

OutwardInterval IntervalPoint(double value) {
    return {value, value};
}

OutwardInterval IntervalAdd(const OutwardInterval& lhs,
                            const OutwardInterval& rhs) {
    return {RoundDown(lhs.lower + rhs.lower),
            RoundUp(lhs.upper + rhs.upper)};
}

OutwardInterval IntervalNegate(const OutwardInterval& value) {
    return {-value.upper, -value.lower};
}

OutwardInterval IntervalSubtract(const OutwardInterval& lhs,
                                 const OutwardInterval& rhs) {
    return IntervalAdd(lhs, IntervalNegate(rhs));
}

OutwardInterval IntervalMultiply(const OutwardInterval& lhs,
                                 const OutwardInterval& rhs) {
    const std::array<double, 4> products = {
        lhs.lower * rhs.lower, lhs.lower * rhs.upper,
        lhs.upper * rhs.lower, lhs.upper * rhs.upper};
    const auto minimum = std::min_element(products.begin(), products.end());
    const auto maximum = std::max_element(products.begin(), products.end());
    return {RoundDown(*minimum), RoundUp(*maximum)};
}

OutwardInterval IntervalDivide(const OutwardInterval& numerator,
                               const OutwardInterval& denominator) {
    if (denominator.lower <= 0.0 && denominator.upper >= 0.0) {
        throw std::runtime_error(
            "Cannot form a conservative Bernstein bound through a zero interval.");
    }
    const std::array<double, 4> quotients = {
        numerator.lower / denominator.lower,
        numerator.lower / denominator.upper,
        numerator.upper / denominator.lower,
        numerator.upper / denominator.upper};
    const auto minimum = std::min_element(quotients.begin(), quotients.end());
    const auto maximum = std::max_element(quotients.begin(), quotients.end());
    return {RoundDown(*minimum), RoundUp(*maximum)};
}

std::vector<OutwardInterval> MultiplyBernsteinByLinear(
    const std::vector<OutwardInterval>& polynomial,
    const std::array<OutwardInterval, 2>& linear) {
    if (polynomial.empty()) {
        throw std::runtime_error("Cannot multiply an empty Bernstein polynomial.");
    }
    const std::size_t degree = polynomial.size() - 1;
    std::vector<OutwardInterval> product(degree + 2);
    for (std::size_t k = 0; k <= degree + 1; ++k) {
        bool initialized = false;
        OutwardInterval coefficient{};
        if (k <= degree) {
            const OutwardInterval ratio = IntervalDivide(
                IntervalPoint(static_cast<double>(degree + 1 - k)),
                IntervalPoint(static_cast<double>(degree + 1)));
            coefficient = IntervalMultiply(
                IntervalMultiply(polynomial[k], linear[0]), ratio);
            initialized = true;
        }
        if (k > 0) {
            const OutwardInterval ratio = IntervalDivide(
                IntervalPoint(static_cast<double>(k)),
                IntervalPoint(static_cast<double>(degree + 1)));
            const OutwardInterval contribution = IntervalMultiply(
                IntervalMultiply(polynomial[k - 1], linear[1]), ratio);
            coefficient = initialized ? IntervalAdd(coefficient, contribution)
                                      : contribution;
        }
        product[k] = coefficient;
    }
    return product;
}

// On every open interval between adjacent interpolation nodes, the signs of
// all Lagrange basis functions are fixed. Therefore the Lebesgue function is
// one ordinary polynomial there. We construct each signed basis directly in
// the Bernstein basis on two half-intervals. Outward-rounded interval
// arithmetic encloses every control coefficient, and the Bernstein convex-hull
// property then supplies a rigorous upper bound without sampling the function.
double ConservativeLebesgueBound(const std::vector<double>& interpolationNodes) {
    const std::size_t count = interpolationNodes.size();
    if (count <= 2) {
        return 1.0;
    }
    double upperBound = 1.0;
    for (std::size_t gap = 0; gap + 1 < count; ++gap) {
        const double lowerNode = interpolationNodes[gap];
        const double upperNode = interpolationNodes[gap + 1];
        const double midpoint = lowerNode + 0.5 * (upperNode - lowerNode);
        const std::array<std::array<double, 2>, 2> subIntervals = {{
            {lowerNode, midpoint}, {midpoint, upperNode}}};
        for (const auto& subInterval : subIntervals) {
            std::vector<OutwardInterval> lebesgueCoefficients(
                count, IntervalPoint(0.0));
            for (std::size_t basisIndex = 0; basisIndex < count; ++basisIndex) {
                std::vector<OutwardInterval> basis(1, IntervalPoint(1.0));
                int sign = 1;
                for (std::size_t other = 0; other < count; ++other) {
                    if (other == basisIndex) {
                        continue;
                    }
                    const OutwardInterval denominator = IntervalSubtract(
                        IntervalPoint(interpolationNodes[basisIndex]),
                        IntervalPoint(interpolationNodes[other]));
                    const std::array<OutwardInterval, 2> factor = {
                        IntervalDivide(
                            IntervalSubtract(IntervalPoint(subInterval[0]),
                                             IntervalPoint(interpolationNodes[other])),
                            denominator),
                        IntervalDivide(
                            IntervalSubtract(IntervalPoint(subInterval[1]),
                                             IntervalPoint(interpolationNodes[other])),
                            denominator)};
                    basis = MultiplyBernsteinByLinear(basis, factor);

                    const bool numeratorNegative = other > gap;
                    const bool denominatorNegative = other > basisIndex;
                    if (numeratorNegative != denominatorNegative) {
                        sign = -sign;
                    }
                }
                for (std::size_t coefficient = 0; coefficient < count;
                     ++coefficient) {
                    const OutwardInterval signedCoefficient =
                        sign > 0 ? basis[coefficient]
                                 : IntervalNegate(basis[coefficient]);
                    lebesgueCoefficients[coefficient] = IntervalAdd(
                        lebesgueCoefficients[coefficient], signedCoefficient);
                }
            }
            for (const OutwardInterval& coefficient : lebesgueCoefficients) {
                upperBound = std::max(upperBound, coefficient.upper);
            }
        }
    }
    if (!std::isfinite(upperBound)) {
        throw std::runtime_error(
            "Cannot construct a finite conservative GLL Lebesgue bound.");
    }
    return RoundUp(upperBound);
}

std::vector<double> GllWeights(const std::vector<double>& nodes) {
    const int count = static_cast<int>(nodes.size());
    const int order = count - 1;
    std::vector<double> weights(nodes.size(), 0.0);
    for (int i = 0; i < count; ++i) {
        const double p = LegendreP(order, nodes[static_cast<std::size_t>(i)]);
        weights[static_cast<std::size_t>(i)] =
            2.0 / (static_cast<double>(order) *
                   static_cast<double>(order + 1) * p * p);
    }
    return weights;
}

std::vector<double> DerivativeMatrix(const std::vector<double>& nodes) {
    const int count = static_cast<int>(nodes.size());
    std::vector<double> bary(static_cast<std::size_t>(count), 1.0);
    for (int j = 0; j < count; ++j) {
        for (int k = 0; k < count; ++k) {
            if (j != k) {
                bary[static_cast<std::size_t>(j)] /=
                    nodes[static_cast<std::size_t>(j)] -
                    nodes[static_cast<std::size_t>(k)];
            }
        }
    }
    std::vector<double> matrix(static_cast<std::size_t>(count) *
                                   static_cast<std::size_t>(count),
                               0.0);
    for (int i = 0; i < count; ++i) {
        double sum = 0.0;
        for (int j = 0; j < count; ++j) {
            if (i == j) {
                continue;
            }
            const double value =
                bary[static_cast<std::size_t>(j)] /
                (bary[static_cast<std::size_t>(i)] *
                 (nodes[static_cast<std::size_t>(i)] -
                  nodes[static_cast<std::size_t>(j)]));
            matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(count) +
                   static_cast<std::size_t>(j)] = value;
            sum += value;
        }
        matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(count) +
               static_cast<std::size_t>(i)] = -sum;
    }
    return matrix;
}

std::vector<double> BarycentricWeights(const std::vector<double>& nodes) {
    std::vector<double> weights(nodes.size(), 1.0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            if (i != j) {
                weights[i] /= nodes[i] - nodes[j];
            }
        }
    }
    return weights;
}

void EvaluateNodalBasis(const std::vector<double>& nodes,
                        const std::vector<double>& barycentricWeights,
                        const std::vector<double>& nodalDerivativeMatrix,
                        double coordinate,
                        std::vector<double>* basis,
                        std::vector<double>* derivative) {
    if (basis == nullptr || derivative == nullptr ||
        nodes.size() != barycentricWeights.size() ||
        nodalDerivativeMatrix.size() != nodes.size() * nodes.size()) {
        throw std::runtime_error("Invalid nodal-basis evaluation workspace.");
    }
    basis->assign(nodes.size(), 0.0);
    derivative->assign(nodes.size(), 0.0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, std::abs(nodes[i]));
        if (std::abs(coordinate - nodes[i]) <= tolerance) {
            (*basis)[i] = 1.0;
            std::copy_n(nodalDerivativeMatrix.data() + i * nodes.size(),
                        nodes.size(), derivative->data());
            return;
        }
    }

    std::vector<double> scaled(nodes.size(), 0.0);
    double denominator = 0.0;
    double squaredDenominator = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double distance = coordinate - nodes[i];
        scaled[i] = barycentricWeights[i] / distance;
        denominator += scaled[i];
        squaredDenominator += scaled[i] / distance;
    }
    if (!std::isfinite(denominator) || denominator == 0.0) {
        throw std::runtime_error("Failed to evaluate a GLL barycentric basis.");
    }
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double value = scaled[i] / denominator;
        (*basis)[i] = value;
        (*derivative)[i] =
            value * (squaredDenominator / denominator -
                     1.0 / (coordinate - nodes[i]));
    }
}

double SquaredNorm3(const std::array<double, 3>& value) {
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

std::size_t NekLocalIndex(int nx, int ny, int i, int j, int k) {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nx) +
           static_cast<std::size_t>(i);
}

struct NekGeometryMetrics {
    std::vector<double> referenceWeights;
    std::vector<double> jacobians;
    std::vector<double> physicalWeights;
};

NekGeometryMetrics ComputeNekGeometryMetrics(vtkDataSet* dataSet,
                                             const Nek5000SpectralMetadata& meta,
                                             bool includeReferenceWeights,
                                             bool includeJacobians,
                                             bool includePhysicalWeights) {
    if (dataSet == nullptr ||
        static_cast<std::size_t>(dataSet->GetNumberOfPoints()) != meta.pointCount()) {
        throw std::runtime_error("Nek5000 geometry and GLL metadata point counts differ.");
    }
    const std::vector<double> nodesX = GllNodes(meta.nx);
    const std::vector<double> nodesY = GllNodes(meta.ny);
    const std::vector<double> nodesZ = meta.is3D() ? GllNodes(meta.nz)
                                                   : std::vector<double>{0.0};
    const std::vector<double> weightsX = GllWeights(nodesX);
    const std::vector<double> weightsY = GllWeights(nodesY);
    const std::vector<double> weightsZ = meta.is3D() ? GllWeights(nodesZ)
                                                     : std::vector<double>{1.0};
    const std::vector<double> dx = DerivativeMatrix(nodesX);
    const std::vector<double> dy = DerivativeMatrix(nodesY);
    const std::vector<double> dz = meta.is3D() ? DerivativeMatrix(nodesZ)
                                               : std::vector<double>();

    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(dataSet);
    vtkDataArray* coordinates =
        pointSet == nullptr || pointSet->GetPoints() == nullptr
            ? nullptr
            : pointSet->GetPoints()->GetData();
    if (coordinates == nullptr || coordinates->GetNumberOfComponents() != 3) {
        throw std::runtime_error(
            "Nek5000 geometry must expose a three-component point-coordinate array.");
    }
    vtkDoubleArray* doubleCoordinates = vtkDoubleArray::SafeDownCast(coordinates);
    const double* coordinateValues =
        doubleCoordinates == nullptr ? nullptr : doubleCoordinates->GetPointer(0);
    const auto coordinate = [coordinates, coordinateValues](std::size_t pointIndex,
                                                            int component) {
        return coordinateValues == nullptr
                   ? coordinates->GetComponent(static_cast<vtkIdType>(pointIndex), component)
                   : coordinateValues[3 * pointIndex + static_cast<std::size_t>(component)];
    };

    NekGeometryMetrics result;
    if (includeReferenceWeights) {
        result.referenceWeights.resize(meta.pointCount(), 0.0);
    }
    if (includeJacobians) {
        result.jacobians.resize(meta.pointCount(), 0.0);
    }
    if (includePhysicalWeights) {
        result.physicalWeights.resize(meta.pointCount(), 0.0);
    }
    const std::size_t nxyz = meta.pointsPerElement();
    for (std::size_t element = 0; element < meta.elementCount; ++element) {
        const std::size_t base = element * nxyz;
        for (int k = 0; k < meta.nz; ++k) {
            for (int j = 0; j < meta.ny; ++j) {
                for (int i = 0; i < meta.nx; ++i) {
                    std::array<double, 3> xr{};
                    std::array<double, 3> xs{};
                    std::array<double, 3> xt{};
                    for (int m = 0; m < meta.nx; ++m) {
                        const double d = dx[static_cast<std::size_t>(i) *
                                                static_cast<std::size_t>(meta.nx) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, m, j, k);
                        for (int component = 0; component < 3; ++component) {
                            xr[static_cast<std::size_t>(component)] +=
                                d * coordinate(index, component);
                        }
                    }
                    for (int m = 0; m < meta.ny; ++m) {
                        const double d = dy[static_cast<std::size_t>(j) *
                                                static_cast<std::size_t>(meta.ny) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, i, m, k);
                        for (int component = 0; component < 3; ++component) {
                            xs[static_cast<std::size_t>(component)] +=
                                d * coordinate(index, component);
                        }
                    }
                    if (meta.is3D()) {
                        for (int m = 0; m < meta.nz; ++m) {
                            const double d = dz[static_cast<std::size_t>(k) *
                                                    static_cast<std::size_t>(meta.nz) +
                                                static_cast<std::size_t>(m)];
                            const std::size_t index =
                                base + NekLocalIndex(meta.nx, meta.ny, i, j, m);
                            for (int component = 0; component < 3; ++component) {
                                xt[static_cast<std::size_t>(component)] +=
                                    d * coordinate(index, component);
                            }
                        }
                    }
                    double jacobian = 0.0;
                    if (meta.is3D()) {
                        jacobian =
                            xr[0] * (xs[1] * xt[2] - xt[1] * xs[2]) -
                            xs[0] * (xr[1] * xt[2] - xt[1] * xr[2]) +
                            xt[0] * (xr[1] * xs[2] - xs[1] * xr[2]);
                    } else {
                        const std::array<double, 3> cross = {
                            xr[1] * xs[2] - xr[2] * xs[1],
                            xr[2] * xs[0] - xr[0] * xs[2],
                            xr[0] * xs[1] - xr[1] * xs[0],
                        };
                        jacobian = std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] +
                                             cross[2] * cross[2]);
                    }
                    if (!std::isfinite(jacobian) || std::abs(jacobian) <= 1.0e-30) {
                        throw std::runtime_error(
                            "Degenerate Nek5000 spectral-element Jacobian at element " +
                            std::to_string(element) + ".");
                    }
                    const double reference =
                        weightsX[static_cast<std::size_t>(i)] *
                        weightsY[static_cast<std::size_t>(j)] *
                        weightsZ[static_cast<std::size_t>(k)];
                    const std::size_t index =
                        base + NekLocalIndex(meta.nx, meta.ny, i, j, k);
                    if (includeReferenceWeights) {
                        result.referenceWeights[index] = reference;
                    }
                    if (includeJacobians) {
                        result.jacobians[index] = jacobian;
                    }
                    if (includePhysicalWeights) {
                        result.physicalWeights[index] =
                            reference * std::abs(jacobian);
                    }
                }
            }
        }
    }
    return result;
}

void AddDoublePointArray(vtkDataSet* dataSet,
                         const char* name,
                         const std::vector<double>& values) {
    vtkNew<vtkDoubleArray> array;
    array->SetName(name);
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(static_cast<vtkIdType>(values.size()));
    std::copy(values.begin(), values.end(), array->GetPointer(0));
    dataSet->GetPointData()->AddArray(array);
}

void AttachNekQuadratureArrays(vtkDataSet* dataSet,
                               const NekGeometryMetrics& metrics) {
    if (!metrics.referenceWeights.empty()) {
        AddDoublePointArray(dataSet, kNek5000ReferenceWeightArray,
                            metrics.referenceWeights);
    }
    if (!metrics.jacobians.empty()) {
        AddDoublePointArray(dataSet, kNek5000JacobianDeterminantArray,
                            metrics.jacobians);
    }
    if (!metrics.physicalWeights.empty()) {
        AddDoublePointArray(dataSet, kNek5000QuadratureWeightArray,
                            metrics.physicalWeights);
    }
}

void AttachNekBasicMetadata(vtkDataSet* dataSet,
                            const Nek5000SpectralMetadata& meta) {
    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(meta.nx, meta.ny, meta.nz);
    dataSet->GetFieldData()->AddArray(dimensions);

    vtkNew<vtkIntArray> count;
    count->SetName(kNek5000ElementCountField);
    count->InsertNextValue(static_cast<int>(meta.elementCount));
    dataSet->GetFieldData()->AddArray(count);

    vtkNew<vtkIntArray> ids;
    ids->SetName(kNek5000ElementIdsField);
    ids->SetNumberOfTuples(static_cast<vtkIdType>(meta.elementIds.size()));
    for (std::size_t i = 0; i < meta.elementIds.size(); ++i) {
        ids->SetValue(static_cast<vtkIdType>(i), meta.elementIds[i]);
    }
    dataSet->GetFieldData()->AddArray(ids);

    vtkNew<vtkStringArray> representation;
    representation->SetName(kNek5000RepresentationField);
    representation->InsertNextValue(
        "element-local GLL nodal values; interfaces intentionally unmerged");
    dataSet->GetFieldData()->AddArray(representation);
}

vtkSmartPointer<vtkDataSet> ReadNek5000(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    const NekDescriptor descriptor = ParseNekDescriptor(frame.path);
    if (frame.sourceIndex < 0 || frame.sourceIndex >= descriptor.numberOfTimeSteps) {
        throw std::runtime_error("Nek5000 frame index is out of range.");
    }
    const std::vector<NekFieldMeta> dataPieces =
        ReadNekPieceSet(descriptor, frame.sourceIndex);
    const auto piecesHaveCoordinates = [](const std::vector<NekFieldMeta>& pieces) {
        return !pieces.empty() &&
               std::all_of(pieces.begin(), pieces.end(),
                           [](const NekFieldMeta& piece) {
                               return piece.code.find('X') != std::string::npos;
                           });
    };
    const bool currentHasCoordinates = piecesHaveCoordinates(dataPieces);
    std::vector<NekFieldMeta> coordinatePieces;
    if (currentHasCoordinates) {
        coordinatePieces = dataPieces;
    } else if (frame.coordinateSourceIndex >= 0) {
        if (frame.coordinateSourceIndex >= descriptor.numberOfTimeSteps) {
            throw std::runtime_error(
                "Nek5000 coordinate frame index is out of range.");
        }
        coordinatePieces =
            ReadNekPieceSet(descriptor, frame.coordinateSourceIndex);
        if (!piecesHaveCoordinates(coordinatePieces)) {
            throw std::runtime_error(
                "Enumerated Nek5000 coordinate frame has no mesh block.");
        }
    } else {
        std::vector<int> candidates;
        candidates.reserve(static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        for (int candidate = frame.sourceIndex; candidate >= 0; --candidate) {
            candidates.push_back(candidate);
        }
        for (int candidate = frame.sourceIndex + 1;
             candidate < descriptor.numberOfTimeSteps; ++candidate) {
            candidates.push_back(candidate);
        }
        for (const int candidate : candidates) {
            std::vector<NekFieldMeta> pieces =
                ReadNekPieceSet(descriptor, candidate);
            if (piecesHaveCoordinates(pieces)) {
                coordinatePieces = std::move(pieces);
                break;
            }
        }
    }
    if (coordinatePieces.empty()) {
        throw std::runtime_error(
            "Nek5000 sequence does not contain a mesh-bearing field frame.");
    }
    ValidateNekPieceCompatible(dataPieces.front(), coordinatePieces.front(), false);

    const NekPieceLookup lookup = BuildNekPieceLookup(dataPieces);
    Nek5000SpectralMetadata metadata;
    metadata.nx = dataPieces.front().nx;
    metadata.ny = dataPieces.front().ny;
    metadata.nz = dataPieces.front().nz;
    metadata.elementCount = lookup.outputElementIds.size();
    metadata.elementIds = lookup.outputElementIds;
    if (metadata.pointCount() >
        static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
        throw std::runtime_error("Nek5000 point count exceeds VTK limits.");
    }

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    ScatterNekCoordinateComponent(coordinatePieces, lookup, "x", &x);
    ScatterNekCoordinateComponent(coordinatePieces, lookup, "y", &y);
    if (metadata.is3D()) {
        ScatterNekCoordinateComponent(coordinatePieces, lookup, "z", &z);
    } else {
        z.assign(metadata.pointCount(), 0.0);
    }

    vtkNew<vtkDoubleArray> pointArray;
    pointArray->SetName("Points");
    pointArray->SetNumberOfComponents(3);
    pointArray->SetNumberOfTuples(static_cast<vtkIdType>(metadata.pointCount()));
    double* tuples = pointArray->GetPointer(0);
    for (std::size_t i = 0; i < metadata.pointCount(); ++i) {
        tuples[3 * i + 0] = x[i];
        tuples[3 * i + 1] = y[i];
        tuples[3 * i + 2] = z[i];
    }
    vtkNew<vtkPoints> points;
    points->SetData(pointArray);
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);

    if (options.buildNek5000LinearizedTopology) {
        const std::size_t cellsPerElement =
            static_cast<std::size_t>(metadata.nx - 1) *
            static_cast<std::size_t>(metadata.ny - 1) *
            static_cast<std::size_t>(metadata.is3D() ? metadata.nz - 1 : 1);
        const std::size_t cellCount = metadata.elementCount * cellsPerElement;
        if (cellCount >
            static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
            throw std::runtime_error(
                "Nek5000 linearized cell count exceeds VTK limits.");
        }
        grid->Allocate(static_cast<vtkIdType>(cellCount));
        vtkNew<vtkIntArray> spectralIds;
        spectralIds->SetName(kNek5000SpectralElementIdArray);
        spectralIds->SetNumberOfTuples(static_cast<vtkIdType>(cellCount));
        vtkIdType cellIndex = 0;
        const std::size_t nxyz = metadata.pointsPerElement();
        for (std::size_t element = 0; element < metadata.elementCount; ++element) {
            const vtkIdType base = static_cast<vtkIdType>(element * nxyz);
            if (metadata.is3D()) {
                for (int k = 0; k < metadata.nz - 1; ++k) {
                    for (int j = 0; j < metadata.ny - 1; ++j) {
                        for (int i = 0; i < metadata.nx - 1; ++i) {
                            vtkIdType ids[8] = {
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j, k)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j, k)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j + 1, k)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j + 1, k)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j, k + 1)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j, k + 1)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j + 1, k + 1)),
                                base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j + 1, k + 1)),
                            };
                            grid->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
                            spectralIds->SetValue(cellIndex++,
                                                  metadata.elementIds[element]);
                        }
                    }
                }
            } else {
                for (int j = 0; j < metadata.ny - 1; ++j) {
                    for (int i = 0; i < metadata.nx - 1; ++i) {
                        vtkIdType ids[4] = {
                            base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j, 0)),
                            base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j, 0)),
                            base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i + 1, j + 1, 0)),
                            base + static_cast<vtkIdType>(NekLocalIndex(metadata.nx, metadata.ny, i, j + 1, 0)),
                        };
                        grid->InsertNextCell(VTK_QUAD, 4, ids);
                        spectralIds->SetValue(cellIndex++,
                                              metadata.elementIds[element]);
                    }
                }
            }
        }
        grid->GetCellData()->AddArray(spectralIds);
    }

    if (options.readAllArrays) {
        const std::string code = dataPieces.front().code;
        if (code.find('U') != std::string::npos) {
            std::vector<double> u;
            std::vector<double> v;
            std::vector<double> w;
            ScatterNekComponent(dataPieces, lookup, "u", &u);
            ScatterNekComponent(dataPieces, lookup, "v", &v);
            if (metadata.is3D()) {
                ScatterNekComponent(dataPieces, lookup, "w", &w);
            } else {
                w.assign(metadata.pointCount(), 0.0);
            }
            vtkNew<vtkDoubleArray> velocity;
            velocity->SetName("Velocity");
            velocity->SetNumberOfComponents(3);
            velocity->SetNumberOfTuples(static_cast<vtkIdType>(metadata.pointCount()));
            for (std::size_t i = 0; i < metadata.pointCount(); ++i) {
                const double tuple[3] = {u[i], v[i], w[i]};
                velocity->SetTuple(static_cast<vtkIdType>(i), tuple);
            }
            grid->GetPointData()->AddArray(velocity);
        }

        const auto addScalar = [&](const std::string& label, const std::string& name) {
            std::vector<double> values;
            ScatterNekComponent(dataPieces, lookup, label, &values);
            vtkNew<vtkDoubleArray> array;
            array->SetName(name.c_str());
            array->SetNumberOfTuples(static_cast<vtkIdType>(values.size()));
            std::copy(values.begin(), values.end(), array->GetPointer(0));
            grid->GetPointData()->AddArray(array);
        };
        if (code.find('P') != std::string::npos) {
            addScalar("p", "Pressure");
        }
        if (code.find('T') != std::string::npos) {
            addScalar("t", "Temperature");
        }
        for (const std::string& label : dataPieces.front().labels) {
            if (!label.empty() && label.front() == 's') {
                std::string name = label;
                name.front() = 'S';
                addScalar(label, name);
            }
        }
    }

    vtkNew<vtkDoubleArray> time;
    time->SetName(kTimeValueArray);
    time->InsertNextValue(dataPieces.front().time);
    grid->GetFieldData()->AddArray(time);

    if (options.attachNek5000SpectralMetadata) {
        AttachNekBasicMetadata(grid, metadata);
        if (options.attachNek5000QuadratureArrays) {
            AttachNekQuadratureArrays(
                grid,
                ComputeNekGeometryMetrics(
                    grid, metadata,
                    options.attachNek5000ReferenceWeightArray,
                    options.attachNek5000JacobianDeterminantArray,
                    options.attachNek5000QuadratureWeightArray));
        }
    }
    return grid;
}

std::vector<DataSetFrame> ParsePvd(const std::filesystem::path& path,
                                   double fallbackPhysicalDt) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Cannot open PVD file: " + path.string());
    }
    std::ostringstream text;
    text << input.rdbuf();
    const std::regex tagPattern(R"(<\s*DataSet\b[^>]*>)", std::regex::icase);
    const std::regex filePattern(R"(file\s*=\s*(['"])(.*?)\1)", std::regex::icase);
    const std::regex timePattern(R"(timestep\s*=\s*(['"])(.*?)\1)", std::regex::icase);
    const std::string xml = text.str();
    std::vector<DataSetFrame> frames;
    for (std::sregex_iterator it(xml.begin(), xml.end(), tagPattern), end;
         it != end; ++it) {
        const std::string tag = it->str();
        std::smatch fileMatch;
        if (!std::regex_search(tag, fileMatch, filePattern)) {
            continue;
        }
        std::filesystem::path framePath(fileMatch[2].str());
        if (framePath.is_relative()) {
            framePath = path.parent_path() / framePath;
        }
        double time = static_cast<double>(frames.size()) * fallbackPhysicalDt;
        std::smatch timeMatch;
        if (std::regex_search(tag, timeMatch, timePattern)) {
            time = std::stod(timeMatch[2].str());
        }
        std::vector<DataSetFrame> expanded =
            EnumerateDataSetFrames(framePath, fallbackPhysicalDt);
        if (expanded.size() == 1) {
            expanded.front().time = time;
            if (expanded.front().kind != DataSetFrameKind::File) {
                expanded.front().readerTime = time;
            }
            frames.push_back(std::move(expanded.front()));
        } else {
            frames.insert(frames.end(), expanded.begin(), expanded.end());
        }
    }
    return frames;
}

long long LastNumericTag(const std::string& name) {
    const std::regex numbers(R"((\d+))");
    long long tag = -1;
    for (std::sregex_iterator it(name.begin(), name.end(), numbers), end;
         it != end; ++it) {
        tag = std::stoll((*it)[1].str());
    }
    return tag;
}

std::vector<DataSetFrame> EnumerateFolder(const std::filesystem::path& folder,
                                          double fallbackPhysicalDt) {
    struct TaggedFrame {
        long long tag = -1;
        std::string name;
        std::size_t localIndex = 0;
        DataSetFrame frame;
    };
    std::vector<TaggedFrame> tagged;
    std::unordered_set<std::string> nek5000Stems;
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
        if (entry.is_regular_file() &&
            LowerExtension(entry.path()) == ".nek5000") {
            nek5000Stems.insert(ToLower(entry.path().stem().string()));
        }
    }
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string extension = LowerExtension(entry.path());
        if (extension == ".nek" &&
            nek5000Stems.find(ToLower(entry.path().stem().string())) !=
                nek5000Stems.end()) {
            continue;
        }
        if (!IsOrdinaryVtkExtension(extension) && extension != ".nek5000" &&
            extension != ".nek" &&
            extension != ".case") {
            continue;
        }
        const std::vector<DataSetFrame> local =
            EnumerateDataSetFrames(entry.path(), fallbackPhysicalDt);
        for (std::size_t i = 0; i < local.size(); ++i) {
            tagged.push_back({LastNumericTag(entry.path().filename().string()),
                              entry.path().filename().string(), i, local[i]});
        }
    }
    std::sort(tagged.begin(), tagged.end(), [](const TaggedFrame& lhs,
                                               const TaggedFrame& rhs) {
        if (lhs.tag != rhs.tag) {
            return lhs.tag < rhs.tag;
        }
        if (lhs.name != rhs.name) {
            return lhs.name < rhs.name;
        }
        return lhs.localIndex < rhs.localIndex;
    });
    std::vector<DataSetFrame> frames;
    frames.reserve(tagged.size());
    for (std::size_t i = 0; i < tagged.size(); ++i) {
        DataSetFrame frame = std::move(tagged[i].frame);
        if (!frame.readerTime.has_value() && frame.kind == DataSetFrameKind::File) {
            frame.time = static_cast<double>(i) * fallbackPhysicalDt;
        }
        frames.push_back(std::move(frame));
    }
    return frames;
}

}  // namespace

struct Nek5000SpectralPointEvaluator::Impl {
    struct ElementBvhNode {
        std::array<double, 6> bounds{};
        double tolerance = 0.0;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = std::numeric_limits<std::size_t>::max();
        std::size_t right = std::numeric_limits<std::size_t>::max();

        bool IsLeaf() const noexcept {
            return left == std::numeric_limits<std::size_t>::max();
        }
    };

    vtkSmartPointer<vtkDataSet> dataSet;
    vtkSmartPointer<vtkDataArray> pointArray;
    Nek5000SpectralMetadata metadata;
    Nek5000PointEvaluationOptions options;
    int components = 0;
    std::array<std::vector<double>, 3> nodes;
    std::array<std::vector<double>, 3> barycentricWeights;
    std::array<std::vector<double>, 3> derivativeMatrices;
    vtkSmartPointer<vtkDataArray> coordinates;
    std::vector<std::array<double, 6>> elementBounds;
    std::vector<std::size_t> bvhElements;
    std::vector<ElementBvhNode> bvhNodes;
    std::size_t bvhRoot = std::numeric_limits<std::size_t>::max();
    vtkSmartPointer<vtkStaticCellLocator> cellLocator;
    std::size_t cellsPerElement = 0;
    std::shared_ptr<const Impl> geometryOwner;

    Impl(vtkDataSet* inputDataSet,
         vtkDataArray* inputArray,
         const Nek5000PointEvaluationOptions& inputOptions)
        : options(inputOptions) {
        if (inputDataSet == nullptr || inputArray == nullptr ||
            inputArray->GetNumberOfTuples() != inputDataSet->GetNumberOfPoints()) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation requires a complete point-data array.");
        }
        if (!std::isfinite(options.referenceTolerance) ||
            options.referenceTolerance < 0.0 ||
            !std::isfinite(options.relativePhysicalTolerance) ||
            options.relativePhysicalTolerance < 0.0 ||
            !std::isfinite(options.absolutePhysicalTolerance) ||
            options.absolutePhysicalTolerance < 0.0 ||
            options.maxNewtonIterations < 1) {
            throw std::runtime_error(
                "Invalid Nek5000 spectral point-evaluation options.");
        }
        dataSet = inputDataSet;
        pointArray = inputArray;
        metadata = GetNek5000SpectralMetadata(inputDataSet);
        components = inputArray->GetNumberOfComponents();
        if (components < 1) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation input has no components.");
        }

        nodes[0] = GllNodes(metadata.nx);
        nodes[1] = GllNodes(metadata.ny);
        nodes[2] = metadata.is3D() ? GllNodes(metadata.nz)
                                   : std::vector<double>{0.0};
        for (int axis = 0; axis < 3; ++axis) {
            barycentricWeights[static_cast<std::size_t>(axis)] =
                BarycentricWeights(nodes[static_cast<std::size_t>(axis)]);
            derivativeMatrices[static_cast<std::size_t>(axis)] =
                nodes[static_cast<std::size_t>(axis)].size() == 1
                    ? std::vector<double>{0.0}
                    : DerivativeMatrix(nodes[static_cast<std::size_t>(axis)]);
        }

        vtkPointSet* pointSet = vtkPointSet::SafeDownCast(inputDataSet);
        if (pointSet == nullptr || pointSet->GetPoints() == nullptr ||
            pointSet->GetPoints()->GetData() == nullptr) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation requires explicit point coordinates.");
        }
        coordinates = pointSet->GetPoints()->GetData();

        const std::size_t nxyz = metadata.pointsPerElement();
        const std::array<std::vector<double>, 3> unitNodes = {
            [&]() {
                std::vector<double> values(nodes[0].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[0][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
            [&]() {
                std::vector<double> values(nodes[1].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[1][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
            [&]() {
                if (!metadata.is3D()) {
                    return std::vector<double>{0.0};
                }
                std::vector<double> values(nodes[2].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[2][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
        };
        const double lebesgueProduct =
            ConservativeLebesgueBound(nodes[0]) *
            ConservativeLebesgueBound(nodes[1]) *
            (metadata.is3D() ? ConservativeLebesgueBound(nodes[2]) : 1.0);
        if (!std::isfinite(lebesgueProduct)) {
            throw std::runtime_error(
                "Cannot construct finite conservative Nek5000 element bounds.");
        }

        elementBounds.resize(metadata.elementCount);
        std::atomic<bool> invalidCoordinates{false};
        vtkSMPTools::For(
            static_cast<vtkIdType>(0),
            static_cast<vtkIdType>(metadata.elementCount),
            static_cast<vtkIdType>(32),
            [&](vtkIdType begin, vtkIdType end) {
                const std::size_t nxy =
                    static_cast<std::size_t>(metadata.nx) *
                    static_cast<std::size_t>(metadata.ny);
                std::vector<std::array<double, 3>> lowerFaces(nxy);
                std::vector<std::array<double, 3>> upperFaces(nxy);
                const auto interpolate = [](double lower, double upper,
                                            double parameter) {
                    if (parameter == 0.0) {
                        return lower;
                    }
                    if (parameter == 1.0) {
                        return upper;
                    }
                    return std::fma(parameter, upper - lower, lower);
                };

                for (vtkIdType elementId = begin; elementId < end; ++elementId) {
                    const std::size_t element =
                        static_cast<std::size_t>(elementId);
                    const std::size_t base = element * nxyz;
                    std::array<std::array<double, 3>, 8> corners{};
                    bool elementValid = true;
                    for (int cornerK = 0; cornerK < 2; ++cornerK) {
                        const int k = metadata.is3D() && cornerK != 0
                                          ? metadata.nz - 1
                                          : 0;
                        for (int cornerJ = 0; cornerJ < 2; ++cornerJ) {
                            const int j = cornerJ == 0 ? 0 : metadata.ny - 1;
                            for (int cornerI = 0; cornerI < 2; ++cornerI) {
                                const int i = cornerI == 0 ? 0 : metadata.nx - 1;
                                const std::size_t corner =
                                    static_cast<std::size_t>(cornerI +
                                                             2 * cornerJ +
                                                             4 * cornerK);
                                const vtkIdType pointIndex = static_cast<vtkIdType>(
                                    base + NekLocalIndex(metadata.nx, metadata.ny,
                                                         i, j, k));
                                coordinates->GetTuple(pointIndex,
                                                      corners[corner].data());
                                for (double coordinate : corners[corner]) {
                                    elementValid = elementValid &&
                                                   std::isfinite(coordinate);
                                }
                            }
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }

                    std::array<std::array<double, 3>, 8> localCorners{};
                    std::array<double, 3> cornerMinimum = {
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
                    std::array<double, 3> cornerMaximum = {
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
                    for (std::size_t corner = 0; corner < corners.size(); ++corner) {
                        for (int axis = 0; axis < 3; ++axis) {
                            const std::size_t component =
                                static_cast<std::size_t>(axis);
                            localCorners[corner][component] =
                                corners[corner][component] - corners[0][component];
                            cornerMinimum[component] =
                                std::min(cornerMinimum[component],
                                         corners[corner][component]);
                            cornerMaximum[component] =
                                std::max(cornerMaximum[component],
                                         corners[corner][component]);
                        }
                    }

                    for (int j = 0; j < metadata.ny; ++j) {
                        const double v = unitNodes[1][static_cast<std::size_t>(j)];
                        for (int i = 0; i < metadata.nx; ++i) {
                            const double u =
                                unitNodes[0][static_cast<std::size_t>(i)];
                            const std::size_t faceIndex =
                                static_cast<std::size_t>(j) *
                                    static_cast<std::size_t>(metadata.nx) +
                                static_cast<std::size_t>(i);
                            for (int axis = 0; axis < 3; ++axis) {
                                const std::size_t component =
                                    static_cast<std::size_t>(axis);
                                const double lowerLow = interpolate(
                                    localCorners[0][component],
                                    localCorners[1][component], u);
                                const double lowerHigh = interpolate(
                                    localCorners[2][component],
                                    localCorners[3][component], u);
                                const double upperLow = interpolate(
                                    localCorners[4][component],
                                    localCorners[5][component], u);
                                const double upperHigh = interpolate(
                                    localCorners[6][component],
                                    localCorners[7][component], u);
                                lowerFaces[faceIndex][component] =
                                    interpolate(lowerLow, lowerHigh, v);
                                upperFaces[faceIndex][component] =
                                    interpolate(upperLow, upperHigh, v);
                            }
                        }
                    }

                    std::array<double, 3> maximumResidual{};
                    std::array<double, 3> maximumLocalMagnitude{};
                    std::array<double, 3> nodeMinimum = {
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
                    std::array<double, 3> nodeMaximum = {
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
                    for (int k = 0; k < metadata.nz && elementValid; ++k) {
                        const double w = metadata.is3D()
                                             ? unitNodes[2][static_cast<std::size_t>(k)]
                                             : 0.0;
                        for (int j = 0; j < metadata.ny && elementValid; ++j) {
                            for (int i = 0; i < metadata.nx; ++i) {
                                const std::size_t faceIndex =
                                    static_cast<std::size_t>(j) *
                                        static_cast<std::size_t>(metadata.nx) +
                                    static_cast<std::size_t>(i);
                                const std::size_t local = NekLocalIndex(
                                    metadata.nx, metadata.ny, i, j, k);
                                std::array<double, 3> point{};
                                coordinates->GetTuple(
                                    static_cast<vtkIdType>(base + local),
                                    point.data());
                                for (int axis = 0; axis < 3; ++axis) {
                                    const std::size_t component =
                                        static_cast<std::size_t>(axis);
                                    const double localPoint =
                                        point[component] - corners[0][component];
                                    const double trilinearPoint = interpolate(
                                        lowerFaces[faceIndex][component],
                                        upperFaces[faceIndex][component], w);
                                    if (!std::isfinite(point[component]) ||
                                        !std::isfinite(localPoint) ||
                                        !std::isfinite(trilinearPoint)) {
                                        elementValid = false;
                                        break;
                                    }
                                    nodeMinimum[component] =
                                        std::min(nodeMinimum[component],
                                                 point[component]);
                                    nodeMaximum[component] =
                                        std::max(nodeMaximum[component],
                                                 point[component]);
                                    maximumLocalMagnitude[component] = std::max(
                                        maximumLocalMagnitude[component],
                                        std::abs(localPoint));
                                    maximumResidual[component] = std::max(
                                        maximumResidual[component],
                                        std::abs(localPoint - trilinearPoint));
                                }
                            }
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }

                    std::array<double, 6> bounds{};
                    // The residual path contains three nested affine
                    // interpolations plus origin shifts. This factor is far
                    // above the corresponding floating-point gamma bound, so
                    // the value added below also encloses their roundoff.
                    constexpr double kResidualRoundoffFactor = 4096.0;
                    for (int axis = 0; axis < 3; ++axis) {
                        const std::size_t component =
                            static_cast<std::size_t>(axis);
                        const double residualPadding =
                            kResidualRoundoffFactor *
                            std::numeric_limits<double>::epsilon() *
                            maximumLocalMagnitude[component];
                        const double safeResidual = RoundUp(
                            maximumResidual[component] + residualPadding);
                        const double margin = RoundUp(
                            lebesgueProduct * safeResidual);
                        const double lower = RoundDown(
                            cornerMinimum[component] - margin);
                        const double upper = RoundUp(
                            cornerMaximum[component] + margin);
                        bounds[2 * component] =
                            std::min(lower, nodeMinimum[component]);
                        bounds[2 * component + 1] =
                            std::max(upper, nodeMaximum[component]);
                        if (!std::isfinite(bounds[2 * component]) ||
                            !std::isfinite(bounds[2 * component + 1])) {
                            elementValid = false;
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }
                    elementBounds[element] = bounds;
                }
            });
        if (invalidCoordinates.load(std::memory_order_relaxed)) {
            throw std::runtime_error(
                "Nek5000 spectral coordinates must be finite and bounded.");
        }
        BuildElementBvh();

        cellsPerElement =
            static_cast<std::size_t>(metadata.nx - 1) *
            static_cast<std::size_t>(metadata.ny - 1) *
            static_cast<std::size_t>(metadata.is3D() ? metadata.nz - 1 : 1);
        if (cellsPerElement > 0 &&
            static_cast<std::size_t>(inputDataSet->GetNumberOfCells()) ==
                metadata.elementCount * cellsPerElement) {
            cellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
            cellLocator->SetDataSet(inputDataSet);
            cellLocator->BuildLocator();
        }
    }

    Impl(std::shared_ptr<const Impl> inputGeometry,
         vtkDataArray* inputArray)
        : geometryOwner(std::move(inputGeometry)) {
        if (geometryOwner == nullptr || inputArray == nullptr ||
            inputArray->GetNumberOfTuples() !=
                geometryOwner->dataSet->GetNumberOfPoints()) {
            throw std::runtime_error(
                "Rebound Nek5000 spectral field has an incompatible point count.");
        }
        pointArray = inputArray;
        components = inputArray->GetNumberOfComponents();
        options = geometryOwner->options;
        if (components < 1) {
            throw std::runtime_error(
                "Rebound Nek5000 spectral field has no components.");
        }
    }

    const Impl& Geometry() const {
        return geometryOwner == nullptr ? *this : *geometryOwner;
    }

    double PhysicalTolerance(std::size_t element) const {
        const auto& bounds = Geometry().elementBounds[element];
        const double dx = bounds[1] - bounds[0];
        const double dy = bounds[3] - bounds[2];
        const double dz = bounds[5] - bounds[4];
        return std::max(options.absolutePhysicalTolerance,
                        options.relativePhysicalTolerance *
                            std::sqrt(dx * dx + dy * dy + dz * dz));
    }

    bool BoundsContain(std::size_t element, const double point[3]) const {
        const auto& bounds = Geometry().elementBounds[element];
        const double tolerance = PhysicalTolerance(element);
        for (int axis = 0; axis < 3; ++axis) {
            if (point[axis] < bounds[2 * static_cast<std::size_t>(axis)] - tolerance ||
                point[axis] > bounds[2 * static_cast<std::size_t>(axis) + 1] + tolerance) {
                return false;
            }
        }
        return true;
    }

    std::size_t BuildElementBvhNode(std::size_t begin, std::size_t end) {
        ElementBvhNode node;
        node.bounds = {
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        std::array<double, 3> minimumCenter = {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
        };
        std::array<double, 3> maximumCenter = {
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        for (std::size_t offset = begin; offset < end; ++offset) {
            const std::size_t element = bvhElements[offset];
            const auto& bounds = elementBounds[element];
            node.tolerance = std::max(node.tolerance,
                                      PhysicalTolerance(element));
            for (int axis = 0; axis < 3; ++axis) {
                const std::size_t lo = 2 * static_cast<std::size_t>(axis);
                const std::size_t hi = lo + 1;
                node.bounds[lo] = std::min(node.bounds[lo], bounds[lo]);
                node.bounds[hi] = std::max(node.bounds[hi], bounds[hi]);
                const double center = 0.5 * (bounds[lo] + bounds[hi]);
                minimumCenter[static_cast<std::size_t>(axis)] =
                    std::min(minimumCenter[static_cast<std::size_t>(axis)], center);
                maximumCenter[static_cast<std::size_t>(axis)] =
                    std::max(maximumCenter[static_cast<std::size_t>(axis)], center);
            }
        }

        const std::size_t nodeIndex = bvhNodes.size();
        bvhNodes.push_back(node);
        constexpr std::size_t kLeafElementCount = 8;
        if (end - begin <= kLeafElementCount) {
            bvhNodes[nodeIndex].begin = begin;
            bvhNodes[nodeIndex].end = end;
            return nodeIndex;
        }

        int splitAxis = 0;
        double splitExtent = maximumCenter[0] - minimumCenter[0];
        for (int axis = 1; axis < 3; ++axis) {
            const double extent = maximumCenter[static_cast<std::size_t>(axis)] -
                                  minimumCenter[static_cast<std::size_t>(axis)];
            if (extent > splitExtent) {
                splitExtent = extent;
                splitAxis = axis;
            }
        }
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(
            bvhElements.begin() + static_cast<std::ptrdiff_t>(begin),
            bvhElements.begin() + static_cast<std::ptrdiff_t>(middle),
            bvhElements.begin() + static_cast<std::ptrdiff_t>(end),
            [&](std::size_t lhs, std::size_t rhs) {
                const std::size_t lo =
                    2 * static_cast<std::size_t>(splitAxis);
                const double lhsCenter =
                    0.5 * (elementBounds[lhs][lo] + elementBounds[lhs][lo + 1]);
                const double rhsCenter =
                    0.5 * (elementBounds[rhs][lo] + elementBounds[rhs][lo + 1]);
                return lhsCenter < rhsCenter ||
                       (lhsCenter == rhsCenter && lhs < rhs);
            });
        const std::size_t left = BuildElementBvhNode(begin, middle);
        const std::size_t right = BuildElementBvhNode(middle, end);
        bvhNodes[nodeIndex].left = left;
        bvhNodes[nodeIndex].right = right;
        return nodeIndex;
    }

    void BuildElementBvh() {
        bvhElements.resize(metadata.elementCount);
        std::iota(bvhElements.begin(), bvhElements.end(), 0);
        if (bvhElements.empty()) {
            return;
        }
        if (bvhElements.size() <=
            std::numeric_limits<std::size_t>::max() / 2) {
            bvhNodes.reserve(2 * bvhElements.size());
        }
        bvhRoot = BuildElementBvhNode(0, bvhElements.size());
    }

    bool BvhNodeContains(const ElementBvhNode& node,
                         const double point[3]) const {
        for (int axis = 0; axis < 3; ++axis) {
            const std::size_t lo = 2 * static_cast<std::size_t>(axis);
            if (point[axis] < node.bounds[lo] - node.tolerance ||
                point[axis] > node.bounds[lo + 1] + node.tolerance) {
                return false;
            }
        }
        return true;
    }

    void CollectBvhCandidates(const double point[3],
                              std::vector<std::size_t>* candidates) const {
        if (candidates == nullptr) {
            return;
        }
        candidates->clear();
        const Impl& geometry = Geometry();
        if (geometry.bvhRoot == std::numeric_limits<std::size_t>::max()) {
            return;
        }
        thread_local std::vector<std::size_t> stack;
        stack.clear();
        stack.push_back(geometry.bvhRoot);
        while (!stack.empty()) {
            const std::size_t nodeIndex = stack.back();
            stack.pop_back();
            const ElementBvhNode& node = geometry.bvhNodes[nodeIndex];
            if (!geometry.BvhNodeContains(node, point)) {
                continue;
            }
            if (node.IsLeaf()) {
                candidates->insert(
                    candidates->end(),
                    geometry.bvhElements.begin() +
                        static_cast<std::ptrdiff_t>(node.begin),
                    geometry.bvhElements.begin() +
                        static_cast<std::ptrdiff_t>(node.end));
                continue;
            }
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }

    void EvaluateTensor(std::size_t element,
                        const std::array<double, 3>& reference,
                        std::array<double, 3>* mapped,
                        std::array<double, 3>* xr,
                        std::array<double, 3>* xs,
                        std::array<double, 3>* xt,
                        double* tuple) const {
        const Impl& geometry = Geometry();
        std::array<std::vector<double>, 3> basis;
        std::array<std::vector<double>, 3> derivative;
        for (int axis = 0; axis < 3; ++axis) {
            EvaluateNodalBasis(geometry.nodes[static_cast<std::size_t>(axis)],
                               geometry.barycentricWeights[static_cast<std::size_t>(axis)],
                               geometry.derivativeMatrices[static_cast<std::size_t>(axis)],
                               reference[static_cast<std::size_t>(axis)],
                               &basis[static_cast<std::size_t>(axis)],
                               &derivative[static_cast<std::size_t>(axis)]);
        }
        if (mapped != nullptr) {
            mapped->fill(0.0);
        }
        if (xr != nullptr) {
            xr->fill(0.0);
        }
        if (xs != nullptr) {
            xs->fill(0.0);
        }
        if (xt != nullptr) {
            xt->fill(0.0);
        }
        if (tuple != nullptr) {
            std::fill_n(tuple, components, 0.0);
        }

        const std::size_t base = element * geometry.metadata.pointsPerElement();
        for (int k = 0; k < geometry.metadata.nz; ++k) {
            for (int j = 0; j < geometry.metadata.ny; ++j) {
                for (int i = 0; i < geometry.metadata.nx; ++i) {
                    const double lx = basis[0][static_cast<std::size_t>(i)];
                    const double ly = basis[1][static_cast<std::size_t>(j)];
                    const double lz = basis[2][static_cast<std::size_t>(k)];
                    const double weight = lx * ly * lz;
                    const double weightR =
                        derivative[0][static_cast<std::size_t>(i)] * ly * lz;
                    const double weightS =
                        lx * derivative[1][static_cast<std::size_t>(j)] * lz;
                    const double weightT =
                        lx * ly * derivative[2][static_cast<std::size_t>(k)];
                    const std::size_t local =
                        NekLocalIndex(geometry.metadata.nx, geometry.metadata.ny,
                                      i, j, k);
                    const std::size_t pointIndex = base + local;
                    for (int axis = 0; axis < 3; ++axis) {
                        const double coordinate =
                            geometry.coordinates->GetComponent(
                                static_cast<vtkIdType>(pointIndex), axis);
                        if (mapped != nullptr) {
                            (*mapped)[static_cast<std::size_t>(axis)] +=
                                weight * coordinate;
                        }
                        if (xr != nullptr) {
                            (*xr)[static_cast<std::size_t>(axis)] +=
                                weightR * coordinate;
                        }
                        if (xs != nullptr) {
                            (*xs)[static_cast<std::size_t>(axis)] +=
                                weightS * coordinate;
                        }
                        if (xt != nullptr) {
                            (*xt)[static_cast<std::size_t>(axis)] +=
                                weightT * coordinate;
                        }
                    }
                    if (tuple != nullptr) {
                        for (int component = 0; component < components; ++component) {
                            tuple[component] +=
                                weight * pointArray->GetComponent(
                                             static_cast<vtkIdType>(pointIndex), component);
                        }
                    }
                }
            }
        }
    }

    static double Dot(const std::array<double, 3>& lhs,
                      const std::array<double, 3>& rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    }

    static std::array<double, 3> Cross(const std::array<double, 3>& lhs,
                                       const std::array<double, 3>& rhs) {
        return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
                lhs[2] * rhs[0] - lhs[0] * rhs[2],
                lhs[0] * rhs[1] - lhs[1] * rhs[0]};
    }

    bool NewtonFrom(std::size_t element,
                    const double point[3],
                    std::array<double, 3> reference,
                    std::array<double, 3>* result,
                    double* finalResidual) const {
        const Nek5000SpectralMetadata& geometryMetadata = Geometry().metadata;
        const double tolerance = PhysicalTolerance(element);
        const double toleranceSquared = tolerance * tolerance;
        double residualSquared = std::numeric_limits<double>::infinity();
        for (int iteration = 0; iteration < options.maxNewtonIterations;
             ++iteration) {
            std::array<double, 3> mapped{};
            std::array<double, 3> xr{};
            std::array<double, 3> xs{};
            std::array<double, 3> xt{};
            EvaluateTensor(element, reference, &mapped, &xr, &xs, &xt, nullptr);
            const std::array<double, 3> residual = {
                mapped[0] - point[0], mapped[1] - point[1], mapped[2] - point[2]};
            residualSquared = SquaredNorm3(residual);
            if (residualSquared <= toleranceSquared) {
                break;
            }

            std::array<double, 3> delta{};
            if (geometryMetadata.is3D()) {
                const double determinant = Dot(xr, Cross(xs, xt));
                if (!std::isfinite(determinant) ||
                    std::abs(determinant) <= 1.0e-30) {
                    return false;
                }
                delta[0] = Dot(residual, Cross(xs, xt)) / determinant;
                delta[1] = Dot(xr, Cross(residual, xt)) / determinant;
                delta[2] = Dot(xr, Cross(xs, residual)) / determinant;
            } else {
                const double g00 = Dot(xr, xr);
                const double g01 = Dot(xr, xs);
                const double g11 = Dot(xs, xs);
                const double determinant = g00 * g11 - g01 * g01;
                if (!std::isfinite(determinant) || determinant <= 1.0e-30) {
                    return false;
                }
                const double rhs0 = Dot(xr, residual);
                const double rhs1 = Dot(xs, residual);
                delta[0] = (g11 * rhs0 - g01 * rhs1) / determinant;
                delta[1] = (-g01 * rhs0 + g00 * rhs1) / determinant;
            }

            double scale = 1.0;
            bool accepted = false;
            for (int lineSearch = 0; lineSearch < 10; ++lineSearch) {
                std::array<double, 3> candidate = {
                    reference[0] - scale * delta[0],
                    reference[1] - scale * delta[1],
                    geometryMetadata.is3D() ? reference[2] - scale * delta[2] : 0.0};
                if (std::abs(candidate[0]) > 4.0 ||
                    std::abs(candidate[1]) > 4.0 ||
                    (geometryMetadata.is3D() && std::abs(candidate[2]) > 4.0)) {
                    scale *= 0.5;
                    continue;
                }
                std::array<double, 3> candidateMapped{};
                EvaluateTensor(element, candidate, &candidateMapped, nullptr,
                               nullptr, nullptr, nullptr);
                const std::array<double, 3> candidateResidual = {
                    candidateMapped[0] - point[0],
                    candidateMapped[1] - point[1],
                    candidateMapped[2] - point[2]};
                if (SquaredNorm3(candidateResidual) < residualSquared) {
                    reference = candidate;
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!accepted) {
                return false;
            }
        }

        std::array<double, 3> mapped{};
        EvaluateTensor(element, reference, &mapped, nullptr, nullptr, nullptr,
                       nullptr);
        const std::array<double, 3> residual = {
            mapped[0] - point[0], mapped[1] - point[1], mapped[2] - point[2]};
        residualSquared = SquaredNorm3(residual);
        if (finalResidual != nullptr) {
            *finalResidual = std::sqrt(residualSquared);
        }
        const double referenceLimit = 1.0 + options.referenceTolerance;
        const bool insideReference =
            std::abs(reference[0]) <= referenceLimit &&
            std::abs(reference[1]) <= referenceLimit &&
            (!geometryMetadata.is3D() || std::abs(reference[2]) <= referenceLimit);
        if (!insideReference ||
            residualSquared > PhysicalTolerance(element) * PhysicalTolerance(element)) {
            return false;
        }
        if (result != nullptr) {
            *result = reference;
        }
        return true;
    }

    bool Locate(std::size_t element,
                const double point[3],
                std::array<double, 3>* reference) const {
        const Impl& geometry = Geometry();
        if (!BoundsContain(element, point)) {
            return false;
        }

        const std::size_t nxyz = geometry.metadata.pointsPerElement();
        const std::size_t base = element * nxyz;
        std::size_t nearest = 0;
        double nearestDistance = std::numeric_limits<double>::infinity();
        for (std::size_t local = 0; local < nxyz; ++local) {
            const vtkIdType pointIndex = static_cast<vtkIdType>(base + local);
            const double dx = geometry.coordinates->GetComponent(pointIndex, 0) - point[0];
            const double dy = geometry.coordinates->GetComponent(pointIndex, 1) - point[1];
            const double dz = geometry.coordinates->GetComponent(pointIndex, 2) - point[2];
            const double distance = dx * dx + dy * dy + dz * dz;
            if (distance < nearestDistance) {
                nearestDistance = distance;
                nearest = local;
            }
        }
        const int nearestI = static_cast<int>(nearest %
                                               static_cast<std::size_t>(geometry.metadata.nx));
        const int nearestJ = static_cast<int>(
            (nearest / static_cast<std::size_t>(geometry.metadata.nx)) %
            static_cast<std::size_t>(geometry.metadata.ny));
        const int nearestK = static_cast<int>(
            nearest / (static_cast<std::size_t>(geometry.metadata.nx) *
                       static_cast<std::size_t>(geometry.metadata.ny)));

        std::array<std::array<double, 3>, 3> initial = {
            std::array<double, 3>{0.0, 0.0, 0.0},
            std::array<double, 3>{geometry.nodes[0][static_cast<std::size_t>(nearestI)],
                                  geometry.nodes[1][static_cast<std::size_t>(nearestJ)],
                                  geometry.metadata.is3D()
                                      ? geometry.nodes[2][static_cast<std::size_t>(nearestK)]
                                      : 0.0},
            std::array<double, 3>{0.0, 0.0, 0.0},
        };
        const auto& bounds = geometry.elementBounds[element];
        for (int axis = 0; axis < (geometry.metadata.is3D() ? 3 : 2); ++axis) {
            const double extent =
                bounds[2 * static_cast<std::size_t>(axis) + 1] -
                bounds[2 * static_cast<std::size_t>(axis)];
            if (extent > 0.0) {
                initial[2][static_cast<std::size_t>(axis)] =
                    std::max(-1.0, std::min(
                                       1.0,
                                       2.0 * (point[axis] -
                                              bounds[2 * static_cast<std::size_t>(axis)]) /
                                               extent -
                                           1.0));
            }
        }

        bool found = false;
        double bestResidual = std::numeric_limits<double>::infinity();
        std::array<double, 3> best{};
        for (const auto& guess : initial) {
            std::array<double, 3> candidate{};
            double residual = 0.0;
            if (NewtonFrom(element, point, guess, &candidate, &residual) &&
                residual < bestResidual) {
                bestResidual = residual;
                best = candidate;
                found = true;
            }
        }
        if (found && reference != nullptr) {
            *reference = best;
        }
        return found;
    }

    bool Evaluate(const double point[3], double* tuple, int* elementHint) const {
        if (point == nullptr || tuple == nullptr || !std::isfinite(point[0]) ||
            !std::isfinite(point[1]) || !std::isfinite(point[2])) {
            return false;
        }
        const Impl& geometry = Geometry();
        const int hint = elementHint == nullptr ? -1 : *elementHint;
        const auto tryElement = [&](std::size_t element) {
            std::array<double, 3> reference{};
            if (!Locate(element, point, &reference)) {
                return false;
            }
            EvaluateTensor(element, reference, nullptr, nullptr, nullptr, nullptr,
                           tuple);
            if (elementHint != nullptr) {
                *elementHint = static_cast<int>(element);
            }
            return true;
        };

        if (hint >= 0 &&
            static_cast<std::size_t>(hint) < geometry.metadata.elementCount &&
            tryElement(static_cast<std::size_t>(hint))) {
            return true;
        }
        int locatorElement = -1;
        if (geometry.cellLocator != nullptr) {
            thread_local vtkSmartPointer<vtkGenericCell> locatorCell =
                vtkSmartPointer<vtkGenericCell>::New();
            thread_local std::vector<double> locatorWeights;
            locatorWeights.resize(static_cast<std::size_t>(
                std::max<vtkIdType>(1, geometry.dataSet->GetMaxCellSize())));
            double query[3] = {point[0], point[1], point[2]};
            double parametricCoordinates[3] = {};
            int subId = 0;
            const vtkIdType cellId = geometry.cellLocator->FindCell(
                query, 0.0, locatorCell, subId, parametricCoordinates,
                locatorWeights.data());
            if (cellId >= 0) {
                locatorElement = static_cast<int>(
                    static_cast<std::size_t>(cellId) / geometry.cellsPerElement);
                if (locatorElement != hint && locatorElement >= 0 &&
                    static_cast<std::size_t>(locatorElement) <
                        geometry.metadata.elementCount &&
                    tryElement(static_cast<std::size_t>(locatorElement))) {
                    return true;
                }
            }
        }
        thread_local std::vector<std::size_t> bvhCandidates;
        geometry.CollectBvhCandidates(point, &bvhCandidates);
        for (const std::size_t element : bvhCandidates) {
            if (static_cast<int>(element) == hint ||
                static_cast<int>(element) == locatorElement) {
                continue;
            }
            if (tryElement(element)) {
                return true;
            }
        }
        if (elementHint != nullptr) {
            *elementHint = -1;
        }
        return false;
    }
};

Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    vtkDataSet* dataSet,
    vtkDataArray* pointArray,
    const Nek5000PointEvaluationOptions& options)
    : impl_(std::make_shared<Impl>(dataSet, pointArray, options)) {}

Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    std::shared_ptr<const Impl> impl)
    : impl_(std::move(impl)) {}

Nek5000SpectralPointEvaluator::~Nek5000SpectralPointEvaluator() = default;
Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    const Nek5000SpectralPointEvaluator&) noexcept = default;
Nek5000SpectralPointEvaluator& Nek5000SpectralPointEvaluator::operator=(
    const Nek5000SpectralPointEvaluator&) noexcept = default;
Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    Nek5000SpectralPointEvaluator&&) noexcept = default;
Nek5000SpectralPointEvaluator& Nek5000SpectralPointEvaluator::operator=(
    Nek5000SpectralPointEvaluator&&) noexcept = default;

int Nek5000SpectralPointEvaluator::numberOfComponents() const noexcept {
    return impl_ == nullptr ? 0 : impl_->components;
}

std::array<double, 6> Nek5000SpectralPointEvaluator::bounds() const noexcept {
    if (impl_ == nullptr) {
        return {};
    }
    const Impl& geometry = impl_->Geometry();
    if (geometry.bvhRoot == std::numeric_limits<std::size_t>::max() ||
        geometry.bvhRoot >= geometry.bvhNodes.size()) {
        return {};
    }
    return geometry.bvhNodes[geometry.bvhRoot].bounds;
}

Nek5000SpectralPointEvaluator
Nek5000SpectralPointEvaluator::WithPointArray(vtkDataArray* pointArray) const {
    if (impl_ == nullptr) {
        throw std::runtime_error(
            "Cannot rebind an empty Nek5000 spectral point evaluator.");
    }
    std::shared_ptr<const Impl> geometry =
        impl_->geometryOwner == nullptr ? impl_ : impl_->geometryOwner;
    return Nek5000SpectralPointEvaluator(
        std::make_shared<Impl>(std::move(geometry), pointArray));
}

bool Nek5000SpectralPointEvaluator::Evaluate(const double physicalPoint[3],
                                             double* tuple,
                                             int* elementHint) const {
    return impl_ != nullptr && impl_->Evaluate(physicalPoint, tuple, elementHint);
}

std::string ToString(DataSetFrameKind kind) {
    switch (kind) {
    case DataSetFrameKind::File:
        return "file";
    case DataSetFrameKind::Nek5000TimeStep:
        return "nek5000_time_step";
    case DataSetFrameKind::EnSightGoldTimeStep:
        return "ensight_gold_time_step";
    case DataSetFrameKind::OpenFoamTimeStep:
        return "openfoam_time_step";
    }
    throw std::runtime_error("Unknown dataset frame kind.");
}

bool IsSupportedDataSource(const std::filesystem::path& sourcePath) {
    if (std::filesystem::is_directory(sourcePath)) {
        if (ResolveOpenFoamCase(sourcePath).has_value()) {
            return true;
        }
        for (const auto& entry :
             std::filesystem::directory_iterator(sourcePath)) {
            if (!entry.is_regular_file()) {
                continue;
            }
            const std::string extension = LowerExtension(entry.path());
            if (IsOrdinaryVtkExtension(extension) || extension == ".pvd" ||
                extension == ".nek5000" || extension == ".nek" ||
                extension == ".case") {
                return true;
            }
        }
        return false;
    }
    if (!std::filesystem::is_regular_file(
            streamcenterplus::vtkhdf::RealPath(sourcePath))) {
        return false;
    }
    const std::string extension = LowerExtension(sourcePath);
    return IsOrdinaryVtkExtension(extension) || extension == ".pvd" ||
           (extension == ".nek5000" || extension == ".nek") ||
           extension == ".case" ||
           extension == ".foam";
}

std::vector<DataSetFrame> EnumerateDataSetFrames(
    const std::filesystem::path& sourcePath,
    double fallbackPhysicalDt) {
    if (!std::isfinite(fallbackPhysicalDt) || fallbackPhysicalDt <= 0.0) {
        throw std::runtime_error("Fallback physical timestep must be positive.");
    }
    if (const auto foam = ResolveOpenFoamCase(sourcePath)) {
        vtkSmartPointer<vtkPOpenFOAMReader> reader = CreateOpenFoamReader(*foam, true);
        const std::vector<double> times = OpenFoamTimeSteps(reader);
        std::vector<DataSetFrame> frames;
        frames.reserve(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            DataSetFrame frame;
            frame.path = *foam;
            frame.time = times[i];
            frame.kind = DataSetFrameKind::OpenFoamTimeStep;
            frame.readerTime = times[i];
            frame.sourceIndex = static_cast<int>(i);
            frames.push_back(std::move(frame));
        }
        return frames;
    }
    if (std::filesystem::is_directory(sourcePath)) {
        std::vector<DataSetFrame> frames =
            EnumerateFolder(sourcePath, fallbackPhysicalDt);
        if (frames.empty()) {
            throw std::runtime_error("No supported dataset frames were found in: " +
                                     sourcePath.string());
        }
        return frames;
    }
    const std::filesystem::path realPath =
        streamcenterplus::vtkhdf::RealPath(sourcePath);
    if (!std::filesystem::is_regular_file(realPath)) {
        throw std::runtime_error("Dataset source does not exist: " +
                                 sourcePath.string());
    }
    const std::string extension = LowerExtension(sourcePath);
    if (extension == ".pvd") {
        return ParsePvd(sourcePath, fallbackPhysicalDt);
    }
    if (extension == ".nek5000" || extension == ".nek") {
        const NekDescriptor descriptor = ParseNekDescriptor(sourcePath);
        std::vector<DataSetFrame> frames;
        frames.reserve(static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        std::vector<bool> meshBearingFrames;
        meshBearingFrames.reserve(
            static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        for (int i = 0; i < descriptor.numberOfTimeSteps; ++i) {
            const NekFieldMeta meta =
                ReadNekHeader(FormatNekFieldPath(descriptor, i, 0));
            DataSetFrame frame;
            frame.path = descriptor.path;
            frame.time = meta.time;
            frame.kind = DataSetFrameKind::Nek5000TimeStep;
            frame.readerTime = meta.time;
            frame.sourceIndex = i;
            frames.push_back(std::move(frame));
            meshBearingFrames.push_back(meta.code.find('X') != std::string::npos);
        }

        int previousMeshFrame = -1;
        for (int i = 0; i < descriptor.numberOfTimeSteps; ++i) {
            if (meshBearingFrames[static_cast<std::size_t>(i)]) {
                previousMeshFrame = i;
            }
            frames[static_cast<std::size_t>(i)].coordinateSourceIndex =
                previousMeshFrame;
        }
        int nextMeshFrame = -1;
        for (int i = descriptor.numberOfTimeSteps - 1; i >= 0; --i) {
            if (meshBearingFrames[static_cast<std::size_t>(i)]) {
                nextMeshFrame = i;
            }
            DataSetFrame& frame = frames[static_cast<std::size_t>(i)];
            if (frame.coordinateSourceIndex < 0) {
                frame.coordinateSourceIndex = nextMeshFrame;
            }
            if (frame.coordinateSourceIndex >= 0) {
                frame.coordinatePath = FormatNekFieldPath(
                    descriptor, frame.coordinateSourceIndex, 0);
            }
        }
        if (nextMeshFrame < 0) {
            throw std::runtime_error(
                "Nek5000 sequence does not contain a mesh-bearing field frame.");
        }
        return frames;
    }
    if (extension == ".case") {
        vtkSmartPointer<vtkGenericEnSightReader> reader =
            CreateEnSightGoldReader(sourcePath, true);
        const std::vector<double> times = EnSightGoldTimeSteps(reader);
        std::vector<DataSetFrame> frames;
        frames.reserve(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            DataSetFrame frame;
            frame.path = std::filesystem::absolute(sourcePath);
            frame.time = times[i];
            frame.kind = DataSetFrameKind::EnSightGoldTimeStep;
            frame.readerTime = times[i];
            frame.sourceIndex = static_cast<int>(i);
            frames.push_back(std::move(frame));
        }
        return frames;
    }
    if (!IsOrdinaryVtkExtension(extension)) {
        throw std::runtime_error("Unsupported dataset source: " + sourcePath.string());
    }
    const std::vector<std::filesystem::path> expanded =
        streamcenterplus::vtkhdf::ExpandTimeSteps(sourcePath);
    std::vector<DataSetFrame> frames;
    frames.reserve(expanded.size());
    for (std::size_t i = 0; i < expanded.size(); ++i) {
        DataSetFrame frame;
        frame.path = expanded[i];
        frame.time = static_cast<double>(i) * fallbackPhysicalDt;
        frame.kind = DataSetFrameKind::File;
        frame.sourceIndex = static_cast<int>(i);
        frames.push_back(std::move(frame));
    }
    return frames;
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const std::filesystem::path& sourcePath,
                                        const ReadOptions& options) {
    std::vector<DataSetFrame> frames = EnumerateDataSetFrames(sourcePath);
    if (frames.empty()) {
        throw std::runtime_error("Dataset source contains no frames: " +
                                 sourcePath.string());
    }
    return ReadDataSet(frames.front(), options);
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    switch (frame.kind) {
    case DataSetFrameKind::File:
        return ReadOrdinaryVtk(frame.path, options);
    case DataSetFrameKind::Nek5000TimeStep:
        return ReadNek5000(frame, options);
    case DataSetFrameKind::EnSightGoldTimeStep:
        return ReadEnSightGold(frame, options);
    case DataSetFrameKind::OpenFoamTimeStep:
        return ReadOpenFoam(frame, options);
    }
    throw std::runtime_error("Unsupported dataset frame kind.");
}

bool HasNek5000SpectralMetadata(vtkDataSet* dataSet) {
    if (dataSet == nullptr || dataSet->GetFieldData() == nullptr) {
        return false;
    }
    vtkDataArray* dimensions =
        dataSet->GetFieldData()->GetArray(kNek5000GllDimensionsField);
    vtkDataArray* count =
        dataSet->GetFieldData()->GetArray(kNek5000ElementCountField);
    return dimensions != nullptr && dimensions->GetNumberOfComponents() == 3 &&
           dimensions->GetNumberOfTuples() == 1 && count != nullptr &&
           count->GetNumberOfTuples() == 1;
}

Nek5000SpectralMetadata GetNek5000SpectralMetadata(vtkDataSet* dataSet) {
    if (!HasNek5000SpectralMetadata(dataSet)) {
        throw std::runtime_error("Dataset has no Streamcenter+ Nek5000 spectral metadata.");
    }
    Nek5000SpectralMetadata metadata;
    vtkDataArray* dimensions =
        dataSet->GetFieldData()->GetArray(kNek5000GllDimensionsField);
    metadata.nx = static_cast<int>(dimensions->GetComponent(0, 0));
    metadata.ny = static_cast<int>(dimensions->GetComponent(0, 1));
    metadata.nz = static_cast<int>(dimensions->GetComponent(0, 2));
    metadata.elementCount = static_cast<std::size_t>(
        dataSet->GetFieldData()->GetArray(kNek5000ElementCountField)
            ->GetComponent(0, 0));
    if (vtkDataArray* ids =
            dataSet->GetFieldData()->GetArray(kNek5000ElementIdsField)) {
        metadata.elementIds.resize(static_cast<std::size_t>(ids->GetNumberOfTuples()));
        for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); ++i) {
            metadata.elementIds[static_cast<std::size_t>(i)] =
                static_cast<int>(ids->GetComponent(i, 0));
        }
    }
    if (metadata.nx < 2 || metadata.ny < 2 || metadata.nz < 1 ||
        metadata.elementCount == 0 || metadata.pointCount() !=
                                            static_cast<std::size_t>(dataSet->GetNumberOfPoints())) {
        throw std::runtime_error("Dataset has inconsistent Nek5000 spectral metadata.");
    }
    return metadata;
}

bool HasSameNek5000SpectralGeometry(vtkDataSet* lhs,
                                    vtkDataSet* rhs,
                                    double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        throw std::runtime_error(
            "Nek5000 spectral geometry tolerance must be finite and non-negative.");
    }
    if (lhs == nullptr || rhs == nullptr ||
        !HasNek5000SpectralMetadata(lhs) ||
        !HasNek5000SpectralMetadata(rhs)) {
        return false;
    }
    const Nek5000SpectralMetadata lhsMetadata =
        GetNek5000SpectralMetadata(lhs);
    const Nek5000SpectralMetadata rhsMetadata =
        GetNek5000SpectralMetadata(rhs);
    if (lhsMetadata.nx != rhsMetadata.nx ||
        lhsMetadata.ny != rhsMetadata.ny ||
        lhsMetadata.nz != rhsMetadata.nz ||
        lhsMetadata.elementCount != rhsMetadata.elementCount ||
        lhsMetadata.elementIds != rhsMetadata.elementIds ||
        lhs->GetNumberOfPoints() != rhs->GetNumberOfPoints()) {
        return false;
    }

    vtkPointSet* lhsPointSet = vtkPointSet::SafeDownCast(lhs);
    vtkPointSet* rhsPointSet = vtkPointSet::SafeDownCast(rhs);
    vtkDoubleArray* lhsCoordinates =
        lhsPointSet == nullptr || lhsPointSet->GetPoints() == nullptr
            ? nullptr
            : vtkDoubleArray::SafeDownCast(lhsPointSet->GetPoints()->GetData());
    vtkDoubleArray* rhsCoordinates =
        rhsPointSet == nullptr || rhsPointSet->GetPoints() == nullptr
            ? nullptr
            : vtkDoubleArray::SafeDownCast(rhsPointSet->GetPoints()->GetData());
    if (lhsCoordinates != nullptr && rhsCoordinates != nullptr &&
        lhsCoordinates->GetNumberOfComponents() == 3 &&
        rhsCoordinates->GetNumberOfComponents() == 3) {
        const std::size_t valueCount =
            3 * static_cast<std::size_t>(lhs->GetNumberOfPoints());
        const double* lhsValues = lhsCoordinates->GetPointer(0);
        const double* rhsValues = rhsCoordinates->GetPointer(0);
        for (std::size_t value = 0; value < valueCount; ++value) {
            if (!std::isfinite(lhsValues[value]) ||
                !std::isfinite(rhsValues[value]) ||
                (tolerance == 0.0
                     ? lhsValues[value] != rhsValues[value]
                     : std::abs(lhsValues[value] - rhsValues[value]) >
                           tolerance)) {
                return false;
            }
        }
        return true;
    }

    double lhsPoint[3] = {};
    double rhsPoint[3] = {};
    for (vtkIdType pointId = 0; pointId < lhs->GetNumberOfPoints(); ++pointId) {
        lhs->GetPoint(pointId, lhsPoint);
        rhs->GetPoint(pointId, rhsPoint);
        for (int axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(lhsPoint[axis]) ||
                !std::isfinite(rhsPoint[axis])) {
                return false;
            }
            if (tolerance == 0.0) {
                if (lhsPoint[axis] != rhsPoint[axis]) {
                    return false;
                }
            } else if (std::abs(lhsPoint[axis] - rhsPoint[axis]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

vtkSmartPointer<vtkDoubleArray> ComputeNek5000SpectralGradient(
    vtkDataSet* dataSet,
    vtkDataArray* pointArray,
    const std::string& outputName) {
    if (dataSet == nullptr || pointArray == nullptr ||
        pointArray->GetNumberOfTuples() != dataSet->GetNumberOfPoints()) {
        throw std::runtime_error(
            "Nek5000 spectral gradient requires a complete point-data array.");
    }
    const Nek5000SpectralMetadata meta = GetNek5000SpectralMetadata(dataSet);
    const int components = pointArray->GetNumberOfComponents();
    if (components <= 0) {
        throw std::runtime_error("Nek5000 spectral gradient input has no components.");
    }
    const std::vector<double> nodesX = GllNodes(meta.nx);
    const std::vector<double> nodesY = GllNodes(meta.ny);
    const std::vector<double> nodesZ = meta.is3D() ? GllNodes(meta.nz)
                                                   : std::vector<double>{0.0};
    const std::vector<double> dx = DerivativeMatrix(nodesX);
    const std::vector<double> dy = DerivativeMatrix(nodesY);
    const std::vector<double> dz = meta.is3D() ? DerivativeMatrix(nodesZ)
                                               : std::vector<double>();

    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(dataSet);
    vtkDataArray* coordinates =
        pointSet == nullptr || pointSet->GetPoints() == nullptr
            ? nullptr
            : pointSet->GetPoints()->GetData();
    if (coordinates == nullptr || coordinates->GetNumberOfComponents() != 3) {
        throw std::runtime_error(
            "Nek5000 geometry must expose a three-component point-coordinate array.");
    }
    vtkDoubleArray* doubleCoordinates = vtkDoubleArray::SafeDownCast(coordinates);
    const double* coordinateValues =
        doubleCoordinates == nullptr ? nullptr : doubleCoordinates->GetPointer(0);
    const auto coordinate = [coordinates, coordinateValues](std::size_t pointIndex,
                                                            int component) {
        return coordinateValues == nullptr
                   ? coordinates->GetComponent(static_cast<vtkIdType>(pointIndex), component)
                   : coordinateValues[3 * pointIndex + static_cast<std::size_t>(component)];
    };

    vtkSmartPointer<vtkDoubleArray> output =
        vtkSmartPointer<vtkDoubleArray>::New();
    const std::string name = outputName.empty()
                                 ? std::string(pointArray->GetName() == nullptr
                                                   ? "SpectralGradient"
                                                   : pointArray->GetName()) +
                                       "_SpectralGradient"
                                 : outputName;
    output->SetName(name.c_str());
    output->SetNumberOfComponents(3 * components);
    output->SetNumberOfTuples(static_cast<vtkIdType>(meta.pointCount()));
    output->FillValue(0.0);

    const std::size_t nxyz = meta.pointsPerElement();
    std::vector<double> qr(static_cast<std::size_t>(components), 0.0);
    std::vector<double> qs(static_cast<std::size_t>(components), 0.0);
    std::vector<double> qt(static_cast<std::size_t>(components), 0.0);
    for (std::size_t element = 0; element < meta.elementCount; ++element) {
        const std::size_t base = element * nxyz;
        for (int k = 0; k < meta.nz; ++k) {
            for (int j = 0; j < meta.ny; ++j) {
                for (int i = 0; i < meta.nx; ++i) {
                    std::array<double, 3> xr{};
                    std::array<double, 3> xs{};
                    std::array<double, 3> xt{};
                    std::fill(qr.begin(), qr.end(), 0.0);
                    std::fill(qs.begin(), qs.end(), 0.0);
                    std::fill(qt.begin(), qt.end(), 0.0);
                    for (int m = 0; m < meta.nx; ++m) {
                        const double d = dx[static_cast<std::size_t>(i) *
                                                static_cast<std::size_t>(meta.nx) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, m, j, k);
                        for (int axis = 0; axis < 3; ++axis) {
                            xr[static_cast<std::size_t>(axis)] +=
                                d * coordinate(index, axis);
                        }
                        for (int component = 0; component < components; ++component) {
                            qr[static_cast<std::size_t>(component)] +=
                                d * pointArray->GetComponent(static_cast<vtkIdType>(index),
                                                             component);
                        }
                    }
                    for (int m = 0; m < meta.ny; ++m) {
                        const double d = dy[static_cast<std::size_t>(j) *
                                                static_cast<std::size_t>(meta.ny) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, i, m, k);
                        for (int axis = 0; axis < 3; ++axis) {
                            xs[static_cast<std::size_t>(axis)] +=
                                d * coordinate(index, axis);
                        }
                        for (int component = 0; component < components; ++component) {
                            qs[static_cast<std::size_t>(component)] +=
                                d * pointArray->GetComponent(static_cast<vtkIdType>(index),
                                                             component);
                        }
                    }
                    if (meta.is3D()) {
                        for (int m = 0; m < meta.nz; ++m) {
                            const double d = dz[static_cast<std::size_t>(k) *
                                                    static_cast<std::size_t>(meta.nz) +
                                                static_cast<std::size_t>(m)];
                            const std::size_t index =
                                base + NekLocalIndex(meta.nx, meta.ny, i, j, m);
                            for (int axis = 0; axis < 3; ++axis) {
                                xt[static_cast<std::size_t>(axis)] +=
                                    d * coordinate(index, axis);
                            }
                            for (int component = 0; component < components; ++component) {
                                qt[static_cast<std::size_t>(component)] +=
                                    d * pointArray->GetComponent(
                                            static_cast<vtkIdType>(index), component);
                            }
                        }
                    }

                    const std::size_t outputIndex =
                        base + NekLocalIndex(meta.nx, meta.ny, i, j, k);
                    if (meta.is3D()) {
                        const double determinant =
                            xr[0] * (xs[1] * xt[2] - xt[1] * xs[2]) -
                            xs[0] * (xr[1] * xt[2] - xt[1] * xr[2]) +
                            xt[0] * (xr[1] * xs[2] - xs[1] * xr[2]);
                        if (!std::isfinite(determinant) ||
                            std::abs(determinant) <= 1.0e-30) {
                            throw std::runtime_error(
                                "Degenerate Nek5000 Jacobian while computing a spectral gradient.");
                        }
                        const double rX = (xs[1] * xt[2] - xt[1] * xs[2]) / determinant;
                        const double sX = (xt[1] * xr[2] - xr[1] * xt[2]) / determinant;
                        const double tX = (xr[1] * xs[2] - xs[1] * xr[2]) / determinant;
                        const double rY = (xt[0] * xs[2] - xs[0] * xt[2]) / determinant;
                        const double sY = (xr[0] * xt[2] - xt[0] * xr[2]) / determinant;
                        const double tY = (xs[0] * xr[2] - xr[0] * xs[2]) / determinant;
                        const double rZ = (xs[0] * xt[1] - xt[0] * xs[1]) / determinant;
                        const double sZ = (xt[0] * xr[1] - xr[0] * xt[1]) / determinant;
                        const double tZ = (xr[0] * xs[1] - xs[0] * xr[1]) / determinant;
                        for (int component = 0; component < components; ++component) {
                            const std::size_t c = static_cast<std::size_t>(component);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 0,
                                                 qr[c] * rX + qs[c] * sX + qt[c] * tX);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 1,
                                                 qr[c] * rY + qs[c] * sY + qt[c] * tY);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 2,
                                                 qr[c] * rZ + qs[c] * sZ + qt[c] * tZ);
                        }
                    } else {
                        const double g00 = xr[0] * xr[0] + xr[1] * xr[1] + xr[2] * xr[2];
                        const double g01 = xr[0] * xs[0] + xr[1] * xs[1] + xr[2] * xs[2];
                        const double g11 = xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2];
                        const double determinant = g00 * g11 - g01 * g01;
                        if (!std::isfinite(determinant) || determinant <= 1.0e-30) {
                            throw std::runtime_error(
                                "Degenerate 2-D Nek5000 metric while computing a spectral gradient.");
                        }
                        for (int component = 0; component < components; ++component) {
                            const std::size_t c = static_cast<std::size_t>(component);
                            const double coefficientR = (g11 * qr[c] - g01 * qs[c]) / determinant;
                            const double coefficientS = (-g01 * qr[c] + g00 * qs[c]) / determinant;
                            for (int axis = 0; axis < 3; ++axis) {
                                output->SetComponent(
                                    static_cast<vtkIdType>(outputIndex),
                                    3 * component + axis,
                                    coefficientR * xr[static_cast<std::size_t>(axis)] +
                                        coefficientS * xs[static_cast<std::size_t>(axis)]);
                            }
                        }
                    }
                }
            }
        }
    }
    return output;
}

}  // namespace streamcenterplus::io
