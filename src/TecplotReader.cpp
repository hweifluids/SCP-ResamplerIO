#include "TecplotReader.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersionMacros.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace streamcenterplus::io::detail {

namespace {

namespace fs = std::filesystem;

constexpr float kTecplotEndOfHeaderMarker = 357.0F;
constexpr float kTecplotZoneHeaderV112Marker = 299.0F;
constexpr float kTecplotZoneHeaderV191Marker = 298.0F;
constexpr float kTecplotUserRecordMarker = 699.0F;
constexpr float kTecplotDataSetAuxMarker = 799.0F;
constexpr float kTecplotVariableAuxMarker = 899.0F;

constexpr int kTecplotFileTypeFull = 0;
constexpr int kTecplotFileTypeGrid = 1;
constexpr int kTecplotFileTypeSolution = 2;

constexpr int kTecplotZoneOrdered = 0;
constexpr int kTecplotZoneLine = 1;
constexpr int kTecplotZoneTriangle = 2;
constexpr int kTecplotZoneQuadrilateral = 3;
constexpr int kTecplotZoneTetrahedron = 4;
constexpr int kTecplotZoneBrick = 5;
constexpr int kTecplotZonePolygon = 6;
constexpr int kTecplotZonePolyhedron = 7;

bool MarkerMatches(float actual, float expected) {
    return std::abs(actual - expected) <= 1.0e-4F;
}

bool HostIsLittleEndian() {
    const std::uint16_t value = 1;
    return *reinterpret_cast<const unsigned char*>(&value) == 1;
}

template <typename T>
T ByteSwap(T value) {
    std::array<unsigned char, sizeof(T)> bytes{};
    std::memcpy(bytes.data(), &value, sizeof(T));
    std::reverse(bytes.begin(), bytes.end());
    std::memcpy(&value, bytes.data(), sizeof(T));
    return value;
}

std::string PathText(const fs::path& path) {
    return path.string();
}

class BinaryInput {
public:
    explicit BinaryInput(const fs::path& path)
        : path_(path), input_(path, std::ios::binary) {
        if (!input_) {
            Fail("Cannot open Tecplot binary file");
        }

        std::array<char, 8> magic{};
        ReadRaw(magic.data(), magic.size(), "Tecplot binary magic");
        const std::string identifier(magic.data(), magic.size());
        if (identifier == "#!TDV112") {
            version_ = 112;
        } else if (identifier == "#!TDV191") {
            version_ = 191;
        } else {
            Fail("Unsupported classic Tecplot binary signature");
        }

        std::int32_t byteOrder = 0;
        ReadRaw(&byteOrder, sizeof(byteOrder), "Tecplot byte-order marker");
        if (byteOrder == 1) {
            swapBytes_ = false;
        } else if (ByteSwap(byteOrder) == 1) {
            swapBytes_ = true;
        } else {
            Fail("Tecplot byte-order marker is invalid");
        }
    }

    int version() const noexcept {
        return version_;
    }

    std::int32_t ReadInt32(const std::string& what) {
        return Read<std::int32_t>(what);
    }

    std::int64_t ReadInt64(const std::string& what) {
        return Read<std::int64_t>(what);
    }

    float ReadFloat(const std::string& what) {
        return Read<float>(what);
    }

    double ReadDouble(const std::string& what) {
        return Read<double>(what);
    }

    std::string ReadString(const std::string& what) {
        constexpr std::size_t kMaximumStringLength = 1u << 20;
        std::string value;
        value.reserve(64);
        while (true) {
            const std::int32_t code = ReadInt32(what);
            if (code == 0) {
                return value;
            }
            if (code < 0 || code > std::numeric_limits<unsigned char>::max()) {
                Fail(what + " contains a non-byte character code");
            }
            if (value.size() >= kMaximumStringLength) {
                Fail(what + " exceeds the supported maximum length");
            }
            value.push_back(static_cast<char>(code));
        }
    }

    vtkSmartPointer<vtkDoubleArray> ReadNumericArray(
        int dataFormat,
        std::size_t valueCount,
        const std::string& name,
        const std::string& what) {
        auto values = vtkSmartPointer<vtkDoubleArray>::New();
        values->SetName(name.c_str());
        values->SetNumberOfComponents(1);
        values->SetNumberOfValues(CheckedVtkIdType(valueCount, what));
        double* destination = values->GetPointer(0);
        switch (dataFormat) {
        case 1:
            ReadConverted<float>(destination, valueCount, what);
            break;
        case 2:
            ReadConverted<double>(destination, valueCount, what);
            break;
        case 3:
            ReadConverted<std::int32_t>(destination, valueCount, what);
            break;
        case 4:
            ReadConverted<std::int16_t>(destination, valueCount, what);
            break;
        case 5:
            ReadConverted<std::uint8_t>(destination, valueCount, what);
            break;
        case 6:
            Fail(what + " uses Tecplot BIT data, which is not supported for flow-field input");
        default:
            Fail(what + " has an unknown Tecplot variable data format");
        }
        return values;
    }

    void SkipInt32Values(std::size_t valueCount, const std::string& what) {
        constexpr std::size_t kChunkSize = 1u << 15;
        std::vector<std::int32_t> buffer;
        std::size_t readCount = 0;
        while (readCount < valueCount) {
            const std::size_t current = std::min(kChunkSize, valueCount - readCount);
            buffer.resize(current);
            ReadRaw(buffer.data(), current * sizeof(std::int32_t), what);
            readCount += current;
        }
    }

    [[noreturn]] void Fail(const std::string& message) const {
        throw std::runtime_error(message + ": " + PathText(path_));
    }

private:
    template <typename T>
    T Read(const std::string& what) {
        T value{};
        ReadRaw(&value, sizeof(value), what);
        return swapBytes_ ? ByteSwap(value) : value;
    }

    template <typename T>
    void ReadConverted(double* destination,
                       std::size_t valueCount,
                       const std::string& what) {
        constexpr std::size_t kChunkSize = 1u << 15;
        std::vector<T> buffer;
        std::size_t readCount = 0;
        while (readCount < valueCount) {
            const std::size_t current = std::min(kChunkSize, valueCount - readCount);
            buffer.resize(current);
            ReadRaw(buffer.data(), current * sizeof(T), what);
            for (std::size_t index = 0; index < current; ++index) {
                T value = buffer[index];
                if (swapBytes_) {
                    value = ByteSwap(value);
                }
                destination[readCount + index] = static_cast<double>(value);
            }
            readCount += current;
        }
    }

    static vtkIdType CheckedVtkIdType(std::size_t value,
                                      const std::string& what) {
        if (value > static_cast<std::size_t>(
                        std::numeric_limits<vtkIdType>::max())) {
            throw std::runtime_error(what + " exceeds VTK's id range");
        }
        return static_cast<vtkIdType>(value);
    }

    void ReadRaw(void* destination,
                 std::size_t byteCount,
                 const std::string& what) {
        if (byteCount > static_cast<std::size_t>(
                            std::numeric_limits<std::streamsize>::max())) {
            Fail(what + " is too large to read");
        }
        input_.read(static_cast<char*>(destination),
                    static_cast<std::streamsize>(byteCount));
        if (!input_) {
            Fail("Cannot read " + what);
        }
    }

    fs::path path_;
    std::ifstream input_;
    bool swapBytes_ = false;
    int version_ = 0;
};

std::size_t PositiveCount(std::int64_t value,
                          const std::string& what,
                          BinaryInput* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot count reader error.");
    }
    if (value <= 0) {
        input->Fail(what + " must be positive");
    }
    const auto unsignedValue = static_cast<std::uint64_t>(value);
    if (unsignedValue > static_cast<std::uint64_t>(
                            std::numeric_limits<std::size_t>::max()) ||
        unsignedValue > static_cast<std::uint64_t>(
                            std::numeric_limits<vtkIdType>::max())) {
        input->Fail(what + " exceeds the supported range");
    }
    return static_cast<std::size_t>(unsignedValue);
}

std::size_t NonNegativeCount(std::int64_t value,
                             const std::string& what,
                             BinaryInput* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot count reader error.");
    }
    if (value < 0) {
        input->Fail(what + " must not be negative");
    }
    const auto unsignedValue = static_cast<std::uint64_t>(value);
    if (unsignedValue > static_cast<std::uint64_t>(
                            std::numeric_limits<std::size_t>::max()) ||
        unsignedValue > static_cast<std::uint64_t>(
                            std::numeric_limits<vtkIdType>::max())) {
        input->Fail(what + " exceeds the supported range");
    }
    return static_cast<std::size_t>(unsignedValue);
}

std::string NormalizeVariableName(std::string value) {
    std::string normalized;
    normalized.reserve(value.size());
    for (const unsigned char character : value) {
        if (std::isalnum(character) != 0) {
            normalized.push_back(static_cast<char>(std::tolower(character)));
        }
    }
    return normalized;
}

bool IsCoordinateName(const std::string& normalized, char axis) {
    const std::string oneCharacter(1, axis);
    const std::string coordinate = std::string("coordinate") + axis;
    const std::string coordinateReverse = oneCharacter + "coordinate";
    const std::string abbreviated = std::string("coord") + axis;
    return normalized == oneCharacter || normalized == coordinate ||
           normalized == coordinateReverse || normalized == abbreviated;
}

struct CoordinateVariables {
    int x = -1;
    int y = -1;
    int z = -1;
};

CoordinateVariables FindCoordinateVariables(
    const std::vector<std::string>& variableNames) {
    CoordinateVariables result;
    for (std::size_t index = 0; index < variableNames.size(); ++index) {
        const std::string normalized = NormalizeVariableName(variableNames[index]);
        if (result.x < 0 && IsCoordinateName(normalized, 'x')) {
            result.x = static_cast<int>(index);
        } else if (result.y < 0 && IsCoordinateName(normalized, 'y')) {
            result.y = static_cast<int>(index);
        } else if (result.z < 0 && IsCoordinateName(normalized, 'z')) {
            result.z = static_cast<int>(index);
        }
    }
    return result;
}

std::string VariableNameList(const std::vector<std::string>& names) {
    std::ostringstream output;
    for (std::size_t index = 0; index < names.size(); ++index) {
        if (index != 0) {
            output << ", ";
        }
        output << '"' << names[index] << '"';
    }
    return output.str();
}

bool IsFiniteElementZone(int zoneType) {
    return zoneType >= kTecplotZoneLine && zoneType <= kTecplotZoneBrick;
}

bool IsFaceBasedZone(int zoneType) {
    return zoneType == kTecplotZonePolygon ||
           zoneType == kTecplotZonePolyhedron;
}

int NodesPerCell(int zoneType) {
    switch (zoneType) {
    case kTecplotZoneLine:
        return 2;
    case kTecplotZoneTriangle:
        return 3;
    case kTecplotZoneQuadrilateral:
    case kTecplotZoneTetrahedron:
        return 4;
    case kTecplotZoneBrick:
        return 8;
    default:
        return 0;
    }
}

int VtkCellType(int zoneType) {
    switch (zoneType) {
    case kTecplotZoneLine:
        return VTK_LINE;
    case kTecplotZoneTriangle:
        return VTK_TRIANGLE;
    case kTecplotZoneQuadrilateral:
        return VTK_QUAD;
    case kTecplotZoneTetrahedron:
        return VTK_TETRA;
    case kTecplotZoneBrick:
        return VTK_HEXAHEDRON;
    default:
        return VTK_EMPTY_CELL;
    }
}

const char* ZoneTypeText(int zoneType) {
    switch (zoneType) {
    case kTecplotZoneOrdered:
        return "ORDERED";
    case kTecplotZoneLine:
        return "FELINESEG";
    case kTecplotZoneTriangle:
        return "FETRIANGLE";
    case kTecplotZoneQuadrilateral:
        return "FEQUADRILATERAL";
    case kTecplotZoneTetrahedron:
        return "FETETRAHEDRON";
    case kTecplotZoneBrick:
        return "FEBRICK";
    case kTecplotZonePolygon:
        return "FEPOLYGON";
    case kTecplotZonePolyhedron:
        return "FEPOLYHEDRON";
    default:
        return "unknown";
    }
}

int FacesPerCell(int zoneType) {
    switch (zoneType) {
    case kTecplotZoneLine:
        return 2;
    case kTecplotZoneTriangle:
        return 3;
    case kTecplotZoneQuadrilateral:
        return 4;
    case kTecplotZoneTetrahedron:
        return 4;
    case kTecplotZoneBrick:
        return 6;
    default:
        return 0;
    }
}

struct TecplotZoneHeader {
    int headerVersion = 0;
    std::string name;
    int zoneType = -1;
    std::vector<int> variableLocations;
    bool rawFaceNeighbors = false;
    int miscellaneousFaceNeighborCount = 0;
    int miscellaneousFaceNeighborMode = -1;
    std::array<int, 3> orderedDimensions = {{0, 0, 0}};
    std::size_t pointCount = 0;
    std::size_t cellCount = 0;
    std::size_t faceCount = 0;
    std::size_t faceNodeCount = 0;
    std::size_t boundaryFaceCount = 0;
    std::size_t boundaryConnectionCount = 0;
};

std::size_t CheckedProduct(const std::array<int, 3>& values,
                           const std::string& what,
                           BinaryInput* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot product reader error.");
    }
    std::size_t product = 1;
    for (const int value : values) {
        if (value <= 0) {
            input->Fail(what + " dimensions must be positive");
        }
        const std::size_t factor = static_cast<std::size_t>(value);
        if (product > std::numeric_limits<std::size_t>::max() / factor ||
            product > static_cast<std::size_t>(
                          std::numeric_limits<vtkIdType>::max()) / factor) {
            input->Fail(what + " exceeds the supported range");
        }
        product *= factor;
    }
    return product;
}

std::size_t OrderedCellCount(const std::array<int, 3>& dimensions,
                             BinaryInput* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot ordered-cell reader error.");
    }
    std::array<int, 3> cellDimensions = {{1, 1, 1}};
    bool hasCells = false;
    for (std::size_t axis = 0; axis < dimensions.size(); ++axis) {
        if (dimensions[axis] <= 0) {
            input->Fail("Tecplot ordered dimensions must be positive");
        }
        if (dimensions[axis] > 1) {
            cellDimensions[axis] = dimensions[axis] - 1;
            hasCells = true;
        }
    }
    return hasCells ? CheckedProduct(cellDimensions, "Tecplot ordered cell count", input)
                    : 0;
}

std::array<int, 3> OrderedCellPayloadDimensions(
    const std::array<int, 3>& dimensions,
    BinaryInput* input) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot ordered-payload reader error.");
    }
    std::array<int, 3> payloadDimensions = dimensions;
    int lastActiveAxis = -1;
    for (std::size_t axis = 0; axis < dimensions.size(); ++axis) {
        if (dimensions[axis] <= 0) {
            input->Fail("Tecplot ordered dimensions must be positive");
        }
        if (dimensions[axis] > 1) {
            lastActiveAxis = static_cast<int>(axis);
        }
    }
    if (lastActiveAxis < 0) {
        return {{0, 0, 0}};
    }
    payloadDimensions[static_cast<std::size_t>(lastActiveAxis)] -= 1;
    return payloadDimensions;
}

std::size_t OrderedCellPayloadCount(const TecplotZoneHeader& header,
                                    BinaryInput* input) {
    const std::array<int, 3> dimensions =
        OrderedCellPayloadDimensions(header.orderedDimensions, input);
    if (dimensions[0] == 0) {
        return 0;
    }
    return CheckedProduct(dimensions, "Tecplot ordered cell payload", input);
}

TecplotZoneHeader ReadZoneHeader(BinaryInput* input,
                                 int markerVersion,
                                 std::size_t variableCount) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot header reader error.");
    }
    TecplotZoneHeader zone;
    zone.headerVersion = markerVersion;
    zone.name = input->ReadString("Tecplot zone name");
    (void)input->ReadInt32("Tecplot zone parent");
    (void)input->ReadInt32("Tecplot zone strand");
    (void)input->ReadDouble("Tecplot zone solution time");
    (void)input->ReadInt32("Tecplot zone color");
    zone.zoneType = input->ReadInt32("Tecplot zone type");
    if (zone.zoneType < kTecplotZoneOrdered ||
        zone.zoneType > kTecplotZonePolyhedron) {
        input->Fail("Tecplot zone type is invalid");
    }

    const int specifiesVariableLocations =
        input->ReadInt32("Tecplot variable-location flag");
    if (specifiesVariableLocations == 0) {
        zone.variableLocations.assign(variableCount, 0);
    } else if (specifiesVariableLocations == 1) {
        zone.variableLocations.reserve(variableCount);
        for (std::size_t index = 0; index < variableCount; ++index) {
            const int location =
                input->ReadInt32("Tecplot variable location");
            if (location != 0 && location != 1) {
                input->Fail("Tecplot variable location is invalid");
            }
            zone.variableLocations.push_back(location);
        }
    } else {
        input->Fail("Tecplot variable-location flag is invalid");
    }

    const int hasRawFaceNeighbors =
        input->ReadInt32("Tecplot raw-face-neighbor flag");
    if (hasRawFaceNeighbors != 0 && hasRawFaceNeighbors != 1) {
        input->Fail("Tecplot raw-face-neighbor flag is invalid");
    }
    if ((zone.zoneType == kTecplotZoneOrdered ||
         IsFaceBasedZone(zone.zoneType)) &&
        hasRawFaceNeighbors != 0) {
        input->Fail(std::string("Tecplot ") + ZoneTypeText(zone.zoneType) +
                    " zones cannot contain raw face neighbors");
    }
    zone.rawFaceNeighbors = hasRawFaceNeighbors != 0;
    zone.miscellaneousFaceNeighborCount =
        input->ReadInt32("Tecplot miscellaneous face-neighbor count");
    if (zone.miscellaneousFaceNeighborCount < 0) {
        input->Fail("Tecplot miscellaneous face-neighbor count is invalid");
    }
    if (IsFaceBasedZone(zone.zoneType) &&
        zone.miscellaneousFaceNeighborCount != 0) {
        input->Fail(std::string("Tecplot ") + ZoneTypeText(zone.zoneType) +
                    " zones cannot contain miscellaneous face neighbors");
    }
    if (zone.miscellaneousFaceNeighborCount > 0) {
        zone.miscellaneousFaceNeighborMode =
            input->ReadInt32("Tecplot miscellaneous face-neighbor mode");
        if (zone.miscellaneousFaceNeighborMode < 0 ||
            zone.miscellaneousFaceNeighborMode > 3) {
            input->Fail("Tecplot miscellaneous face-neighbor mode is invalid");
        }
        if (zone.zoneType != kTecplotZoneOrdered) {
            (void)input->ReadInt32(
                "Tecplot finite-element face-neighbor completion flag");
        }
    }

    if (zone.zoneType == kTecplotZoneOrdered) {
        zone.orderedDimensions[0] =
            input->ReadInt32("Tecplot ordered I dimension");
        zone.orderedDimensions[1] =
            input->ReadInt32("Tecplot ordered J dimension");
        zone.orderedDimensions[2] =
            input->ReadInt32("Tecplot ordered K dimension");
        zone.pointCount = CheckedProduct(
            zone.orderedDimensions, "Tecplot ordered point count", input);
        zone.cellCount = OrderedCellCount(zone.orderedDimensions, input);
    } else {
        zone.pointCount = PositiveCount(
            input->ReadInt32("Tecplot finite-element point count"),
            "Tecplot finite-element point count", input);
        if (IsFaceBasedZone(zone.zoneType)) {
            const std::int64_t rawFaceCount = markerVersion == 191
                ? input->ReadInt64("Tecplot face count")
                : input->ReadInt32("Tecplot face count");
            const std::int64_t rawFaceNodeCount = markerVersion == 191
                ? input->ReadInt64("Tecplot face-node count")
                : input->ReadInt32("Tecplot face-node count");
            zone.faceCount = PositiveCount(rawFaceCount,
                                           "Tecplot face count", input);
            zone.faceNodeCount = PositiveCount(rawFaceNodeCount,
                                               "Tecplot face-node count", input);
            zone.boundaryFaceCount = NonNegativeCount(
                input->ReadInt32("Tecplot boundary-face count"),
                "Tecplot boundary-face count", input);
            zone.boundaryConnectionCount = NonNegativeCount(
                input->ReadInt32("Tecplot boundary-connection count"),
                "Tecplot boundary-connection count", input);
            if (zone.zoneType == kTecplotZonePolygon) {
                if (zone.faceCount >
                    std::numeric_limits<std::size_t>::max() / 2u ||
                    zone.faceNodeCount != zone.faceCount * 2u) {
                    input->Fail("Tecplot FEPOLYGON face-node count is invalid");
                }
            }
        }
        zone.cellCount = PositiveCount(
            input->ReadInt32("Tecplot finite-element cell count"),
            "Tecplot finite-element cell count", input);
        (void)input->ReadInt32("Tecplot finite-element I cell dimension");
        (void)input->ReadInt32("Tecplot finite-element J cell dimension");
        (void)input->ReadInt32("Tecplot finite-element K cell dimension");
    }

    while (true) {
        const int hasAuxiliaryData =
            input->ReadInt32("Tecplot zone auxiliary-data flag");
        if (hasAuxiliaryData == 0) {
            break;
        }
        if (hasAuxiliaryData != 1) {
            input->Fail("Tecplot zone auxiliary-data flag is invalid");
        }
        (void)input->ReadString("Tecplot zone auxiliary-data name");
        const int format =
            input->ReadInt32("Tecplot zone auxiliary-data format");
        if (format != 0) {
            input->Fail("Tecplot zone auxiliary-data format is unsupported");
        }
        (void)input->ReadString("Tecplot zone auxiliary-data value");
    }

    if (zone.zoneType != kTecplotZoneOrdered &&
        !IsFiniteElementZone(zone.zoneType) &&
        !IsFaceBasedZone(zone.zoneType)) {
        input->Fail(std::string("Tecplot binary ") + ZoneTypeText(zone.zoneType) +
                    " zones are not supported");
    }
    return zone;
}

void SkipHeaderRecord(BinaryInput* input, float marker) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot header-record reader error.");
    }
    if (MarkerMatches(marker, kTecplotUserRecordMarker)) {
        (void)input->ReadString("Tecplot user record");
        return;
    }
    if (MarkerMatches(marker, kTecplotDataSetAuxMarker)) {
        (void)input->ReadString("Tecplot dataset auxiliary-data name");
        const int format =
            input->ReadInt32("Tecplot dataset auxiliary-data format");
        if (format != 0) {
            input->Fail("Tecplot dataset auxiliary-data format is unsupported");
        }
        (void)input->ReadString("Tecplot dataset auxiliary-data value");
        return;
    }
    if (MarkerMatches(marker, kTecplotVariableAuxMarker)) {
        (void)input->ReadInt32("Tecplot variable auxiliary-data index");
        (void)input->ReadString("Tecplot variable auxiliary-data name");
        const int format =
            input->ReadInt32("Tecplot variable auxiliary-data format");
        if (format != 0) {
            input->Fail("Tecplot variable auxiliary-data format is unsupported");
        }
        (void)input->ReadString("Tecplot variable auxiliary-data value");
        return;
    }
    input->Fail("Tecplot binary file contains an unsupported non-data header record");
}

struct DecodedZone {
    TecplotZoneHeader header;
    std::vector<vtkSmartPointer<vtkDoubleArray>> variables;
    vtkSmartPointer<vtkDataSet> grid;
};

void RequirePointCoordinate(const TecplotZoneHeader& header,
                            int variableIndex,
                            const char* axis,
                            BinaryInput* input) {
    if (variableIndex < 0) {
        return;
    }
    const std::size_t index = static_cast<std::size_t>(variableIndex);
    if (index >= header.variableLocations.size() ||
        header.variableLocations[index] != 0) {
        input->Fail(std::string("Tecplot ") + axis +
                    " coordinate variable must be point-centered");
    }
}

vtkSmartPointer<vtkDoubleArray> CompactOrderedCellValues(
    const TecplotZoneHeader& header,
    vtkDoubleArray* source,
    BinaryInput* input) {
    if (source == nullptr || input == nullptr) {
        throw std::runtime_error("Internal Tecplot ordered cell-array reader error.");
    }
    const std::size_t payloadCount = OrderedCellPayloadCount(header, input);
    if (source->GetNumberOfValues() != static_cast<vtkIdType>(payloadCount)) {
        input->Fail("Tecplot ordered cell payload has an unexpected size");
    }

    auto values = vtkSmartPointer<vtkDoubleArray>::New();
    values->SetName(source->GetName());
    values->SetNumberOfComponents(1);
    values->SetNumberOfValues(static_cast<vtkIdType>(header.cellCount));
    if (header.cellCount == 0) {
        return values;
    }

    const std::array<int, 3> payloadDimensions =
        OrderedCellPayloadDimensions(header.orderedDimensions, input);
    std::array<int, 3> cellDimensions = {{1, 1, 1}};
    for (std::size_t axis = 0; axis < cellDimensions.size(); ++axis) {
        if (header.orderedDimensions[axis] > 1) {
            cellDimensions[axis] = header.orderedDimensions[axis] - 1;
        }
    }

    vtkIdType targetIndex = 0;
    for (int k = 0; k < cellDimensions[2]; ++k) {
        for (int j = 0; j < cellDimensions[1]; ++j) {
            for (int i = 0; i < cellDimensions[0]; ++i) {
                const std::size_t sourceIndex =
                    static_cast<std::size_t>(i) +
                    static_cast<std::size_t>(payloadDimensions[0]) *
                        (static_cast<std::size_t>(j) +
                         static_cast<std::size_t>(payloadDimensions[1]) *
                             static_cast<std::size_t>(k));
                if (sourceIndex >= payloadCount) {
                    input->Fail("Tecplot ordered cell payload layout is invalid");
                }
                values->SetValue(targetIndex++,
                                 source->GetValue(static_cast<vtkIdType>(sourceIndex)));
            }
        }
    }
    if (targetIndex != static_cast<vtkIdType>(header.cellCount)) {
        input->Fail("Tecplot ordered cell count does not match its payload layout");
    }
    return values;
}

vtkSmartPointer<vtkDataSet> CreateGrid(
    const TecplotZoneHeader& header,
    const CoordinateVariables& coordinates,
    const std::vector<std::string>& variableNames,
    const std::vector<vtkSmartPointer<vtkDoubleArray>>& values,
    BinaryInput* input) {
    if (coordinates.x < 0 || coordinates.y < 0) {
        input->Fail("Tecplot binary input requires point variables named X and Y "
                    "(or CoordinateX and CoordinateY); found " +
                    VariableNameList(variableNames));
    }
    RequirePointCoordinate(header, coordinates.x, "X", input);
    RequirePointCoordinate(header, coordinates.y, "Y", input);
    RequirePointCoordinate(header, coordinates.z, "Z", input);

    const auto coordinateArray = [&](int index,
                                     const char* axis) -> vtkDoubleArray* {
        if (index < 0) {
            return nullptr;
        }
        const std::size_t variableIndex = static_cast<std::size_t>(index);
        if (variableIndex >= values.size() || values[variableIndex] == nullptr ||
            values[variableIndex]->GetNumberOfValues() !=
                static_cast<vtkIdType>(header.pointCount)) {
            input->Fail(std::string("Tecplot ") + axis +
                        " coordinate variable is unavailable or incomplete");
        }
        return values[variableIndex];
    };

    vtkDoubleArray* x = coordinateArray(coordinates.x, "X");
    vtkDoubleArray* y = coordinateArray(coordinates.y, "Y");
    vtkDoubleArray* z = coordinateArray(coordinates.z, "Z");

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(static_cast<vtkIdType>(header.pointCount));
    for (std::size_t index = 0; index < header.pointCount; ++index) {
        const vtkIdType pointId = static_cast<vtkIdType>(index);
        points->SetPoint(pointId,
                         x->GetValue(pointId),
                         y->GetValue(pointId),
                         z == nullptr ? 0.0 : z->GetValue(pointId));
    }

    vtkSmartPointer<vtkDataSet> grid;
    if (header.zoneType == kTecplotZoneOrdered) {
        auto structured = vtkSmartPointer<vtkStructuredGrid>::New();
        structured->SetDimensions(header.orderedDimensions.data());
        structured->SetPoints(points);
        grid = structured;
    } else {
        auto unstructured = vtkSmartPointer<vtkUnstructuredGrid>::New();
        unstructured->SetPoints(points);
        grid = unstructured;
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (static_cast<int>(index) == coordinates.x ||
            static_cast<int>(index) == coordinates.y ||
            static_cast<int>(index) == coordinates.z ||
            values[index] == nullptr) {
            continue;
        }
        if (header.variableLocations[index] == 0) {
            grid->GetPointData()->AddArray(values[index]);
        } else {
            grid->GetCellData()->AddArray(values[index]);
        }
    }
    return grid;
}

void SkipFaceNeighborData(BinaryInput* input,
                          const TecplotZoneHeader& header) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot face-neighbor reader error.");
    }
    if (header.rawFaceNeighbors) {
        const int facesPerCell = FacesPerCell(header.zoneType);
        if (facesPerCell <= 0 ||
            header.cellCount > std::numeric_limits<std::size_t>::max() /
                                   static_cast<std::size_t>(facesPerCell)) {
            input->Fail("Tecplot raw face-neighbor array size is invalid");
        }
        input->SkipInt32Values(
            header.cellCount * static_cast<std::size_t>(facesPerCell),
            "Tecplot raw face-neighbor array");
    }
    for (int record = 0;
         record < header.miscellaneousFaceNeighborCount;
         ++record) {
        const int currentCell =
            input->ReadInt32("Tecplot miscellaneous face-neighbor cell");
        const int currentFace =
            input->ReadInt32("Tecplot miscellaneous face-neighbor face");
        (void)currentCell;
        (void)currentFace;
        switch (header.miscellaneousFaceNeighborMode) {
        case 0:
            (void)input->ReadInt32(
                "Tecplot local one-to-one face-neighbor cell");
            break;
        case 1: {
            (void)input->ReadInt32(
                "Tecplot local one-to-many face-neighbor obscuration");
            const int count = input->ReadInt32(
                "Tecplot local one-to-many face-neighbor count");
            if (count < 0) {
                input->Fail(
                    "Tecplot local one-to-many face-neighbor count is invalid");
            }
            input->SkipInt32Values(
                static_cast<std::size_t>(count),
                "Tecplot local one-to-many face-neighbor cells");
            break;
        }
        case 2:
            (void)input->ReadInt32(
                "Tecplot global one-to-one face-neighbor zone");
            (void)input->ReadInt32(
                "Tecplot global one-to-one face-neighbor cell");
            break;
        case 3: {
            (void)input->ReadInt32(
                "Tecplot global one-to-many face-neighbor obscuration");
            const int count = input->ReadInt32(
                "Tecplot global one-to-many face-neighbor count");
            if (count < 0 ||
                static_cast<std::size_t>(count) >
                    std::numeric_limits<std::size_t>::max() / 2u) {
                input->Fail(
                    "Tecplot global one-to-many face-neighbor count is invalid");
            }
            input->SkipInt32Values(
                static_cast<std::size_t>(count) * 2u,
                "Tecplot global one-to-many face-neighbor pairs");
            break;
        }
        default:
            input->Fail("Tecplot miscellaneous face-neighbor mode is invalid");
        }
    }
}

struct FaceOwner {
    std::size_t faceIndex = 0;
    bool reverse = false;
};

struct PolyBoundaryMap {
    std::vector<std::size_t> offsets;
    std::vector<std::int32_t> elements;
    std::vector<std::int32_t> zones;
};

std::vector<std::size_t> ReadFaceNodeCounts(
    BinaryInput* input,
    const TecplotZoneHeader& header) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot face-map reader error.");
    }
    std::vector<std::size_t> counts(header.faceCount, 2u);
    if (header.zoneType == kTecplotZonePolygon) {
        return counts;
    }

    if (header.headerVersion == 112) {
        if (header.faceCount == std::numeric_limits<std::size_t>::max()) {
            input->Fail("Tecplot face count is too large for its offset array");
        }
        std::vector<std::size_t> offsets(header.faceCount + 1u, 0u);
        for (std::size_t index = 0; index < offsets.size(); ++index) {
            const std::int32_t rawOffset =
                input->ReadInt32("Tecplot face-node offset");
            if (rawOffset < 0) {
                input->Fail("Tecplot face-node offset is invalid");
            }
            offsets[index] = static_cast<std::size_t>(rawOffset);
            if ((index == 0 && offsets[index] != 0u) ||
                (index > 0 && offsets[index] < offsets[index - 1u]) ||
                offsets[index] > header.faceNodeCount) {
                input->Fail("Tecplot face-node offsets are invalid");
            }
        }
        if (offsets.back() != header.faceNodeCount) {
            input->Fail("Tecplot face-node offsets do not match the header");
        }
        for (std::size_t index = 0; index < counts.size(); ++index) {
            counts[index] = offsets[index + 1u] - offsets[index];
            if (counts[index] < 3u) {
                input->Fail("Tecplot FEPOLYHEDRON face has fewer than three nodes");
            }
        }
        return counts;
    }

    std::size_t total = 0;
    for (std::size_t index = 0; index < counts.size(); ++index) {
        const std::int32_t rawCount =
            input->ReadInt32("Tecplot face-node count");
        if (rawCount < 3) {
            input->Fail("Tecplot FEPOLYHEDRON face has fewer than three nodes");
        }
        counts[index] = static_cast<std::size_t>(rawCount);
        if (counts[index] > header.faceNodeCount - total) {
            input->Fail("Tecplot face-node counts do not match the header");
        }
        total += counts[index];
    }
    if (total != header.faceNodeCount) {
        input->Fail("Tecplot face-node counts do not match the header");
    }
    return counts;
}

std::vector<std::vector<vtkIdType>> ReadFaceNodes(
    BinaryInput* input,
    const TecplotZoneHeader& header,
    const std::vector<std::size_t>& faceNodeCounts) {
    if (input == nullptr || faceNodeCounts.size() != header.faceCount) {
        throw std::runtime_error("Internal Tecplot face-node reader error.");
    }
    std::vector<std::vector<vtkIdType>> faces;
    faces.reserve(header.faceCount);
    for (const std::size_t nodeCount : faceNodeCounts) {
        std::vector<vtkIdType> nodes;
        nodes.reserve(nodeCount);
        for (std::size_t node = 0; node < nodeCount; ++node) {
            const std::int32_t rawNode =
                input->ReadInt32("Tecplot face-map node");
            if (rawNode < 0 ||
                static_cast<std::uint64_t>(rawNode) >=
                    static_cast<std::uint64_t>(header.pointCount)) {
                input->Fail("Tecplot face-map uses an invalid zero-based node index");
            }
            nodes.push_back(static_cast<vtkIdType>(rawNode));
        }
        faces.push_back(std::move(nodes));
    }
    return faces;
}

PolyBoundaryMap ReadPolyBoundaryMap(BinaryInput* input,
                                    const TecplotZoneHeader& header) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot boundary-map reader error.");
    }
    PolyBoundaryMap boundary;
    if (header.boundaryFaceCount == 0u) {
        if (header.boundaryConnectionCount != 0u) {
            input->Fail("Tecplot boundary connections require boundary faces");
        }
        return boundary;
    }
    if (header.boundaryFaceCount == std::numeric_limits<std::size_t>::max()) {
        input->Fail("Tecplot boundary-face count is too large");
    }
    boundary.offsets.resize(header.boundaryFaceCount + 1u);
    for (std::size_t index = 0; index < boundary.offsets.size(); ++index) {
        const std::int32_t rawOffset =
            input->ReadInt32("Tecplot boundary-face connection offset");
        if (rawOffset < 0) {
            input->Fail("Tecplot boundary-face connection offset is invalid");
        }
        boundary.offsets[index] = static_cast<std::size_t>(rawOffset);
        if ((index == 0 && boundary.offsets[index] != 0u) ||
            (index == 1 && boundary.offsets[index] != 0u) ||
            (index > 0 && boundary.offsets[index] < boundary.offsets[index - 1u]) ||
            boundary.offsets[index] > header.boundaryConnectionCount) {
            input->Fail("Tecplot boundary-face connection offsets are invalid");
        }
    }
    if (boundary.offsets.back() != header.boundaryConnectionCount) {
        input->Fail("Tecplot boundary-face offsets do not match the header");
    }

    boundary.elements.resize(header.boundaryConnectionCount);
    boundary.zones.resize(header.boundaryConnectionCount);
    for (std::size_t index = 0; index < boundary.elements.size(); ++index) {
        boundary.elements[index] =
            input->ReadInt32("Tecplot boundary-face connection element");
    }
    for (std::size_t index = 0; index < boundary.zones.size(); ++index) {
        boundary.zones[index] =
            input->ReadInt32("Tecplot boundary-face connection zone");
        if (boundary.zones[index] < -1) {
            input->Fail("Tecplot boundary-face connection zone is invalid");
        }
    }
    return boundary;
}

void AppendFaceOwner(std::vector<std::vector<FaceOwner>>* cellFaces,
                     const TecplotZoneHeader& header,
                     std::int32_t rawCell,
                     std::size_t faceIndex,
                     bool reverse,
                     BinaryInput* input) {
    if (cellFaces == nullptr || input == nullptr) {
        throw std::runtime_error("Internal Tecplot face-owner reader error.");
    }
    if (rawCell < 0 || static_cast<std::uint64_t>(rawCell) >=
                           static_cast<std::uint64_t>(header.cellCount)) {
        input->Fail("Tecplot face-map uses an invalid zero-based cell index");
    }
    (*cellFaces)[static_cast<std::size_t>(rawCell)].push_back(
        FaceOwner{faceIndex, reverse});
}

void AppendFaceSideOwners(std::vector<std::vector<FaceOwner>>* cellFaces,
                          const TecplotZoneHeader& header,
                          const PolyBoundaryMap& boundary,
                          std::int32_t rawCell,
                          std::size_t faceIndex,
                          bool reverse,
                          BinaryInput* input) {
    if (rawCell >= 0) {
        AppendFaceOwner(cellFaces, header, rawCell, faceIndex, reverse, input);
        return;
    }
    if (rawCell == -1) {
        return;
    }
    const std::int64_t rawOffset = -static_cast<std::int64_t>(rawCell);
    if (rawOffset < 2 ||
        static_cast<std::uint64_t>(rawOffset) >=
            static_cast<std::uint64_t>(boundary.offsets.size())) {
        input->Fail("Tecplot face-map boundary reference is invalid");
    }
    const std::size_t offset = static_cast<std::size_t>(rawOffset);
    const std::size_t begin = boundary.offsets[offset - 1u];
    const std::size_t end = boundary.offsets[offset];
    for (std::size_t index = begin; index < end; ++index) {
        const std::int32_t boundaryCell = boundary.elements[index];
        if (boundaryCell < -1) {
            input->Fail("Tecplot boundary-face connection element is invalid");
        }
        if (boundary.zones[index] == -1 && boundaryCell >= 0) {
            AppendFaceOwner(cellFaces, header, boundaryCell, faceIndex,
                            reverse, input);
        }
    }
}

std::vector<vtkIdType> MakePolygonLoop(
    const std::vector<FaceOwner>& owners,
    const std::vector<std::vector<vtkIdType>>& faces,
    BinaryInput* input) {
    if (input == nullptr || owners.size() < 3u) {
        if (input != nullptr) {
            input->Fail("Tecplot FEPOLYGON cell does not have a closed boundary");
        }
        throw std::runtime_error("Internal Tecplot polygon reader error.");
    }
    std::vector<bool> used(owners.size(), false);
    const FaceOwner& firstOwner = owners.front();
    if (firstOwner.faceIndex >= faces.size() ||
        faces[firstOwner.faceIndex].size() != 2u) {
        input->Fail("Tecplot FEPOLYGON face is malformed");
    }
    const std::vector<vtkIdType>& firstFace = faces[firstOwner.faceIndex];
    vtkIdType first = firstFace[firstOwner.reverse ? 1u : 0u];
    vtkIdType current = firstFace[firstOwner.reverse ? 0u : 1u];
    std::vector<vtkIdType> loop = {first, current};
    used[0] = true;

    for (std::size_t edgeCount = 1; edgeCount < owners.size(); ++edgeCount) {
        bool found = false;
        for (std::size_t candidate = 1; candidate < owners.size(); ++candidate) {
            if (used[candidate]) {
                continue;
            }
            const FaceOwner& owner = owners[candidate];
            if (owner.faceIndex >= faces.size() || faces[owner.faceIndex].size() != 2u) {
                input->Fail("Tecplot FEPOLYGON face is malformed");
            }
            const std::vector<vtkIdType>& face = faces[owner.faceIndex];
            vtkIdType start = face[owner.reverse ? 1u : 0u];
            vtkIdType end = face[owner.reverse ? 0u : 1u];
            if (start != current && end != current) {
                continue;
            }
            const vtkIdType next = start == current ? end : start;
            if ((edgeCount + 1u != owners.size() && next == first) ||
                next == current) {
                continue;
            }
            used[candidate] = true;
            loop.push_back(next);
            current = next;
            found = true;
            break;
        }
        if (!found) {
            input->Fail("Tecplot FEPOLYGON faces do not form a closed cell");
        }
    }
    if (current != first || loop.size() != owners.size() + 1u) {
        input->Fail("Tecplot FEPOLYGON faces do not form a closed cell");
    }
    loop.pop_back();
    return loop;
}

void ReadFaceBasedConnectivity(BinaryInput* input,
                               vtkUnstructuredGrid* grid,
                               const TecplotZoneHeader& header) {
    if (input == nullptr || grid == nullptr) {
        throw std::runtime_error("Internal Tecplot face-map reader error.");
    }
    const std::vector<std::size_t> faceNodeCounts =
        ReadFaceNodeCounts(input, header);
    const std::vector<std::vector<vtkIdType>> faces =
        ReadFaceNodes(input, header, faceNodeCounts);
    std::vector<std::int32_t> left(header.faceCount);
    std::vector<std::int32_t> right(header.faceCount);
    for (std::size_t index = 0; index < header.faceCount; ++index) {
        left[index] = input->ReadInt32("Tecplot face-map left element");
    }
    for (std::size_t index = 0; index < header.faceCount; ++index) {
        right[index] = input->ReadInt32("Tecplot face-map right element");
    }
    const PolyBoundaryMap boundary = ReadPolyBoundaryMap(input, header);

    std::vector<std::vector<FaceOwner>> cellFaces(header.cellCount);
    for (std::size_t face = 0; face < header.faceCount; ++face) {
        // Tecplot orders the face nodes so that the face normal points to the
        // right element. VTK polyhedron faces need to point out of each cell.
        AppendFaceSideOwners(&cellFaces, header, boundary, left[face], face,
                             false, input);
        AppendFaceSideOwners(&cellFaces, header, boundary, right[face], face,
                             true, input);
    }

    for (std::size_t cell = 0; cell < cellFaces.size(); ++cell) {
        const std::vector<FaceOwner>& owners = cellFaces[cell];
        if (header.zoneType == kTecplotZonePolygon) {
            const std::vector<vtkIdType> points =
                MakePolygonLoop(owners, faces, input);
            grid->InsertNextCell(VTK_POLYGON,
                                 static_cast<vtkIdType>(points.size()),
                                 points.data());
            continue;
        }

        if (owners.size() < 4u) {
            input->Fail("Tecplot FEPOLYHEDRON cell has fewer than four faces");
        }
        std::vector<vtkIdType> pointIds;
        std::vector<vtkIdType> faceStream;
        for (const FaceOwner& owner : owners) {
            if (owner.faceIndex >= faces.size() ||
                faces[owner.faceIndex].size() < 3u) {
                input->Fail("Tecplot FEPOLYHEDRON face is malformed");
            }
            const std::vector<vtkIdType>& face = faces[owner.faceIndex];
            faceStream.push_back(static_cast<vtkIdType>(face.size()));
            if (owner.reverse) {
                for (auto node = face.rbegin(); node != face.rend(); ++node) {
                    faceStream.push_back(*node);
                    pointIds.push_back(*node);
                }
            } else {
                for (const vtkIdType node : face) {
                    faceStream.push_back(node);
                    pointIds.push_back(node);
                }
            }
        }
        std::sort(pointIds.begin(), pointIds.end());
        pointIds.erase(std::unique(pointIds.begin(), pointIds.end()),
                       pointIds.end());
        if (pointIds.size() < 4u) {
            input->Fail("Tecplot FEPOLYHEDRON cell has too few unique nodes");
        }
        grid->InsertNextCell(VTK_POLYHEDRON,
                             static_cast<vtkIdType>(pointIds.size()),
                             pointIds.data(),
                             static_cast<vtkIdType>(owners.size()),
                             faceStream.data());
    }
}

void ReadConnectivity(BinaryInput* input,
                      DecodedZone* zone,
                      const std::vector<DecodedZone>& priorZones,
                      int sharedConnectivityZone) {
    if (input == nullptr || zone == nullptr || zone->grid == nullptr) {
        throw std::runtime_error("Internal Tecplot connectivity reader error.");
    }
    const TecplotZoneHeader& header = zone->header;
    if (header.zoneType == kTecplotZoneOrdered) {
        if (sharedConnectivityZone >= 0) {
            const std::size_t sourceIndex =
                static_cast<std::size_t>(sharedConnectivityZone);
            if (sourceIndex >= priorZones.size()) {
                input->Fail("Tecplot ordered face neighbors reference a zone that has not been read");
            }
            const DecodedZone& source = priorZones[sourceIndex];
            if (source.header.zoneType != kTecplotZoneOrdered ||
                source.header.orderedDimensions != header.orderedDimensions) {
                input->Fail("Tecplot ordered face-neighbor sharing is incompatible with its zone");
            }
            return;
        }
        SkipFaceNeighborData(input, header);
        return;
    }

    vtkUnstructuredGrid* grid =
        vtkUnstructuredGrid::SafeDownCast(zone->grid);
    if (grid == nullptr) {
        input->Fail("Tecplot finite-element zone did not create an unstructured grid");
    }
    if (sharedConnectivityZone >= 0) {
        const std::size_t sourceIndex =
            static_cast<std::size_t>(sharedConnectivityZone);
        if (sourceIndex >= priorZones.size() ||
            priorZones[sourceIndex].grid == nullptr) {
            input->Fail("Tecplot connectivity references a zone that has not been read");
        }
        const DecodedZone& source = priorZones[sourceIndex];
        vtkUnstructuredGrid* sourceGrid =
            vtkUnstructuredGrid::SafeDownCast(source.grid);
        if (source.header.zoneType != header.zoneType ||
            source.header.pointCount != header.pointCount ||
            source.header.cellCount != header.cellCount ||
            source.header.faceCount != header.faceCount ||
            source.header.faceNodeCount != header.faceNodeCount ||
            sourceGrid == nullptr) {
            input->Fail("Tecplot shared connectivity is incompatible with its zone");
        }
        if (
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 5, 0)
            sourceGrid->GetPolyhedronFaceLocations() != nullptr &&
            sourceGrid->GetPolyhedronFaces() != nullptr
#else
            sourceGrid->GetFaceLocations() != nullptr &&
            sourceGrid->GetFaces() != nullptr
#endif
        ) {
#if VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 5, 0)
            grid->SetPolyhedralCells(sourceGrid->GetCellTypesArray(),
                                     sourceGrid->GetCells(),
                                     sourceGrid->GetPolyhedronFaceLocations(),
                                     sourceGrid->GetPolyhedronFaces());
#else
            grid->SetCells(sourceGrid->GetCellTypesArray(), sourceGrid->GetCells(),
                           sourceGrid->GetFaceLocations(), sourceGrid->GetFaces());
#endif
        } else {
            grid->SetCells(sourceGrid->GetCellTypesArray(), sourceGrid->GetCells());
        }
        return;
    }

    if (IsFaceBasedZone(header.zoneType)) {
        ReadFaceBasedConnectivity(input, grid, header);
        return;
    }

    const int nodeCount = NodesPerCell(header.zoneType);
    const int cellType = VtkCellType(header.zoneType);
    if (nodeCount <= 0 || cellType == VTK_EMPTY_CELL) {
        input->Fail("Tecplot finite-element zone type is unsupported");
    }
    std::vector<vtkIdType> nodeIds(static_cast<std::size_t>(nodeCount));
    for (std::size_t cellIndex = 0; cellIndex < header.cellCount; ++cellIndex) {
        for (int localIndex = 0; localIndex < nodeCount; ++localIndex) {
            const std::int32_t rawNodeId =
                input->ReadInt32("Tecplot finite-element connectivity");
            if (rawNodeId < 0 ||
                static_cast<std::uint64_t>(rawNodeId) >=
                    static_cast<std::uint64_t>(header.pointCount)) {
                input->Fail("Tecplot finite-element connectivity uses an invalid "
                            "zero-based node index");
            }
            nodeIds[static_cast<std::size_t>(localIndex)] =
                static_cast<vtkIdType>(rawNodeId);
        }
        grid->InsertNextCell(
            cellType, static_cast<vtkIdType>(nodeIds.size()), nodeIds.data());
    }
    SkipFaceNeighborData(input, header);
}

DecodedZone ReadZoneData(
    BinaryInput* input,
    const TecplotZoneHeader& header,
    const std::vector<std::string>& variableNames,
    const CoordinateVariables& coordinates,
    const std::vector<DecodedZone>& priorZones) {
    if (input == nullptr) {
        throw std::runtime_error("Internal Tecplot zone-data reader error.");
    }
    const float marker = input->ReadFloat("Tecplot data-zone marker");
    const float expectedMarker =
        header.headerVersion == 191 ? kTecplotZoneHeaderV191Marker
                                    : kTecplotZoneHeaderV112Marker;
    if (!MarkerMatches(marker, expectedMarker)) {
        input->Fail("Tecplot data-zone marker does not match the zone header");
    }

    const std::size_t variableCount = variableNames.size();
    std::vector<int> dataFormats(variableCount, 0);
    for (std::size_t index = 0; index < variableCount; ++index) {
        dataFormats[index] =
            input->ReadInt32("Tecplot variable data format");
        if (dataFormats[index] < 1 || dataFormats[index] > 6) {
            input->Fail("Tecplot variable data format is invalid");
        }
    }

    const int hasPassiveVariables =
        input->ReadInt32("Tecplot passive-variable flag");
    if (hasPassiveVariables != 0 && hasPassiveVariables != 1) {
        input->Fail("Tecplot passive-variable flag is invalid");
    }
    std::vector<bool> passive(variableCount, false);
    if (hasPassiveVariables != 0) {
        for (std::size_t index = 0; index < variableCount; ++index) {
            const int value =
                input->ReadInt32("Tecplot passive-variable value");
            if (value != 0 && value != 1) {
                input->Fail("Tecplot passive-variable value is invalid");
            }
            passive[index] = value != 0;
        }
    }

    const int hasVariableSharing =
        input->ReadInt32("Tecplot variable-sharing flag");
    if (hasVariableSharing != 0 && hasVariableSharing != 1) {
        input->Fail("Tecplot variable-sharing flag is invalid");
    }
    std::vector<int> sharedVariableZones(variableCount, -1);
    if (hasVariableSharing != 0) {
        for (std::size_t index = 0; index < variableCount; ++index) {
            sharedVariableZones[index] =
                input->ReadInt32("Tecplot variable-sharing zone");
            if (sharedVariableZones[index] < -1) {
                input->Fail("Tecplot variable-sharing zone is invalid");
            }
        }
    }
    const int sharedConnectivityZone =
        input->ReadInt32("Tecplot connectivity-sharing zone");
    if (sharedConnectivityZone < -1) {
        input->Fail("Tecplot connectivity-sharing zone is invalid");
    }

    for (std::size_t index = 0; index < variableCount; ++index) {
        if (!passive[index] && sharedVariableZones[index] < 0) {
            (void)input->ReadDouble("Tecplot variable minimum");
            (void)input->ReadDouble("Tecplot variable maximum");
        }
    }

    DecodedZone zone;
    zone.header = header;
    zone.variables.resize(variableCount);
    for (std::size_t index = 0; index < variableCount; ++index) {
        if (passive[index]) {
            continue;
        }
        const int sharedZone = sharedVariableZones[index];
        const std::size_t expectedValueCount =
            header.variableLocations[index] == 0 ? header.pointCount
                                                 : header.cellCount;
        if (sharedZone >= 0) {
            const std::size_t sourceIndex =
                static_cast<std::size_t>(sharedZone);
            if (sourceIndex >= priorZones.size() ||
                index >= priorZones[sourceIndex].variables.size() ||
                priorZones[sourceIndex].variables[index] == nullptr) {
                input->Fail("Tecplot variable references a zone that has not been read");
            }
            if (priorZones[sourceIndex].header.variableLocations[index] !=
                header.variableLocations[index]) {
                input->Fail("Tecplot shared variable changes point/cell association");
            }
            if (priorZones[sourceIndex].variables[index]->GetNumberOfValues() !=
                static_cast<vtkIdType>(expectedValueCount)) {
                input->Fail("Tecplot shared variable has an incompatible tuple count");
            }
            zone.variables[index] = priorZones[sourceIndex].variables[index];
            continue;
        }
        const bool orderedCellValues =
            header.zoneType == kTecplotZoneOrdered &&
            header.variableLocations[index] != 0;
        const std::size_t storedValueCount = orderedCellValues
            ? OrderedCellPayloadCount(header, input)
            : expectedValueCount;
        vtkSmartPointer<vtkDoubleArray> values = input->ReadNumericArray(
            dataFormats[index], storedValueCount, variableNames[index],
            "Tecplot variable \"" + variableNames[index] + "\"");
        if (orderedCellValues) {
            values = CompactOrderedCellValues(header, values, input);
        }
        zone.variables[index] = values;
    }

    zone.grid = CreateGrid(header, coordinates, variableNames,
                           zone.variables, input);
    ReadConnectivity(input, &zone, priorZones, sharedConnectivityZone);
    return zone;
}

}  // namespace

bool IsClassicTecplotBinaryFile(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        return false;
    }
    std::array<char, 8> magic{};
    input.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    if (input.gcount() != static_cast<std::streamsize>(magic.size())) {
        return false;
    }
    const std::string identifier(magic.data(), magic.size());
    return identifier == "#!TDV112" || identifier == "#!TDV191";
}

vtkSmartPointer<vtkDataObject> ReadClassicTecplotBinary(const fs::path& path) {
    BinaryInput input(path);
    const int fileType = input.ReadInt32("Tecplot file type");
    if (fileType != kTecplotFileTypeFull && fileType != kTecplotFileTypeGrid &&
        fileType != kTecplotFileTypeSolution) {
        input.Fail("Tecplot file type is invalid");
    }
    if (fileType == kTecplotFileTypeSolution) {
        input.Fail("Tecplot solution-only PLT input is not self-contained; "
                   "load or export a FULL PLT file with its grid");
    }
    (void)input.ReadString("Tecplot dataset title");
    const std::int32_t rawVariableCount =
        input.ReadInt32("Tecplot variable count");
    const std::size_t variableCount = PositiveCount(
        rawVariableCount, "Tecplot variable count", &input);
    std::vector<std::string> variableNames;
    variableNames.reserve(variableCount);
    for (std::size_t index = 0; index < variableCount; ++index) {
        variableNames.push_back(
            input.ReadString("Tecplot variable name"));
    }
    const CoordinateVariables coordinates =
        FindCoordinateVariables(variableNames);

    std::vector<TecplotZoneHeader> headers;
    while (true) {
        const float marker = input.ReadFloat("Tecplot header marker");
        if (MarkerMatches(marker, kTecplotEndOfHeaderMarker)) {
            break;
        }
        if (MarkerMatches(marker, kTecplotZoneHeaderV112Marker) ||
            MarkerMatches(marker, kTecplotZoneHeaderV191Marker)) {
            const int markerVersion =
                MarkerMatches(marker, kTecplotZoneHeaderV191Marker) ? 191 : 112;
            if (markerVersion == 191 && input.version() != 191) {
                input.Fail("Tecplot v112 files cannot contain v191 zone headers");
            }
            headers.push_back(
                ReadZoneHeader(&input, markerVersion, variableCount));
            continue;
        }
        SkipHeaderRecord(&input, marker);
    }
    if (headers.empty()) {
        input.Fail("Tecplot binary file contains no zones");
    }

    std::vector<DecodedZone> zones;
    zones.reserve(headers.size());
    for (const TecplotZoneHeader& header : headers) {
        zones.push_back(ReadZoneData(&input, header, variableNames,
                                     coordinates, zones));
    }

    auto output = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    output->SetNumberOfBlocks(static_cast<unsigned int>(zones.size()));
    for (std::size_t index = 0; index < zones.size(); ++index) {
        output->SetBlock(static_cast<unsigned int>(index), zones[index].grid);
    }
    return output;
}

}  // namespace streamcenterplus::io::detail
