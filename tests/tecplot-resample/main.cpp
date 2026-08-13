#include "FlowDataInspector.h"
#include <streamcenterplus/io/DataSetReader.h>
#include <streamcenterplus_resample/Core.h>

#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>

#include <QDir>
#include <QString>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace io = streamcenterplus::io;
namespace resample = streamcenterplus::resample;

namespace {

constexpr std::array<std::array<double, 3>, 8> kPoints = {{
    {{0.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}},
    {{1.0, 1.0, 0.0}}, {{0.0, 1.0, 0.0}},
    {{0.0, 0.0, 1.0}}, {{1.0, 0.0, 1.0}},
    {{1.0, 1.0, 1.0}}, {{0.0, 1.0, 1.0}},
}};

constexpr std::array<std::array<double, 3>, 4> kPolygonPoints = {{
    {{0.0, 0.0, 0.0}}, {{1.0, 0.0, 0.0}},
    {{1.0, 1.0, 0.0}}, {{0.0, 1.0, 0.0}},
}};

constexpr std::array<const char*, 7> kVariableNames = {
    "X", "Y", "Z", "U", "V", "W", "Pressure",
};

constexpr std::array<const char*, 8> kV191VariableNames = {
    "X", "Y", "Z", "U", "V", "W", "Pressure", "CellMarker",
};

constexpr std::array<std::int32_t, 8> kV191VariableLocations = {
    0, 0, 0, 0, 0, 0, 0, 1,
};

constexpr std::array<const char*, 8> kOrderedVariableNames = {
    "X", "Y", "Z", "U", "V", "W", "Pressure", "CellMarker",
};

constexpr std::array<std::int32_t, 8> kOrderedVariableLocations = {
    0, 0, 0, 0, 0, 0, 0, 1,
};

constexpr std::array<std::int32_t, 3> kOrderedDimensions = {{3, 2, 2}};

constexpr std::array<float, 6> kOrderedCellPayload = {{
    17.0F, 23.0F, 0.0F, 0.0F, 0.0F, 0.0F,
}};

constexpr double kCellMarker = 17.0;

constexpr int kTecplotZonePolygon = 6;
constexpr int kTecplotZonePolyhedron = 7;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool NearlyEqual(double lhs, double rhs, double tolerance = 1.0e-8) {
    return std::abs(lhs - rhs) <= tolerance *
           std::max({1.0, std::abs(lhs), std::abs(rhs)});
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

template <typename T>
void WriteEndian(std::ofstream& output, T value, bool littleEndian) {
    if (HostIsLittleEndian() != littleEndian) {
        value = ByteSwap(value);
    }
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template <typename T>
void WriteLittleEndian(std::ofstream& output, T value) {
    WriteEndian(output, value, true);
}

void WriteTecplotString(std::ofstream& output,
                        const std::string& value,
                        bool littleEndian = true) {
    for (const unsigned char ch : value) {
        WriteEndian<std::int32_t>(output, static_cast<std::int32_t>(ch),
                                  littleEndian);
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);
}

std::array<double, 7> ValuesForPoint(const std::array<double, 3>& point) {
    const double x = point[0];
    const double y = point[1];
    const double z = point[2];
    return {x, y, z, x + 1.0, y + 2.0, z + 3.0, x + y + z};
}

void WriteAsciiTecplot(const fs::path& path) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot create ASCII Tecplot fixture.");
    output << "TITLE = \"SCP Tecplot ASCII fixture\"\n"
           << "VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"Pressure\"\n"
           << "ZONE T=\"volume\", N=8, E=1, DATAPACKING=POINT, ZONETYPE=FEBRICK\n";
    for (const std::array<double, 3>& point : kPoints) {
        const std::array<double, 7> values = ValuesForPoint(point);
        for (std::size_t index = 0; index < values.size(); ++index) {
            if (index != 0) {
                output << ' ';
            }
            output << values[index];
        }
        output << '\n';
    }
    output << "1 2 3 4 5 6 7 8\n";
    Require(static_cast<bool>(output), "Failed while writing ASCII Tecplot fixture.");
}

void WriteAsciiTecplotPolygon(const fs::path& path) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot create ASCII polygon Tecplot fixture.");
    output << "TITLE = \"SCP Tecplot ASCII polygon fixture\"\n"
           << "VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"Pressure\"\n"
           << "ZONE T=\"polygon\", N=4, E=1, FACES=4, DATAPACKING=BLOCK, "
              "ZONETYPE=FEPOLYGON\n";
    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        for (std::size_t point = 0; point < kPolygonPoints.size(); ++point) {
            if (point != 0) {
                output << ' ';
            }
            output << ValuesForPoint(kPolygonPoints[point])[variable];
        }
        output << '\n';
    }
    output << "1 2 2 3 3 4 4 1\n"  // face nodes are ASCII one-based
           << "0 0 0 0\n"          // no left element
           << "1 1 1 1\n";         // right element 1
    Require(static_cast<bool>(output),
            "Failed while writing ASCII polygon Tecplot fixture.");
}

void WriteAsciiTecplotOrdered(const fs::path& path) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot create ASCII ordered Tecplot fixture.");
    output << "TITLE = \"SCP Tecplot ASCII ordered fixture\"\n"
           << "VARIABLES = \"X\" \"Y\" \"Z\" \"U\" \"V\" \"W\" \"Pressure\"\n"
           << "ZONE T=\"ordered\", I=3, J=2, K=2, DATAPACKING=POINT, "
              "ZONETYPE=ORDERED\n";
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 3; ++i) {
                const std::array<double, 7> values = ValuesForPoint(
                    {{0.5 * static_cast<double>(i), static_cast<double>(j),
                      static_cast<double>(k)}});
                for (std::size_t variable = 0; variable < values.size(); ++variable) {
                    if (variable != 0) {
                        output << ' ';
                    }
                    output << values[variable];
                }
                output << '\n';
            }
        }
    }
    Require(static_cast<bool>(output),
            "Failed while writing ASCII ordered Tecplot fixture.");
}

void WriteBinaryTecplot(const fs::path& path) {
    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create binary Tecplot fixture.");

    // This is the documented classic v112 PLT layout for one FEBRICK zone.
    // The data section is variable-blocked and uses zero-based connectivity.
    constexpr char magic[] = "#!TDV112";
    output.write(magic, sizeof(magic) - 1);
    WriteLittleEndian<std::int32_t>(output, 1);       // byte-order marker
    WriteLittleEndian<std::int32_t>(output, 0);       // FULL file type
    WriteTecplotString(output, "SCP Tecplot binary fixture");
    WriteLittleEndian<std::int32_t>(output,
                                    static_cast<std::int32_t>(kVariableNames.size()));
    for (const char* variable : kVariableNames) {
        WriteTecplotString(output, variable);
    }

    WriteLittleEndian<float>(output, 299.0f);         // v112 zone header marker
    WriteTecplotString(output, "volume");
    WriteLittleEndian<std::int32_t>(output, 0);       // parent zone
    WriteLittleEndian<std::int32_t>(output, -1);      // static strand
    WriteLittleEndian<double>(output, 0.0);           // solution time
    WriteLittleEndian<std::int32_t>(output, -1);      // default zone color
    WriteLittleEndian<std::int32_t>(output, 5);       // FEBRICK
    WriteLittleEndian<std::int32_t>(output, 0);       // all variables nodal
    WriteLittleEndian<std::int32_t>(output, 0);       // no raw face neighbors
    WriteLittleEndian<std::int32_t>(output, 0);       // no misc face neighbors
    WriteLittleEndian<std::int32_t>(output, 8);       // point count
    WriteLittleEndian<std::int32_t>(output, 1);       // element count
    WriteLittleEndian<std::int32_t>(output, 0);       // ICellDim
    WriteLittleEndian<std::int32_t>(output, 0);       // JCellDim
    WriteLittleEndian<std::int32_t>(output, 0);       // KCellDim
    WriteLittleEndian<std::int32_t>(output, 0);       // no zone auxiliary data

    WriteLittleEndian<float>(output, 357.0f);         // end-of-header marker
    WriteLittleEndian<float>(output, 299.0f);         // zone data marker
    for (std::size_t index = 0; index < kVariableNames.size(); ++index) {
        WriteLittleEndian<std::int32_t>(output, 1);   // float32 variable data
    }
    WriteLittleEndian<std::int32_t>(output, 0);       // no passive variables
    WriteLittleEndian<std::int32_t>(output, 0);       // no variable sharing
    WriteLittleEndian<std::int32_t>(output, -1);      // no connectivity sharing

    std::array<std::array<double, 7>, 8> values{};
    for (std::size_t point = 0; point < kPoints.size(); ++point) {
        values[point] = ValuesForPoint(kPoints[point]);
    }
    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        for (const auto& pointValues : values) {
            minimum = std::min(minimum, pointValues[variable]);
            maximum = std::max(maximum, pointValues[variable]);
        }
        WriteLittleEndian<double>(output, minimum);
        WriteLittleEndian<double>(output, maximum);
    }
    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        for (const auto& pointValues : values) {
            WriteLittleEndian<float>(output,
                                     static_cast<float>(pointValues[variable]));
        }
    }
    for (std::int32_t point = 0; point < 8; ++point) {
        WriteLittleEndian<std::int32_t>(output, point);
    }
    Require(static_cast<bool>(output), "Failed while writing binary Tecplot fixture.");
}

void WriteV191ZoneHeader(std::ofstream& output,
                         const std::string& name,
                         double solutionTime) {
    constexpr bool kLittleEndian = false;
    WriteEndian<float>(output, 298.0f, kLittleEndian);  // v191 zone header marker
    WriteTecplotString(output, name, kLittleEndian);
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // parent zone
    WriteEndian<std::int32_t>(output, -1, kLittleEndian);   // static strand
    WriteEndian<double>(output, solutionTime, kLittleEndian);
    WriteEndian<std::int32_t>(output, -1, kLittleEndian);   // default zone color
    WriteEndian<std::int32_t>(output, 5, kLittleEndian);    // FEBRICK
    WriteEndian<std::int32_t>(output, 1, kLittleEndian);    // explicit locations
    for (const std::int32_t location : kV191VariableLocations) {
        WriteEndian<std::int32_t>(output, location, kLittleEndian);
    }
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // no raw face neighbors
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // no misc face neighbors
    WriteEndian<std::int32_t>(output, 8, kLittleEndian);    // point count
    WriteEndian<std::int32_t>(output, 1, kLittleEndian);    // element count
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // ICellDim
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // JCellDim
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // KCellDim
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);    // no zone auxiliary data
}

void WriteV191InitialZoneData(std::ofstream& output) {
    constexpr bool kLittleEndian = false;
    WriteEndian<float>(output, 298.0f, kLittleEndian);  // v191 zone data marker
    for (std::size_t index = 0; index < kV191VariableNames.size(); ++index) {
        WriteEndian<std::int32_t>(output, 1, kLittleEndian);  // float32 data
    }
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);   // no passive variables
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);   // no variable sharing
    WriteEndian<std::int32_t>(output, -1, kLittleEndian);  // no connectivity sharing

    for (std::size_t variable = 0; variable < kV191VariableNames.size(); ++variable) {
        double minimum = kCellMarker;
        double maximum = kCellMarker;
        if (variable < kVariableNames.size()) {
            minimum = std::numeric_limits<double>::infinity();
            maximum = -std::numeric_limits<double>::infinity();
            for (const std::array<double, 3>& point : kPoints) {
                const double value = ValuesForPoint(point)[variable];
                minimum = std::min(minimum, value);
                maximum = std::max(maximum, value);
            }
        }
        WriteEndian<double>(output, minimum, kLittleEndian);
        WriteEndian<double>(output, maximum, kLittleEndian);
    }
    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        for (const std::array<double, 3>& point : kPoints) {
            WriteEndian<float>(output,
                               static_cast<float>(ValuesForPoint(point)[variable]),
                               kLittleEndian);
        }
    }
    WriteEndian<float>(output, static_cast<float>(kCellMarker), kLittleEndian);
    for (std::int32_t point = 0; point < 8; ++point) {
        WriteEndian<std::int32_t>(output, point, kLittleEndian);
    }
}

void WriteV191SharedZoneData(std::ofstream& output) {
    constexpr bool kLittleEndian = false;
    WriteEndian<float>(output, 298.0f, kLittleEndian);  // v191 zone data marker
    for (std::size_t index = 0; index < kV191VariableNames.size(); ++index) {
        WriteEndian<std::int32_t>(output, 1, kLittleEndian);  // float32 data
    }
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);  // no passive variables
    WriteEndian<std::int32_t>(output, 1, kLittleEndian);  // variables shared from zone 0
    for (std::size_t index = 0; index < kV191VariableNames.size(); ++index) {
        WriteEndian<std::int32_t>(output, 0, kLittleEndian);
    }
    // The second zone has neither variable payload nor connectivity payload:
    // both are referenced from the first zone by zero-based zone index.
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);
}

void WriteV191BigEndianSharedTecplot(const fs::path& path) {
    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create v191 Tecplot fixture.");

    constexpr bool kLittleEndian = false;
    constexpr char magic[] = "#!TDV191";
    output.write(magic, sizeof(magic) - 1);
    WriteEndian<std::int32_t>(output, 1, kLittleEndian);  // byte-order marker
    WriteEndian<std::int32_t>(output, 0, kLittleEndian);  // FULL file type
    WriteTecplotString(output, "SCP Tecplot v191 big-endian fixture", kLittleEndian);
    WriteEndian<std::int32_t>(
        output, static_cast<std::int32_t>(kV191VariableNames.size()), kLittleEndian);
    for (const char* variable : kV191VariableNames) {
        WriteTecplotString(output, variable, kLittleEndian);
    }

    WriteV191ZoneHeader(output, "volume-0", 0.0);
    WriteV191ZoneHeader(output, "volume-1-shared", 0.25);
    WriteEndian<float>(output, 357.0f, kLittleEndian);  // end-of-header marker
    WriteV191InitialZoneData(output);
    WriteV191SharedZoneData(output);
    Require(static_cast<bool>(output), "Failed while writing v191 Tecplot fixture.");
}

std::array<double, 3> OrderedPoint(int i, int j, int k) {
    return {{0.5 * static_cast<double>(i), static_cast<double>(j),
             static_cast<double>(k)}};
}

void WriteOrderedZoneHeader(std::ofstream& output,
                            int version,
                            bool littleEndian) {
    WriteEndian<float>(output, version == 191 ? 298.0F : 299.0F, littleEndian);
    WriteTecplotString(output, "ordered-volume", littleEndian);
    WriteEndian<std::int32_t>(output, 0, littleEndian);    // parent zone
    WriteEndian<std::int32_t>(output, -1, littleEndian);   // static strand
    WriteEndian<double>(output, 0.0, littleEndian);        // solution time
    WriteEndian<std::int32_t>(output, -1, littleEndian);   // default zone color
    WriteEndian<std::int32_t>(output, 0, littleEndian);    // ORDERED
    WriteEndian<std::int32_t>(output, 1, littleEndian);    // explicit locations
    for (const std::int32_t location : kOrderedVariableLocations) {
        WriteEndian<std::int32_t>(output, location, littleEndian);
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);    // no raw face neighbors
    WriteEndian<std::int32_t>(output, 0, littleEndian);    // no misc face neighbors
    for (const std::int32_t dimension : kOrderedDimensions) {
        WriteEndian<std::int32_t>(output, dimension, littleEndian);
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);    // no zone auxiliary data
}

void WriteOrderedZoneData(std::ofstream& output,
                          int version,
                          bool littleEndian) {
    WriteEndian<float>(output, version == 191 ? 298.0F : 299.0F, littleEndian);
    for (std::size_t index = 0; index < kOrderedVariableNames.size(); ++index) {
        WriteEndian<std::int32_t>(output, 1, littleEndian);  // float32 variable data
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);   // no passive variables
    WriteEndian<std::int32_t>(output, 0, littleEndian);   // no variable sharing
    WriteEndian<std::int32_t>(output, -1, littleEndian);  // no connectivity sharing

    std::array<std::array<double, 7>, 12> values{};
    std::size_t pointIndex = 0;
    for (int k = 0; k < kOrderedDimensions[2]; ++k) {
        for (int j = 0; j < kOrderedDimensions[1]; ++j) {
            for (int i = 0; i < kOrderedDimensions[0]; ++i) {
                values[pointIndex++] = ValuesForPoint(OrderedPoint(i, j, k));
            }
        }
    }

    for (std::size_t variable = 0; variable < kOrderedVariableNames.size(); ++variable) {
        double minimum = std::numeric_limits<double>::infinity();
        double maximum = -std::numeric_limits<double>::infinity();
        if (variable < kVariableNames.size()) {
            for (const auto& pointValues : values) {
                minimum = std::min(minimum, pointValues[variable]);
                maximum = std::max(maximum, pointValues[variable]);
            }
        } else {
            minimum = 17.0;
            maximum = 23.0;
        }
        WriteEndian<double>(output, minimum, littleEndian);
        WriteEndian<double>(output, maximum, littleEndian);
    }

    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        for (const auto& pointValues : values) {
            WriteEndian<float>(output,
                               static_cast<float>(pointValues[variable]),
                               littleEndian);
        }
    }
    for (const float value : kOrderedCellPayload) {
        WriteEndian<float>(output, value, littleEndian);
    }
}

void WriteOrderedTecplot(const fs::path& path,
                         int version,
                         bool littleEndian) {
    Require(version == 112 || version == 191,
            "Unsupported ordered Tecplot fixture version.");
    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create ordered Tecplot fixture.");

    const char* magic = version == 191 ? "#!TDV191" : "#!TDV112";
    output.write(magic, 8);
    WriteEndian<std::int32_t>(output, 1, littleEndian);  // byte-order marker
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // FULL file type
    WriteTecplotString(output, "SCP ordered Tecplot fixture", littleEndian);
    WriteEndian<std::int32_t>(
        output, static_cast<std::int32_t>(kOrderedVariableNames.size()), littleEndian);
    for (const char* variable : kOrderedVariableNames) {
        WriteTecplotString(output, variable, littleEndian);
    }

    WriteOrderedZoneHeader(output, version, littleEndian);
    WriteEndian<float>(output, 357.0F, littleEndian);  // end-of-header marker
    WriteOrderedZoneData(output, version, littleEndian);
    Require(static_cast<bool>(output), "Failed while writing ordered Tecplot fixture.");
}

void WriteFaceBasedTecplot(const fs::path& path,
                           int version,
                           bool littleEndian,
                           int zoneType) {
    Require((version == 112 || version == 191) &&
                (zoneType == kTecplotZonePolygon ||
                 zoneType == kTecplotZonePolyhedron),
            "Unsupported face-based Tecplot fixture.");
    const bool polyhedron = zoneType == kTecplotZonePolyhedron;
    const bool hasBoundaryMap = polyhedron && version == 191;
    std::vector<std::array<double, 3>> points;
    if (polyhedron) {
        points.assign(kPoints.begin(), kPoints.end());
    } else {
        points.assign(kPolygonPoints.begin(), kPolygonPoints.end());
    }

    const std::size_t faceCount = polyhedron ? 6u : 4u;
    const std::size_t faceNodeCount = polyhedron ? 24u : 8u;
    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create face-based Tecplot fixture.");

    const char* magic = version == 191 ? "#!TDV191" : "#!TDV112";
    output.write(magic, 8);
    WriteEndian<std::int32_t>(output, 1, littleEndian);  // byte-order marker
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // FULL file type
    WriteTecplotString(output, "SCP face-based Tecplot fixture", littleEndian);
    WriteEndian<std::int32_t>(
        output, static_cast<std::int32_t>(kV191VariableNames.size()), littleEndian);
    for (const char* variable : kV191VariableNames) {
        WriteTecplotString(output, variable, littleEndian);
    }

    WriteEndian<float>(output, version == 191 ? 298.0F : 299.0F, littleEndian);
    WriteTecplotString(output, polyhedron ? "polyhedron" : "polygon", littleEndian);
    WriteEndian<std::int32_t>(output, 0, littleEndian);   // parent zone
    WriteEndian<std::int32_t>(output, -1, littleEndian);  // static strand
    WriteEndian<double>(output, 0.0, littleEndian);       // solution time
    WriteEndian<std::int32_t>(output, -1, littleEndian);  // default zone color
    WriteEndian<std::int32_t>(output, zoneType, littleEndian);
    WriteEndian<std::int32_t>(output, 1, littleEndian);  // explicit locations
    for (const std::int32_t location : kV191VariableLocations) {
        WriteEndian<std::int32_t>(output, location, littleEndian);
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // no raw face neighbors
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // no misc face neighbors
    WriteEndian<std::int32_t>(output, static_cast<std::int32_t>(points.size()),
                              littleEndian);
    if (version == 191) {
        WriteEndian<std::int64_t>(output, static_cast<std::int64_t>(faceCount),
                                  littleEndian);
        WriteEndian<std::int64_t>(output, static_cast<std::int64_t>(faceNodeCount),
                                  littleEndian);
    } else {
        WriteEndian<std::int32_t>(output, static_cast<std::int32_t>(faceCount),
                                  littleEndian);
        WriteEndian<std::int32_t>(output,
                                  static_cast<std::int32_t>(faceNodeCount),
                                  littleEndian);
    }
    WriteEndian<std::int32_t>(output, hasBoundaryMap ? 2 : 0, littleEndian);
    WriteEndian<std::int32_t>(output, hasBoundaryMap ? 1 : 0, littleEndian);
    WriteEndian<std::int32_t>(output, 1, littleEndian);  // element count
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // ICellDim
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // JCellDim
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // KCellDim
    WriteEndian<std::int32_t>(output, 0, littleEndian);  // no zone auxiliary data

    WriteEndian<float>(output, 357.0F, littleEndian);  // end-of-header marker
    WriteEndian<float>(output,  version == 191 ? 298.0F : 299.0F, littleEndian);
    for (std::size_t variable = 0; variable < kV191VariableNames.size(); ++variable) {
        WriteEndian<std::int32_t>(output, 1, littleEndian);  // float32 variable data
    }
    WriteEndian<std::int32_t>(output, 0, littleEndian);   // no passive variables
    WriteEndian<std::int32_t>(output, 0, littleEndian);   // no variable sharing
    WriteEndian<std::int32_t>(output, -1, littleEndian);  // no connectivity sharing
    for (std::size_t variable = 0; variable < kV191VariableNames.size(); ++variable) {
        double minimum = kCellMarker;
        double maximum = kCellMarker;
        if (variable < kVariableNames.size()) {
            minimum = std::numeric_limits<double>::infinity();
            maximum = -std::numeric_limits<double>::infinity();
            for (const std::array<double, 3>& point : points) {
                const double value = ValuesForPoint(point)[variable];
                minimum = std::min(minimum, value);
                maximum = std::max(maximum, value);
            }
        }
        WriteEndian<double>(output, minimum, littleEndian);
        WriteEndian<double>(output, maximum, littleEndian);
    }
    for (std::size_t variable = 0; variable < kVariableNames.size(); ++variable) {
        for (const std::array<double, 3>& point : points) {
            WriteEndian<float>(output,
                               static_cast<float>(ValuesForPoint(point)[variable]),
                               littleEndian);
        }
    }
    WriteEndian<float>(output, static_cast<float>(kCellMarker), littleEndian);

    if (polyhedron) {
        if (version == 112) {
            for (std::int32_t offset = 0; offset <= 24; offset += 4) {
                WriteEndian<std::int32_t>(output, offset, littleEndian);
            }
        } else {
            for (int face = 0; face < 6; ++face) {
                WriteEndian<std::int32_t>(output, 4, littleEndian);
            }
        }
        constexpr std::array<std::int32_t, 24> faceNodes = {{
            0, 3, 2, 1,  // z-min, outward
            4, 5, 6, 7,  // z-max, outward
            0, 1, 5, 4,  // y-min, outward
            1, 2, 6, 5,  // x-max, outward
            2, 3, 7, 6,  // y-max, outward
            3, 0, 4, 7,  // x-min, outward
        }};
        for (const std::int32_t node : faceNodes) {
            WriteEndian<std::int32_t>(output, node, littleEndian);
        }
        for (int face = 0; face < 6; ++face) {
            WriteEndian<std::int32_t>(output,
                                      hasBoundaryMap && face == 0 ? -2 : 0,
                                      littleEndian);
        }
        for (int face = 0; face < 6; ++face) {
            WriteEndian<std::int32_t>(output, -1, littleEndian);
        }
        if (hasBoundaryMap) {
            for (const std::int32_t offset : {0, 0, 1}) {
                WriteEndian<std::int32_t>(output, offset, littleEndian);
            }
            WriteEndian<std::int32_t>(output, 0, littleEndian);   // local cell 0
            WriteEndian<std::int32_t>(output, -1, littleEndian);  // current zone
        }
    } else {
        constexpr std::array<std::int32_t, 8> faceNodes = {{
            0, 1, 1, 2, 2, 3, 3, 0,
        }};
        for (const std::int32_t node : faceNodes) {
            WriteEndian<std::int32_t>(output, node, littleEndian);
        }
        for (int face = 0; face < 4; ++face) {
            WriteEndian<std::int32_t>(output, -1, littleEndian);
        }
        for (int face = 0; face < 4; ++face) {
            WriteEndian<std::int32_t>(output, 0, littleEndian);
        }
    }
    Require(static_cast<bool>(output), "Failed while writing face-based Tecplot fixture.");
}

QString ToQString(const fs::path& path) {
    return QDir::fromNativeSeparators(QString::fromStdWString(path.wstring()));
}

void ValidateDataSet(const fs::path& path) {
    Require(io::IsSupportedDataSource(path),
            "Tecplot source was not recognized by the shared reader.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(path);
    Require(frames.size() == 1 && frames.front().kind == io::DataSetFrameKind::File,
            "Tecplot source did not enumerate as a single file frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfPoints() == 8 &&
                data->GetNumberOfCells() == 1,
            "Tecplot geometry was not converted to one unstructured brick.");
    for (const char* name : {"U", "V", "W", "Pressure"}) {
        vtkDataArray* array = data->GetPointData()->GetArray(name);
        Require(array != nullptr && array->GetNumberOfComponents() == 1 &&
                    array->GetNumberOfTuples() == data->GetNumberOfPoints(),
                std::string("Tecplot nodal array is missing or malformed: ") + name);
    }
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    Require(NearlyEqual(pressure->GetComponent(6, 0), 3.0),
            "Tecplot nodal Pressure values are incorrect.");
}

void ValidateOrderedDataSet(const fs::path& path,
                            bool expectCellMarker = true) {
    Require(io::IsSupportedDataSource(path),
            "Ordered Tecplot source was not recognized by the shared reader.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(path);
    Require(frames.size() == 1 && frames.front().kind == io::DataSetFrameKind::File,
            "Ordered Tecplot source did not enumerate as a single file frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    vtkStructuredGrid* structured = vtkStructuredGrid::SafeDownCast(data);
    Require(structured != nullptr && structured->GetNumberOfPoints() == 12 &&
                structured->GetNumberOfCells() == 2,
            "Tecplot ORDERED geometry was not converted to a structured grid.");
    int dimensions[3] = {};
    structured->GetDimensions(dimensions);
    Require(dimensions[0] == 3 && dimensions[1] == 2 && dimensions[2] == 2,
            "Tecplot ORDERED dimensions are incorrect.");

    vtkDataArray* pressure = structured->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr && pressure->GetNumberOfTuples() == 12 &&
                NearlyEqual(pressure->GetComponent(11, 0), 3.0),
            "Tecplot ORDERED nodal data is incorrect.");
    if (expectCellMarker) {
        vtkDataArray* cellMarker = structured->GetCellData()->GetArray("CellMarker");
        Require(cellMarker != nullptr && cellMarker->GetNumberOfTuples() == 2 &&
                    NearlyEqual(cellMarker->GetComponent(0, 0), 17.0) &&
                    NearlyEqual(cellMarker->GetComponent(1, 0), 23.0),
                "Tecplot ORDERED cell data was not compacted from its binary payload.");
    }
}

void ValidatePolygonDataSet(const fs::path& path,
                            bool expectCellMarker = true) {
    Require(io::IsSupportedDataSource(path),
            "Polygon Tecplot source was not recognized by the shared reader.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(path);
    Require(frames.size() == 1 && frames.front().kind == io::DataSetFrameKind::File,
            "Polygon Tecplot source did not enumerate as a single file frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfPoints() == 4 &&
                data->GetNumberOfCells() == 1 &&
                data->GetCellType(0) == VTK_POLYGON,
            "Tecplot FEPOLYGON geometry was not converted to a VTK polygon.");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr && pressure->GetNumberOfTuples() == 4 &&
                NearlyEqual(pressure->GetComponent(2, 0), 2.0),
            "Tecplot FEPOLYGON nodal data is incorrect.");
    if (expectCellMarker) {
        vtkDataArray* cellMarker = data->GetCellData()->GetArray("CellMarker");
        Require(cellMarker != nullptr && cellMarker->GetNumberOfTuples() == 1 &&
                    NearlyEqual(cellMarker->GetComponent(0, 0), kCellMarker),
                "Tecplot FEPOLYGON cell data is incorrect.");
    }
}

void ValidatePolyhedronDataSet(const fs::path& path) {
    Require(io::IsSupportedDataSource(path),
            "Polyhedron Tecplot source was not recognized by the shared reader.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(path);
    Require(frames.size() == 1 && frames.front().kind == io::DataSetFrameKind::File,
            "Polyhedron Tecplot source did not enumerate as a single file frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfPoints() == 8 &&
                data->GetNumberOfCells() == 1 &&
                data->GetCellType(0) == VTK_POLYHEDRON,
            "Tecplot FEPOLYHEDRON geometry was not converted to a VTK polyhedron.");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr && pressure->GetNumberOfTuples() == 8 &&
                NearlyEqual(pressure->GetComponent(6, 0), 3.0),
            "Tecplot FEPOLYHEDRON nodal data is incorrect.");
    vtkDataArray* cellMarker = data->GetCellData()->GetArray("CellMarker");
    Require(cellMarker != nullptr && cellMarker->GetNumberOfTuples() == 1 &&
                NearlyEqual(cellMarker->GetComponent(0, 0), kCellMarker),
            "Tecplot FEPOLYHEDRON cell data is incorrect.");
}

void ValidateV191BigEndianSharedDataSet(const fs::path& path) {
    Require(io::IsSupportedDataSource(path),
            "Tecplot v191 source was not recognized by the shared reader.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(path);
    Require(frames.size() == 1 && frames.front().kind == io::DataSetFrameKind::File,
            "Tecplot v191 source did not enumerate as a single file frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfPoints() == 16 &&
                data->GetNumberOfCells() == 2,
            "Tecplot v191 shared zones were not converted to two brick blocks.");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr && pressure->GetNumberOfTuples() == 16 &&
                NearlyEqual(pressure->GetComponent(6, 0), 3.0) &&
                NearlyEqual(pressure->GetComponent(14, 0), 3.0),
            "Tecplot v191 shared nodal variables are incorrect.");
    vtkDataArray* cellMarker = data->GetCellData()->GetArray("CellMarker");
    Require(cellMarker != nullptr && cellMarker->GetNumberOfComponents() == 1 &&
                cellMarker->GetNumberOfTuples() == 2 &&
                NearlyEqual(cellMarker->GetComponent(0, 0), kCellMarker) &&
                NearlyEqual(cellMarker->GetComponent(1, 0), kCellMarker),
            "Tecplot v191 shared cell-centered variable is incorrect.");
}

void ValidateResampling(const fs::path& path, bool expectCellMarker = false) {
    resample::Parameters parameters;
    parameters.inputPath = path;
    parameters.dimensions = {3, 3, 3};
    parameters.physicalDt = 0.25;

    const resample::SequenceInfo sequence = resample::CollectInputSequence(parameters);
    Require(sequence.frames.size() == 1 && sequence.seriesBoundsPrecomputed &&
                NearlyEqual(sequence.bounds[0], 0.0) && NearlyEqual(sequence.bounds[1], 1.0) &&
                NearlyEqual(sequence.bounds[2], 0.0) && NearlyEqual(sequence.bounds[3], 1.0) &&
                NearlyEqual(sequence.bounds[4], 0.0) && NearlyEqual(sequence.bounds[5], 1.0),
            "Resampler did not discover the Tecplot file bounds.");

    vtkSmartPointer<vtkDataSet> input = resample::ReadDataSet(sequence.frames.front());
    vtkSmartPointer<vtkDataSet> prepared = resample::PreparePointData(input);
    if (expectCellMarker) {
        vtkDataArray* cellMarker = prepared->GetPointData()->GetArray("CellMarker");
        Require(cellMarker != nullptr &&
                    cellMarker->GetNumberOfTuples() == prepared->GetNumberOfPoints() &&
                    NearlyEqual(cellMarker->GetComponent(0, 0), kCellMarker),
                "Resampler did not promote the Tecplot cell-centered variable.");
    }
    vtkSmartPointer<vtkImageData> image =
        resample::ResampleFrame(prepared, parameters, sequence.bounds);
    Require(image != nullptr && image->GetNumberOfPoints() == 27,
            "Resampler did not produce the requested uniform image.");
    vtkDataArray* pressure = image->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr && NearlyEqual(pressure->GetComponent(13, 0), 1.5),
            "Resampled Tecplot pressure is incorrect at the brick centre.");
    if (expectCellMarker) {
        vtkDataArray* cellMarker = image->GetPointData()->GetArray("CellMarker");
        Require(cellMarker != nullptr &&
                    NearlyEqual(cellMarker->GetComponent(13, 0), kCellMarker),
                "Resampled Tecplot cell-centered value is incorrect at the brick centre.");
    }
}

void ValidateInspector(const fs::path& path, const QString& expectedFileType) {
    const FlowDataScanResult scan = FlowDataInspector::scanFolder(
        ToQString(path), FlowDataKind::Raw, true);
    Require(scan.passed && scan.files.size() == 1 &&
                scan.entry.dataName == QString::fromStdString(path.stem().string()) &&
                scan.entry.fileType == expectedFileType,
            "FlowDataInspector did not validate the Tecplot file input.");
}

void ValidateFolder(const fs::path& directory) {
    resample::Parameters parameters;
    parameters.inputPath = directory;
    parameters.physicalDt = 0.25;
    const resample::SequenceInfo sequence = resample::CollectInputSequence(parameters);
    Require(sequence.frames.size() == 3 &&
                sequence.frames[0].path.extension() == ".dat" &&
                sequence.frames[1].path.extension() == ".plt" &&
                sequence.frames[2].path.extension() == ".plt" &&
                NearlyEqual(sequence.frames[0].time, 0.0) &&
                NearlyEqual(sequence.frames[1].time, 0.25) &&
                NearlyEqual(sequence.frames[2].time, 0.5),
            "Tecplot folder enumeration did not preserve numeric frame order.");

    const FlowDataScanResult scan = FlowDataInspector::scanFolder(
        ToQString(directory), FlowDataKind::Raw, true);
    Require(scan.passed && scan.files.size() == 3 &&
                scan.entry.fileType == QStringLiteral("tecplot series"),
            "FlowDataInspector did not validate the mixed Tecplot folder.");
}

}  // namespace

int main() {
    try {
        const fs::path root = fs::path(SCP_RESAMPLER_IO_TEST_OUTPUT_ROOT) /
                              "tecplot-resample";
        fs::remove_all(root);
        fs::create_directories(root);
        const fs::path asciiPath = root / "frame_000.dat";
        const fs::path binaryPath = root / "frame_001.plt";
        const fs::path v191BigEndianPath = root / "frame_002_v191_big_endian.plt";
        const fs::path orderedDirectory = root / "ordered";
        fs::create_directories(orderedDirectory);
        const fs::path orderedAsciiPath = orderedDirectory / "ordered_ascii.dat";
        const fs::path orderedV112Path = orderedDirectory / "ordered_v112.plt";
        const fs::path orderedV191BigEndianPath =
            orderedDirectory / "ordered_v191_big_endian.plt";
        const fs::path polyDirectory = root / "poly";
        fs::create_directories(polyDirectory);
        const fs::path polygonAsciiPath = polyDirectory / "polygon_ascii.dat";
        const fs::path polygonV112Path = polyDirectory / "polygon_v112.plt";
        const fs::path polyhedronV112Path = polyDirectory / "polyhedron_v112.plt";
        const fs::path polyhedronV191BigEndianPath =
            polyDirectory / "polyhedron_v191_big_endian.plt";
        WriteAsciiTecplot(asciiPath);
        WriteBinaryTecplot(binaryPath);
        WriteV191BigEndianSharedTecplot(v191BigEndianPath);
        WriteAsciiTecplotOrdered(orderedAsciiPath);
        WriteOrderedTecplot(orderedV112Path, 112, true);
        WriteOrderedTecplot(orderedV191BigEndianPath, 191, false);
        WriteAsciiTecplotPolygon(polygonAsciiPath);
        WriteFaceBasedTecplot(polygonV112Path, 112, true, kTecplotZonePolygon);
        WriteFaceBasedTecplot(polyhedronV112Path, 112, true,
                              kTecplotZonePolyhedron);
        WriteFaceBasedTecplot(polyhedronV191BigEndianPath, 191, false,
                              kTecplotZonePolyhedron);

        ValidateDataSet(asciiPath);
        ValidateDataSet(binaryPath);
        ValidateV191BigEndianSharedDataSet(v191BigEndianPath);
        ValidateOrderedDataSet(orderedAsciiPath, false);
        ValidateOrderedDataSet(orderedV112Path);
        ValidateOrderedDataSet(orderedV191BigEndianPath);
        ValidatePolygonDataSet(polygonAsciiPath, false);
        ValidatePolygonDataSet(polygonV112Path);
        ValidatePolyhedronDataSet(polyhedronV112Path);
        ValidatePolyhedronDataSet(polyhedronV191BigEndianPath);
        ValidateResampling(asciiPath);
        ValidateResampling(binaryPath);
        ValidateResampling(v191BigEndianPath, true);
        ValidateResampling(orderedAsciiPath);
        ValidateResampling(orderedV112Path);
        ValidateResampling(orderedV191BigEndianPath);
        ValidateResampling(polyhedronV112Path, true);
        ValidateResampling(polyhedronV191BigEndianPath, true);
        ValidateInspector(asciiPath, QStringLiteral("tecplot ASCII series"));
        ValidateInspector(binaryPath, QStringLiteral("tecplot binary series"));
        ValidateInspector(v191BigEndianPath, QStringLiteral("tecplot binary series"));
        ValidateInspector(orderedAsciiPath, QStringLiteral("tecplot ASCII series"));
        ValidateInspector(orderedV112Path, QStringLiteral("tecplot binary series"));
        ValidateInspector(orderedV191BigEndianPath,
                          QStringLiteral("tecplot binary series"));
        ValidateInspector(polygonAsciiPath, QStringLiteral("tecplot ASCII series"));
        ValidateInspector(polygonV112Path, QStringLiteral("tecplot binary series"));
        ValidateInspector(polyhedronV112Path,
                          QStringLiteral("tecplot binary series"));
        ValidateInspector(polyhedronV191BigEndianPath,
                          QStringLiteral("tecplot binary series"));
        ValidateFolder(root);
        std::cout << "Tecplot Resample tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Tecplot Resample tests failed: " << error.what() << '\n';
        return 1;
    }
}
