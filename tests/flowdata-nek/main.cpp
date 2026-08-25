#include <FlowDataInspector.h>
#include <streamcenterplus/io/DataSetReader.h>

#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkXMLStructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace io = streamcenterplus::io;

namespace {

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

double rangeUpperBound(const QString& range) {
    const int comma = range.indexOf(QLatin1Char(','));
    bool ok = false;
    const double value = comma >= 0
        ? range.mid(comma + 1).trimmed().chopped(1).toDouble(&ok)
        : 0.0;
    Require(ok, "Cannot parse FlowData range: " + range.toStdString());
    return value;
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
void WriteValue(std::ofstream& output, T value, bool littleEndian) {
    if (littleEndian != HostIsLittleEndian()) {
        value = ByteSwap(value);
    }
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

std::size_t LocalIndex(int i, int j, int k) {
    return static_cast<std::size_t>((k * 2 + j) * 2 + i);
}

std::vector<double> Gll4() {
    const double inner = 1.0 / std::sqrt(5.0);
    return {-1.0, -inner, inner, 1.0};
}

std::size_t LocalIndex4(int i, int j, int k) {
    return static_cast<std::size_t>((k * 4 + j) * 4 + i);
}

void WriteNekField(const fs::path& path, int frameIndex, bool coordinates,
                   double physicalTime) {
    std::string header = "#std 8 2 2 2 1 1 "
        + std::to_string(physicalTime) + " " + std::to_string(frameIndex)
        + " 0 1 " + (coordinates ? "XUP" : "UP") + " 0.0 F";
    Require(header.size() <= 132, "Synthetic Nek5000 header is too long.");
    header.resize(132, ' ');

    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create synthetic Nek5000 field.");
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    WriteValue<float>(output, 6.54321f, true);
    WriteValue<std::int32_t>(output, 42, true);

    std::array<std::vector<double>, 7> arrays;
    for (std::vector<double>& array : arrays) {
        array.resize(8, 0.0);
    }
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                const std::size_t index = LocalIndex(i, j, k);
                const double x = static_cast<double>(i);
                const double y = static_cast<double>(j);
                const double z = static_cast<double>(k);
                arrays[0][index] = x;
                arrays[1][index] = y;
                arrays[2][index] = z;
                arrays[3][index] = x + physicalTime;
                arrays[4][index] = y + physicalTime;
                arrays[5][index] = z + physicalTime;
                arrays[6][index] = x + y + z + physicalTime;
            }
        }
    }
    const auto writeArray = [&](const std::vector<double>& values) {
        for (const double value : values) {
            WriteValue<double>(output, value, true);
        }
    };
    if (coordinates) {
        writeArray(arrays[0]);
        writeArray(arrays[1]);
        writeArray(arrays[2]);
    }
    writeArray(arrays[3]);
    writeArray(arrays[4]);
    writeArray(arrays[5]);
    writeArray(arrays[6]);
    Require(static_cast<bool>(output), "Failed while writing synthetic Nek5000 field.");
}

void WriteCurvedNekField(const fs::path& path) {
    std::string header = "#std 8 4 4 4 1 1 0.0 0 0 1 XUP 0.0 F";
    Require(header.size() <= 132, "Curved Nek5000 header is too long.");
    header.resize(132, ' ');

    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create curved Nek5000 field.");
    output.write(header.data(), static_cast<std::streamsize>(header.size()));
    WriteValue<float>(output, 6.54321f, true);
    WriteValue<std::int32_t>(output, 77, true);

    std::array<std::vector<double>, 7> arrays;
    for (std::vector<double>& array : arrays) {
        array.resize(64, 0.0);
    }
    const std::vector<double> nodes = Gll4();
    for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i) {
                const std::size_t index = LocalIndex4(i, j, k);
                const double r = nodes[static_cast<std::size_t>(i)];
                const double s = nodes[static_cast<std::size_t>(j)];
                const double t = nodes[static_cast<std::size_t>(k)];
                // The exact surface maximum x=2 occurs at (r,s,t)=(1,0,0),
                // which is not a GLL node for order four.  This exposes the
                // difference between a nodal AABB and spectral geometry.
                const double x = r + (1.0 - s * s) * (1.0 - t * t);
                arrays[0][index] = x;
                arrays[1][index] = s;
                arrays[2][index] = t;
                arrays[3][index] = x;
                arrays[4][index] = s;
                arrays[5][index] = t;
                arrays[6][index] = x + s + t;
            }
        }
    }
    const auto writeArray = [&](const std::vector<double>& values) {
        for (const double value : values) {
            WriteValue<double>(output, value, true);
        }
    };
    for (int component = 0; component < 3; ++component) {
        writeArray(arrays[static_cast<std::size_t>(component)]);
    }
    for (int component = 3; component < 7; ++component) {
        writeArray(arrays[static_cast<std::size_t>(component)]);
    }
    Require(static_cast<bool>(output), "Failed while writing curved Nek5000 field.");
}

fs::path CreateFixture(const fs::path& root) {
    const fs::path directory = root / "flowdata-nek5000";
    fs::create_directories(directory);
    const fs::path descriptor = directory / "flow.nek5000";
    std::ofstream descriptorOutput(descriptor);
    Require(static_cast<bool>(descriptorOutput), "Cannot create Nek5000 descriptor.");
    descriptorOutput << "filetemplate: flow%01d.f%05d\n"
                     << "firsttimestep: 0\n"
                     << "numtimesteps: 2\n";
    descriptorOutput.close();
    fs::copy_file(descriptor, directory / "flow.nek",
                  fs::copy_options::overwrite_existing);

    // Two logical frames deliberately share their physical time.  The
    // descriptor source index must keep both frames distinct.
    WriteNekField(directory / "flow0.f00000", 0, true, 0.25);
    WriteNekField(directory / "flow0.f00001", 1, false, 0.25);
    return descriptor;
}

fs::path CreateCurvedFixture(const fs::path& root) {
    const fs::path directory = root / "flowdata-nek5000-curved";
    fs::create_directories(directory);
    const fs::path descriptor = directory / "curved.nek5000";
    std::ofstream descriptorOutput(descriptor);
    Require(static_cast<bool>(descriptorOutput), "Cannot create curved Nek5000 descriptor.");
    descriptorOutput << "filetemplate: curved%01d.f%05d\n"
                     << "firsttimestep: 0\n"
                     << "numtimesteps: 1\n";
    descriptorOutput.close();
    WriteCurvedNekField(directory / "curved0.f00000");
    return descriptor;
}

void ValidateCanonicalFlowDataScan(const fs::path& descriptor) {
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(descriptor);
    Require(frames.size() == 2 && frames[0].time == frames[1].time,
            "Synthetic Nek5000 descriptor did not retain repeated logical frames.");

    io::ReadOptions options;
    vtkSmartPointer<vtkDataSet> canonical = io::ReadDataSet(frames.at(1), options);
    Require(canonical != nullptr && io::HasNek5000SpectralMetadata(canonical)
                && canonical->GetNumberOfPoints() == 8
                && canonical->GetNumberOfCells() == 0,
            "Synthetic Native Nek5000 frame was not preserved as cell-free canonical GLL data.");

    const QString descriptorPath = QString::fromStdWString(descriptor.wstring());
    const FlowDataScanResult direct = FlowDataInspector::scanNek5000File(
        descriptorPath, FlowDataKind::Raw, true);
    Require(direct.passed && direct.entry.fileType == QStringLiteral("nek5000 series")
                && direct.entry.meshCount == QStringLiteral("2")
                && direct.entry.dt.isEmpty()
                && direct.termCandidates.contains(QStringLiteral("Velocity"))
                && direct.termCandidates.contains(QStringLiteral("Pressure")),
            "Direct FlowData Nek5000 scan did not retain canonical frames and fields.");

    const FlowDataScanResult genericFile = FlowDataInspector::scanFolder(
        descriptorPath, FlowDataKind::Raw, true);
    Require(genericFile.passed && genericFile.entry.meshCount == QStringLiteral("2"),
            "Generic single-file FlowData scan bypassed the canonical Nek5000 path.");

    const QString directoryPath = QString::fromStdWString(descriptor.parent_path().wstring());
    const FlowDataScanResult genericFolder = FlowDataInspector::scanFolder(
        directoryPath, FlowDataKind::Raw, true);
    Require(genericFolder.passed && genericFolder.entry.meshCount == QStringLiteral("2"),
            "FlowData folder scan did not collapse .nek/.nek5000 aliases into one native series.");

    fs::path alias = descriptor;
    alias.replace_extension(".nek");
    const FlowDataScanResult aliasScan = FlowDataInspector::scanNek5000File(
        QString::fromStdWString(alias.wstring()), FlowDataKind::Raw, false);
    Require(aliasScan.passed && aliasScan.entry.meshCount == QStringLiteral("2"),
            "The .nek descriptor alias did not use the canonical FlowData reader path.");
}

void ValidateSpectralBounds(const fs::path& descriptor) {
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(descriptor);
    io::ReadOptions options;
    vtkSmartPointer<vtkDataSet> canonical = io::ReadDataSet(frames.front(), options);
    Require(canonical != nullptr && canonical->GetNumberOfCells() == 0,
            "Curved Native Nek5000 test frame is not canonical.");
    double nodalBounds[6] = {};
    canonical->GetBounds(nodalBounds);
    Require(nodalBounds[1] < 1.9,
            "Curved fixture no longer distinguishes nodal and spectral bounds.");

    const FlowDataScanResult scan = FlowDataInspector::scanNek5000File(
        QString::fromStdWString(descriptor.wstring()), FlowDataKind::Raw, false);
    const double spectralUpper = rangeUpperBound(scan.entry.xRange);
    Require(scan.passed && !scan.warningMessage.isEmpty()
                && spectralUpper >= 2.0 && spectralUpper > nodalBounds[1] + 0.1,
            "FlowData did not conservatively enclose canonical spectral geometry for Nek5000 bounds: "
                + scan.entry.xRange.toStdString());
}

void WriteStructuredSnapshot(const fs::path& path, double interiorXOffset) {
    vtkNew<vtkStructuredGrid> grid;
    grid->SetDimensions(3, 3, 3);
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(27);
    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(27);
    vtkIdType pointId = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i, ++pointId) {
                const double x = static_cast<double>(i) +
                    (i == 1 && j == 1 && k == 1 ? interiorXOffset : 0.0);
                points->SetPoint(pointId, x, static_cast<double>(j), static_cast<double>(k));
                velocity->SetTuple3(pointId, x, static_cast<double>(j), static_cast<double>(k));
            }
        }
    }
    grid->SetPoints(points);
    grid->GetPointData()->AddArray(velocity);
    vtkNew<vtkXMLStructuredGridWriter> writer;
    writer->SetFileName(path.string().c_str());
    writer->SetInputData(grid);
    Require(writer->Write() == 1,
            "Cannot write FlowData coordinate-validation fixture: " + path.string());
}

void ValidateCoordinateValidation(const fs::path& root) {
    const fs::path directory = root / "flowdata-coordinate-validation";
    fs::create_directories(directory);
    WriteStructuredSnapshot(directory / "snapshot0.vts", 0.0);
    WriteStructuredSnapshot(directory / "snapshot1.vts", 0.25);
    const FlowDataScanResult checked = FlowDataInspector::scanFolder(
        QString::fromStdWString(directory.wstring()), FlowDataKind::Raw, true);
    Require(!checked.passed && checked.errorMessage.contains(QStringLiteral("grid coordinates differ")),
            "Coordinate validation accepted an interior-node displacement with unchanged bounds.");
    const FlowDataScanResult unchecked = FlowDataInspector::scanFolder(
        QString::fromStdWString(directory.wstring()), FlowDataKind::Raw, false);
    Require(unchecked.passed,
            "Disabling coordinate validation unexpectedly rejected the valid field schema.");
}

}  // namespace

int main() {
    try {
        const fs::path root = fs::path(SCP_RESAMPLER_IO_TEST_OUTPUT_ROOT) /
                              "flowdata-nek";
        fs::remove_all(root);
        fs::create_directories(root);
        ValidateCanonicalFlowDataScan(CreateFixture(root));
        ValidateSpectralBounds(CreateCurvedFixture(root));
        ValidateCoordinateValidation(root);
        std::cout << "FlowData Native Nek5000 canonical scan test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FlowData Native Nek5000 canonical scan test failed: "
                  << error.what() << '\n';
        return 1;
    }
}
