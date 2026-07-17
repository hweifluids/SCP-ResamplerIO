#include <streamcenterplus/io/DataSetReader.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

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
#include <sstream>
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

bool NearlyEqual(double lhs, double rhs, double tolerance = 1.0e-10) {
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
void WriteValue(std::ofstream& output, T value, bool littleEndian) {
    if (littleEndian != HostIsLittleEndian()) {
        value = ByteSwap(value);
    }
    output.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

std::vector<double> Gll4() {
    const double inner = 1.0 / std::sqrt(5.0);
    return {-1.0, -inner, inner, 1.0};
}

std::vector<double> Gll8() {
    return {-1.0,
            -0.8717401485096066,
            -0.5917001814331423,
            -0.2092992179024789,
            0.2092992179024789,
            0.5917001814331423,
            0.8717401485096066,
            1.0};
}

std::size_t Index(int i, int j, int k, int nx, int ny) {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nx) +
           static_cast<std::size_t>(i);
}

void WriteNekField(const fs::path& path,
                   int frame,
                   bool coordinates,
                   bool littleEndian,
                   int elementId = 42,
                   int fileIndex = 0,
                   int fileCount = 1,
                   int globalElementCount = 1,
                   double xShift = 0.0,
                   int wordBytes = 8) {
    const int n = 4;
    const std::vector<double> nodes = Gll4();
    const std::string code = coordinates ? "XUP" : "UP";
    std::ostringstream header;
    header << "#std " << wordBytes << ' ' << n << ' ' << n << ' ' << n
           << " 1 " << globalElementCount << ' '
           << (0.25 * static_cast<double>(frame))
           << ' ' << frame << ' ' << fileIndex << ' ' << fileCount << ' '
           << code << " 0.0 F";
    std::string bytes = header.str();
    Require(bytes.size() <= 132, "Synthetic Nek header is too long.");
    bytes.resize(132, ' ');

    std::ofstream output(path, std::ios::binary);
    Require(static_cast<bool>(output), "Cannot create synthetic Nek field.");
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    WriteValue<float>(output, 6.54321f, littleEndian);
    WriteValue<std::int32_t>(output, elementId, littleEndian);

    const std::size_t pointCount = static_cast<std::size_t>(n * n * n);
    std::vector<double> x(pointCount, 0.0);
    std::vector<double> y(pointCount, 0.0);
    std::vector<double> z(pointCount, 0.0);
    std::vector<double> u(pointCount, 0.0);
    std::vector<double> v(pointCount, 0.0);
    std::vector<double> w(pointCount, 0.0);
    std::vector<double> p(pointCount, 0.0);
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                const std::size_t index = Index(i, j, k, n, n);
                x[index] = 1.0 + xShift + nodes[static_cast<std::size_t>(i)];
                y[index] = 2.0 + nodes[static_cast<std::size_t>(j)];
                z[index] = 3.0 + nodes[static_cast<std::size_t>(k)];
                const double offset = static_cast<double>(frame);
                u[index] = x[index] * x[index] + 2.0 * y[index] +
                           3.0 * z[index] + offset;
                v[index] = x[index] * y[index] - offset;
                w[index] = z[index] * z[index] * z[index] + 0.5 * offset;
                p[index] = x[index] + y[index] + z[index] + offset;
            }
        }
    }
    const auto writeArray = [&](const std::vector<double>& values) {
        for (double value : values) {
            if (wordBytes == 4) {
                WriteValue<float>(output, static_cast<float>(value), littleEndian);
            } else {
                WriteValue<double>(output, value, littleEndian);
            }
        }
    };
    if (coordinates) {
        writeArray(x);
        writeArray(y);
        writeArray(z);
    }
    writeArray(u);
    writeArray(v);
    writeArray(w);
    writeArray(p);
    Require(static_cast<bool>(output), "Failed while writing synthetic Nek field.");
}

fs::path CreateNekFloatFixture(const fs::path& root) {
    const fs::path directory = root / "nek5000-float";
    fs::create_directories(directory);
    {
        std::ofstream descriptor(directory / "float.nek5000");
        descriptor << "filetemplate: float%01d.f%05d\n"
                   << "firsttimestep: 0\n"
                   << "numtimesteps: 1\n";
    }
    WriteNekField(directory / "float0.f00000", 0, true, true,
                  11, 0, 1, 1, 0.0, 4);
    return directory / "float.nek5000";
}

void ValidateNekFloat(const fs::path& descriptor) {
    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(descriptor);
    const io::Nek5000SpectralMetadata metadata =
        io::GetNek5000SpectralMetadata(data);
    vtkDataArray* velocity = data->GetPointData()->GetArray("Velocity");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    Require(metadata.elementIds == std::vector<int>({11}) &&
                data->GetNumberOfPoints() == 64 && velocity != nullptr &&
                pressure != nullptr,
            "Four-byte little-endian Nek5000 metadata or arrays are incorrect.");
    double point[3] = {};
    data->GetPoint(37, point);
    Require(NearlyEqual(velocity->GetComponent(37, 0),
                        point[0] * point[0] + 2.0 * point[1] + 3.0 * point[2],
                        2.0e-5) &&
                NearlyEqual(pressure->GetComponent(37, 0),
                            point[0] + point[1] + point[2], 2.0e-5),
            "Four-byte little-endian Nek5000 values were not decoded correctly.");
}

fs::path CreateNekFixture(const fs::path& root) {
    const fs::path directory = root / "nek5000";
    fs::create_directories(directory);
    std::ofstream descriptor(directory / "tiny.nek5000");
    descriptor << "filetemplate: tiny%01d.f%05d\n"
               << "firsttimestep: 0\n"
               << "numtimesteps: 2\n";
    descriptor.close();
    fs::copy_file(directory / "tiny.nek5000", directory / "tiny.nek",
                  fs::copy_options::overwrite_existing);
    WriteNekField(directory / "tiny0.f00000", 0, true, false);
    WriteNekField(directory / "tiny0.f00001", 1, false, false);
    return directory / "tiny.nek5000";
}

void ValidateNek(const fs::path& descriptor) {
    const std::vector<io::DataSetFrame> frames =
        io::EnumerateDataSetFrames(descriptor);
    Require(frames.size() == 2, "Synthetic Nek descriptor must expose two frames.");
    fs::path alias = descriptor;
    alias.replace_extension(".nek");
    Require(io::EnumerateDataSetFrames(alias).size() == frames.size(),
            "The .nek Nek5000 descriptor alias was not recognized.");
    const std::vector<io::DataSetFrame> folderFrames =
        io::EnumerateDataSetFrames(descriptor.parent_path());
    Require(folderFrames.size() == frames.size(),
            "Folder enumeration duplicated same-stem .nek/.nek5000 aliases.");
    Require(frames[1].kind == io::DataSetFrameKind::Nek5000TimeStep,
            "Synthetic Nek frame kind is incorrect.");
    Require(NearlyEqual(frames[1].time, 0.25),
            "Synthetic Nek physical time was not preserved.");
    Require(frames[0].coordinateSourceIndex == 0 &&
                frames[1].coordinateSourceIndex == 0 &&
                frames[1].coordinatePath.filename() == "tiny0.f00000",
            "Nek enumeration did not cache the mesh-bearing coordinate frame.");

    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames[1]);
    Require(data != nullptr && data->GetNumberOfPoints() == 64 &&
                data->GetNumberOfCells() == 27,
            "Synthetic Nek topology has unexpected point/cell counts.");
    Require(io::HasNek5000SpectralMetadata(data),
            "Synthetic Nek dataset is missing spectral metadata.");
    const io::Nek5000SpectralMetadata metadata =
        io::GetNek5000SpectralMetadata(data);
    Require(metadata.nx == 4 && metadata.ny == 4 && metadata.nz == 4 &&
                metadata.elementCount == 1 && metadata.elementIds == std::vector<int>{42},
            "Synthetic Nek GLL metadata is incorrect.");

    vtkDataArray* velocity = data->GetPointData()->GetArray("Velocity");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    vtkDataArray* weights =
        data->GetPointData()->GetArray(io::kNek5000QuadratureWeightArray);
    Require(velocity != nullptr && velocity->GetDataType() == VTK_DOUBLE &&
                velocity->GetNumberOfComponents() == 3,
            "Synthetic Nek velocity was not retained as double3 point data.");
    Require(data->GetPointData()->GetArray("Velocity Magnitude") == nullptr,
            "The Nek reader generated a derived velocity magnitude array.");
    Require(pressure != nullptr && weights != nullptr,
            "Synthetic Nek pressure or quadrature weights are missing.");

    double weightSum = 0.0;
    for (vtkIdType i = 0; i < weights->GetNumberOfTuples(); ++i) {
        weightSum += weights->GetComponent(i, 0);
    }
    Require(NearlyEqual(weightSum, 8.0, 2.0e-12),
            "Nek GLL/Jacobian quadrature does not integrate the affine cube volume.");

    vtkSmartPointer<vtkDoubleArray> gradient =
        io::ComputeNek5000SpectralGradient(data, pressure, "pressure_gradient");
    Require(gradient != nullptr && gradient->GetNumberOfComponents() == 3,
            "Nek spectral pressure gradient has an unexpected layout.");
    for (vtkIdType i = 0; i < gradient->GetNumberOfTuples(); ++i) {
        Require(NearlyEqual(gradient->GetComponent(i, 0), 1.0, 2.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 1), 1.0, 2.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 2), 1.0, 2.0e-11),
                "Nek spectral gradient is not exact for a linear physical field.");
    }

    double firstPoint[3] = {};
    data->GetPoint(0, firstPoint);
    Require(NearlyEqual(firstPoint[0], 0.0) && NearlyEqual(firstPoint[1], 1.0) &&
                NearlyEqual(firstPoint[2], 2.0),
            "Nek frame without X did not reuse the mesh-bearing first frame.");
    Require(NearlyEqual(pressure->GetComponent(0, 0), 4.0),
            "Nek physical field value was not read at its original GLL point.");

    io::Nek5000SpectralPointEvaluator evaluator(data, velocity);
    const double query[3] = {0.37, 1.72, 2.44};
    double interpolated[3] = {};
    int elementHint = -1;
    Require(evaluator.Evaluate(query, interpolated, &elementHint) &&
                elementHint == 0,
            "Nek element-local spectral point location failed inside the cube.");
    Require(NearlyEqual(interpolated[0],
                        query[0] * query[0] + 2.0 * query[1] +
                            3.0 * query[2] + 1.0,
                        2.0e-11) &&
                NearlyEqual(interpolated[1], query[0] * query[1] - 1.0,
                            2.0e-11) &&
                NearlyEqual(interpolated[2],
                            query[2] * query[2] * query[2] + 0.5,
                            2.0e-11),
            "Nek GLL spectral interpolation is not exact for represented polynomials.");

    io::Nek5000SpectralPointEvaluator pressureEvaluator =
        evaluator.WithPointArray(pressure);
    double interpolatedPressure = 0.0;
    Require(pressureEvaluator.numberOfComponents() == 1 &&
                pressureEvaluator.Evaluate(query, &interpolatedPressure,
                                           &elementHint) &&
                NearlyEqual(interpolatedPressure,
                            query[0] + query[1] + query[2] + 1.0,
                            2.0e-11),
            "A Nek evaluator rebound to another GLL field did not reuse the "
            "spectral geometry correctly.");
    const double outside[3] = {3.0, 1.72, 2.44};
    Require(!evaluator.Evaluate(outside, interpolated, &elementHint),
            "Nek spectral point evaluator accepted a point outside the mesh.");

    io::ReadOptions noTopologyOptions;
    noTopologyOptions.buildNek5000LinearizedTopology = false;
    vtkSmartPointer<vtkDataSet> noTopology =
        io::ReadDataSet(frames[1], noTopologyOptions);
    Require(noTopology->GetNumberOfPoints() == 64 &&
                noTopology->GetNumberOfCells() == 0 &&
                io::HasNek5000SpectralMetadata(noTopology) &&
                noTopology->GetFieldData()->GetArray(
                    io::kNek5000ElementIdsField) != nullptr &&
                noTopology->GetFieldData()->GetAbstractArray(
                    io::kNek5000RepresentationField) != nullptr &&
                noTopology->GetPointData()->GetArray(
                    io::kNek5000QuadratureWeightArray) != nullptr,
            "Disabling Nek linearized topology did not preserve points, basic "
            "spectral metadata and optional quadrature arrays.");
    vtkDataArray* noTopologyVelocity =
        noTopology->GetPointData()->GetArray("Velocity");
    io::Nek5000SpectralPointEvaluator noTopologyEvaluator(
        noTopology, noTopologyVelocity);
    elementHint = -1;
    Require(noTopologyEvaluator.Evaluate(query, interpolated, &elementHint) &&
                NearlyEqual(interpolated[0],
                            query[0] * query[0] + 2.0 * query[1] +
                                3.0 * query[2] + 1.0,
                            2.0e-11) &&
                !noTopologyEvaluator.Evaluate(outside, interpolated,
                                              &elementHint),
            "The element-AABB locator failed without linearized Nek cells.");

    io::ReadOptions noQuadratureOptions;
    noQuadratureOptions.attachNek5000QuadratureArrays = false;
    vtkSmartPointer<vtkDataSet> noQuadrature =
        io::ReadDataSet(frames[1], noQuadratureOptions);
    Require(noQuadrature->GetNumberOfCells() == 27 &&
                io::HasNek5000SpectralMetadata(noQuadrature) &&
                noQuadrature->GetFieldData()->GetArray(
                    io::kNek5000ElementIdsField) != nullptr &&
                noQuadrature->GetPointData()->GetArray(
                    io::kNek5000ReferenceWeightArray) == nullptr &&
                noQuadrature->GetPointData()->GetArray(
                    io::kNek5000JacobianDeterminantArray) == nullptr &&
                noQuadrature->GetPointData()->GetArray(
                    io::kNek5000QuadratureWeightArray) == nullptr,
            "Disabling Nek quadrature arrays removed basic metadata or retained "
            "a per-point quadrature array.");

    io::ReadOptions physicalWeightOnlyOptions;
    physicalWeightOnlyOptions.buildNek5000LinearizedTopology = false;
    physicalWeightOnlyOptions.attachNek5000ReferenceWeightArray = false;
    physicalWeightOnlyOptions.attachNek5000JacobianDeterminantArray = false;
    vtkSmartPointer<vtkDataSet> physicalWeightOnly =
        io::ReadDataSet(frames[1], physicalWeightOnlyOptions);
    Require(physicalWeightOnly->GetNumberOfCells() == 0 &&
                physicalWeightOnly->GetPointData()->GetArray(
                    io::kNek5000ReferenceWeightArray) == nullptr &&
                physicalWeightOnly->GetPointData()->GetArray(
                    io::kNek5000JacobianDeterminantArray) == nullptr &&
                physicalWeightOnly->GetPointData()->GetArray(
                    io::kNek5000QuadratureWeightArray) != nullptr,
            "Nek physical-weight-only mode retained unused large arrays or "
            "removed the physical quadrature measure.");

    vtkSmartPointer<vtkUnstructuredGrid> geometryCopy =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    geometryCopy->DeepCopy(data);
    Require(io::HasSameNek5000SpectralGeometry(data, geometryCopy),
            "Identical Nek spectral geometries were reported as different.");
    double movedPoint[3] = {};
    geometryCopy->GetPoint(0, movedPoint);
    movedPoint[0] += 1.0e-6;
    geometryCopy->GetPoints()->SetPoint(0, movedPoint);
    Require(!io::HasSameNek5000SpectralGeometry(data, geometryCopy) &&
                io::HasSameNek5000SpectralGeometry(data, geometryCopy,
                                                   2.0e-6),
            "Nek spectral geometry comparison did not detect/tolerate a "
            "coordinate change as configured.");
}

fs::path CreateMultiFileNekFixture(const fs::path& root) {
    const fs::path directory = root / "nek5000-multifile";
    fs::create_directories(directory);
    {
        std::ofstream descriptor(directory / "multi.nek5000");
        descriptor << "filetemplate: multi%01d.f%05d\n"
                   << "firsttimestep: 0\n"
                   << "numtimesteps: 2\n";
    }
    // Frame zero deliberately stores the higher global id in piece zero.
    WriteNekField(directory / "multi0.f00000", 0, true, true,
                  200, 0, 2, 2, 2.0);
    WriteNekField(directory / "multi1.f00000", 0, true, true,
                  100, 1, 2, 2, 0.0);
    // Frame one reverses ownership/order and omits coordinates. Both frames
    // must still expose the same canonical global-element ordering.
    WriteNekField(directory / "multi0.f00001", 1, false, true,
                  100, 0, 2, 2, 0.0);
    WriteNekField(directory / "multi1.f00001", 1, false, true,
                  200, 1, 2, 2, 2.0);
    return directory / "multi.nek5000";
}

void ValidateMultiFileNek(const fs::path& descriptor) {
    const std::vector<io::DataSetFrame> frames =
        io::EnumerateDataSetFrames(descriptor);
    Require(frames.size() == 2,
            "Multi-file Nek fixture did not expose both frames.");
    vtkSmartPointer<vtkDataSet> first = io::ReadDataSet(frames[0]);
    vtkSmartPointer<vtkDataSet> second = io::ReadDataSet(frames[1]);
    for (vtkDataSet* data : {first.GetPointer(), second.GetPointer()}) {
        const io::Nek5000SpectralMetadata metadata =
            io::GetNek5000SpectralMetadata(data);
        Require(metadata.elementCount == 2 &&
                    metadata.elementIds == std::vector<int>({100, 200}) &&
                    data->GetNumberOfPoints() == 128 &&
                    data->GetNumberOfCells() == 54,
                "Multi-file Nek element maps were not canonicalized by global id.");
        vtkDataArray* weights = data->GetPointData()->GetArray(
            io::kNek5000QuadratureWeightArray);
        double sum = 0.0;
        for (vtkIdType i = 0; i < weights->GetNumberOfTuples(); ++i) {
            sum += weights->GetComponent(i, 0);
        }
        Require(NearlyEqual(sum, 16.0, 3.0e-12),
                "Multi-file Nek quadrature does not cover both elements.");
        double point[3] = {};
        data->GetPoint(0, point);
        Require(NearlyEqual(point[0], 0.0),
                "Canonical Nek element 100 was not placed first.");
        data->GetPoint(64, point);
        Require(NearlyEqual(point[0], 2.0),
                "Canonical Nek element 200 was not placed second.");
    }
    Require(io::HasSameNek5000SpectralGeometry(first, second),
            "Repartitioned multi-file Nek frames lost canonical spectral geometry.");
    vtkDataArray* secondPressure = second->GetPointData()->GetArray("Pressure");
    Require(secondPressure != nullptr &&
                NearlyEqual(secondPressure->GetComponent(0, 0), 4.0) &&
                NearlyEqual(secondPressure->GetComponent(64, 0), 6.0),
            "Repartitioned multi-file Nek fields were not scattered by global id.");
}

fs::path CreateNek2DFixture(const fs::path& root) {
    const fs::path directory = root / "nek5000-2d";
    fs::create_directories(directory);
    {
        std::ofstream descriptor(directory / "tiny2d.nek5000");
        descriptor << "filetemplate: tiny2d%01d.f%05d\n"
                   << "firsttimestep: 0\n"
                   << "numtimesteps: 1\n";
    }
    std::ostringstream header;
    header << "#std 8 4 3 1 1 1 0 0 0 1 XUP 0.0 F";
    std::string headerBytes = header.str();
    headerBytes.resize(132, ' ');
    std::ofstream output(directory / "tiny2d0.f00000", std::ios::binary);
    output.write(headerBytes.data(), static_cast<std::streamsize>(headerBytes.size()));
    WriteValue<float>(output, 6.54321f, false);
    WriteValue<std::int32_t>(output, 7, false);

    const std::array<double, 4> xNodes = {
        -1.0, -1.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0), 1.0};
    const std::array<double, 3> yNodes = {-1.0, 0.0, 1.0};
    std::array<std::vector<double>, 5> fields;
    for (auto& field : fields) {
        field.resize(12);
    }
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 4; ++i) {
            const std::size_t index = static_cast<std::size_t>(j * 4 + i);
            const double x = 1.0 + xNodes[static_cast<std::size_t>(i)];
            const double y = 2.0 + yNodes[static_cast<std::size_t>(j)];
            fields[0][index] = x;
            fields[1][index] = y;
            fields[2][index] = x * x * x + y * y;
            fields[3][index] = x * y;
            fields[4][index] = x * x + 2.0 * y;
        }
    }
    for (const auto& field : fields) {
        for (const double value : field) {
            WriteValue<double>(output, value, false);
        }
    }
    Require(static_cast<bool>(output), "Failed while writing synthetic 2-D Nek field.");
    return directory / "tiny2d.nek5000";
}

void ValidateNek2D(const fs::path& descriptor) {
    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(descriptor);
    const io::Nek5000SpectralMetadata metadata =
        io::GetNek5000SpectralMetadata(data);
    Require(metadata.nx == 4 && metadata.ny == 3 && metadata.nz == 1 &&
                data->GetNumberOfPoints() == 12 && data->GetNumberOfCells() == 6,
            "Synthetic 2-D Nek topology or tensor dimensions are incorrect.");
    vtkDataArray* velocity = data->GetPointData()->GetArray("Velocity");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    vtkDataArray* weights =
        data->GetPointData()->GetArray(io::kNek5000QuadratureWeightArray);
    Require(velocity != nullptr && pressure != nullptr && weights != nullptr,
            "Synthetic 2-D Nek arrays are missing.");
    double weightSum = 0.0;
    for (vtkIdType i = 0; i < weights->GetNumberOfTuples(); ++i) {
        weightSum += weights->GetComponent(i, 0);
        Require(NearlyEqual(velocity->GetComponent(i, 2), 0.0),
                "2-D Nek velocity did not receive a zero z component.");
    }
    Require(NearlyEqual(weightSum, 4.0, 2.0e-12),
            "2-D Nek surface GLL quadrature did not integrate the rectangle area.");

    vtkSmartPointer<vtkDoubleArray> gradient =
        io::ComputeNek5000SpectralGradient(data, pressure, "pressure_gradient_2d");
    double point[3] = {};
    for (vtkIdType i = 0; i < gradient->GetNumberOfTuples(); ++i) {
        data->GetPoint(i, point);
        Require(NearlyEqual(gradient->GetComponent(i, 0), 2.0 * point[0], 2.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 1), 2.0, 2.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 2), 0.0, 2.0e-11),
                "2-D Nek tangent-metric spectral gradient is incorrect.");
    }
    io::Nek5000SpectralPointEvaluator evaluator(data, velocity);
    const double query[3] = {0.37, 1.72, 0.0};
    double tuple[3] = {};
    Require(evaluator.Evaluate(query, tuple) &&
                NearlyEqual(tuple[0], query[0] * query[0] * query[0] +
                                          query[1] * query[1],
                            2.0e-11) &&
                NearlyEqual(tuple[1], query[0] * query[1], 2.0e-11),
            "2-D Nek tensor-product GLL interpolation is incorrect.");
}

void ValidateCurvedNekMetrics() {
    const std::vector<double> nodes = Gll4();
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(64);
    vtkNew<vtkDoubleArray> scalar;
    scalar->SetName("CurvedLinearPhysicalField");
    scalar->SetNumberOfTuples(64);
    vtkIdType pointId = 0;
    for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i, ++pointId) {
                const double r = nodes[static_cast<std::size_t>(i)];
                const double s = nodes[static_cast<std::size_t>(j)];
                const double t = nodes[static_cast<std::size_t>(k)];
                const double x = 1.0 + r + 0.10 * r * s + 0.03 * s * t;
                const double y = 2.0 + s + 0.08 * r * t;
                const double z = 3.0 + t + 0.05 * r * s;
                points->SetPoint(pointId, x, y, z);
                scalar->SetValue(pointId, 2.0 * x - 3.0 * y + 4.0 * z);
            }
        }
    }
    vtkSmartPointer<vtkUnstructuredGrid> data =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    data->SetPoints(points);
    data->GetPointData()->AddArray(scalar);
    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(io::kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(4, 4, 4);
    data->GetFieldData()->AddArray(dimensions);
    vtkNew<vtkIntArray> count;
    count->SetName(io::kNek5000ElementCountField);
    count->InsertNextValue(1);
    data->GetFieldData()->AddArray(count);

    vtkSmartPointer<vtkDoubleArray> gradient =
        io::ComputeNek5000SpectralGradient(data, scalar, "curved_gradient");
    for (vtkIdType i = 0; i < gradient->GetNumberOfTuples(); ++i) {
        Require(NearlyEqual(gradient->GetComponent(i, 0), 2.0, 5.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 1), -3.0, 5.0e-11) &&
                    NearlyEqual(gradient->GetComponent(i, 2), 4.0, 5.0e-11),
                "Curved Nek inverse-Jacobian spectral gradient is incorrect.");
    }

    const double r = 0.23;
    const double s = -0.31;
    const double t = 0.17;
    const double query[3] = {
        1.0 + r + 0.10 * r * s + 0.03 * s * t,
        2.0 + s + 0.08 * r * t,
        3.0 + t + 0.05 * r * s,
    };
    io::Nek5000SpectralPointEvaluator evaluator(data, scalar);
    double value = 0.0;
    Require(evaluator.Evaluate(query, &value) &&
                NearlyEqual(value,
                            2.0 * query[0] - 3.0 * query[1] + 4.0 * query[2],
                            5.0e-11),
            "Curved Nek reference-coordinate inversion/interpolation is incorrect.");
    const double outside[3] = {8.0, 8.0, 8.0};
    Require(!evaluator.Evaluate(outside, &value),
            "Cell-free curved Nek AABB locator accepted an exterior point.");
}

void ValidateConservativeCurvedNekBounds() {
    const std::vector<double> nodes = Gll8();
    const int n = static_cast<int>(nodes.size());
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(n * n * n);
    vtkNew<vtkDoubleArray> scalar;
    scalar->SetName("ReferenceLinearField");
    scalar->SetNumberOfTuples(n * n * n);
    double nodeMaximumX = -std::numeric_limits<double>::infinity();
    vtkIdType pointId = 0;
    for (int k = 0; k < n; ++k) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i, ++pointId) {
                const double r = nodes[static_cast<std::size_t>(i)];
                const double s = nodes[static_cast<std::size_t>(j)];
                const double t = nodes[static_cast<std::size_t>(k)];
                // This triangular map has determinant one and is globally
                // invertible, but x reaches its maximum at s=0, between the
                // two interior GLL nodes. A node-only AABB therefore misses a
                // valid part of the spectral element.
                const double x = r + 0.20 * (1.0 - s * s);
                const double y = s;
                const double z = t;
                points->SetPoint(pointId, x, y, z);
                scalar->SetValue(pointId, 2.0 * r - 3.0 * s + 4.0 * t);
                nodeMaximumX = std::max(nodeMaximumX, x);
            }
        }
    }
    vtkSmartPointer<vtkUnstructuredGrid> data =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    data->SetPoints(points);
    data->GetPointData()->AddArray(scalar);
    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(io::kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(n, n, n);
    data->GetFieldData()->AddArray(dimensions);
    vtkNew<vtkIntArray> count;
    count->SetName(io::kNek5000ElementCountField);
    count->InsertNextValue(1);
    data->GetFieldData()->AddArray(count);

    const double r = 0.999;
    const double s = 0.0;
    const double t = 0.13;
    const double query[3] = {r + 0.20 * (1.0 - s * s), s, t};
    Require(query[0] > nodeMaximumX,
            "Curved-bound regression point is not outside the GLL node AABB.");

    io::Nek5000SpectralPointEvaluator evaluator(data, scalar);
    const std::array<double, 6> bounds = evaluator.bounds();
    Require(bounds[0] <= query[0] && bounds[1] >= 1.20 &&
                bounds[2] <= query[1] && bounds[3] >= query[1] &&
                bounds[4] <= query[2] && bounds[5] >= query[2],
            "High-order Nek evaluator bounds do not enclose the curved map.");
    double value = 0.0;
    Require(evaluator.Evaluate(query, &value) &&
                NearlyEqual(value, 2.0 * r - 3.0 * s + 4.0 * t, 5.0e-11),
            "Conservative Nek AABB did not expose a valid point beyond the node AABB.");
}

void ValidateCellFreeElementBvh() {
    constexpr int elementCount = 10;
    constexpr int pointsPerElement = 8;
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(elementCount * pointsPerElement);
    vtkNew<vtkDoubleArray> scalar;
    scalar->SetName("LinearPhysicalField");
    scalar->SetNumberOfTuples(elementCount * pointsPerElement);
    for (int element = 0; element < elementCount; ++element) {
        const double center = 4.0 * static_cast<double>(element);
        for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
                for (int i = 0; i < 2; ++i) {
                    const vtkIdType pointId =
                        static_cast<vtkIdType>(element * pointsPerElement +
                                               k * 4 + j * 2 + i);
                    const double x = center + (i == 0 ? -1.0 : 1.0);
                    const double y = j == 0 ? -1.0 : 1.0;
                    const double z = k == 0 ? -1.0 : 1.0;
                    points->SetPoint(pointId, x, y, z);
                    scalar->SetValue(pointId, x + 2.0 * y + 3.0 * z);
                }
            }
        }
    }
    vtkSmartPointer<vtkUnstructuredGrid> data =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    data->SetPoints(points);
    data->GetPointData()->AddArray(scalar);
    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(io::kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(2, 2, 2);
    data->GetFieldData()->AddArray(dimensions);
    vtkNew<vtkIntArray> count;
    count->SetName(io::kNek5000ElementCountField);
    count->InsertNextValue(elementCount);
    data->GetFieldData()->AddArray(count);
    vtkNew<vtkIntArray> ids;
    ids->SetName(io::kNek5000ElementIdsField);
    for (int element = 0; element < elementCount; ++element) {
        ids->InsertNextValue(100 + element);
    }
    data->GetFieldData()->AddArray(ids);

    io::Nek5000SpectralPointEvaluator evaluator(data, scalar);
    const double query[3] = {36.2, 0.1, -0.3};
    double value = 0.0;
    int hint = -1;
    Require(evaluator.Evaluate(query, &value, &hint) && hint == 9 &&
                NearlyEqual(value,
                            query[0] + 2.0 * query[1] + 3.0 * query[2],
                            2.0e-12),
            "Cell-free Nek element-AABB BVH did not locate the expected leaf.");
    const double outside[3] = {100.0, 0.0, 0.0};
    Require(!evaluator.Evaluate(outside, &value, &hint),
            "Cell-free Nek element-AABB BVH accepted an exterior point.");
}

fs::path CreateEnSightFixture(const fs::path& root) {
    const fs::path directory = root / "ensight";
    fs::create_directories(directory);
    {
        std::ofstream file(directory / "tiny.case");
        file << "FORMAT\n"
             << "type: ensight gold\n\n"
             << "GEOMETRY\n"
             << "model: tiny.geo\n\n"
             << "VARIABLE\n"
             << "vector per node: Velocity tiny.vel\n"
             << "scalar per node: Pressure tiny.pre\n";
    }
    const std::array<std::array<double, 3>, 8> points = {{{0, 0, 0}, {1, 0, 0},
                                                           {1, 1, 0}, {0, 1, 0},
                                                           {0, 0, 1}, {1, 0, 1},
                                                           {1, 1, 1}, {0, 1, 1}}};
    {
        std::ofstream file(directory / "tiny.geo");
        file << "EnSight Gold geometry\nGenerated by SCP test\n"
             << "node id assign\nelement id assign\n"
             << "part\n1\nvolume\ncoordinates\n8\n";
        for (int component = 0; component < 3; ++component) {
            for (const auto& point : points) {
                file << point[static_cast<std::size_t>(component)] << '\n';
            }
        }
        file << "hexa8\n1\n1 2 3 4 5 6 7 8\n";
    }
    {
        std::ofstream file(directory / "tiny.vel");
        file << "Velocity\npart\n1\ncoordinates\n";
        for (int component = 0; component < 3; ++component) {
            for (const auto& point : points) {
                file << point[static_cast<std::size_t>(component)] + 1.0 << '\n';
            }
        }
    }
    {
        std::ofstream file(directory / "tiny.pre");
        file << "Pressure\npart\n1\ncoordinates\n";
        for (const auto& point : points) {
            file << point[0] + point[1] + point[2] << '\n';
        }
    }
    return directory / "tiny.case";
}

void ValidateEnSight(const fs::path& casePath) {
    const std::vector<io::DataSetFrame> frames =
        io::EnumerateDataSetFrames(casePath);
    Require(frames.size() == 1 &&
                frames.front().kind == io::DataSetFrameKind::EnSightGoldTimeStep,
            "Static EnSight Gold fixture must expose one native frame.");
    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfPoints() == 8 &&
                data->GetNumberOfCells() == 1,
            "EnSight Gold fixture topology is incorrect.");
    Require(data->GetPointData()->GetArray("Velocity") != nullptr &&
                data->GetPointData()->GetArray("Pressure") != nullptr,
            "EnSight Gold point variables were not read.");
}

void WriteFoamHeader(std::ofstream& file,
                     const std::string& objectClass,
                     const std::string& location,
                     const std::string& object) {
    file << "FoamFile\n{\n"
         << "    version 2.0;\n    format ascii;\n"
         << "    class " << objectClass << ";\n"
         << "    location \"" << location << "\";\n"
         << "    object " << object << ";\n}\n";
}

fs::path CreateOpenFoamFixture(const fs::path& root) {
    const fs::path directory = root / "openfoam";
    fs::create_directories(directory / "system");
    fs::create_directories(directory / "constant" / "polyMesh");
    fs::create_directories(directory / "0");
    std::ofstream(directory / "tiny.foam").close();
    {
        std::ofstream file(directory / "system" / "controlDict");
        WriteFoamHeader(file, "dictionary", "system", "controlDict");
        file << "application simpleFoam;\nstartFrom startTime;\nstartTime 0;\n"
             << "stopAt endTime;\nendTime 0;\ndeltaT 1;\n"
             << "writeControl timeStep;\nwriteInterval 1;\n";
    }
    {
        std::ofstream file(directory / "constant" / "polyMesh" / "points");
        WriteFoamHeader(file, "vectorField", "constant/polyMesh", "points");
        file << "8\n(\n(0 0 0)\n(1 0 0)\n(1 1 0)\n(0 1 0)\n"
             << "(0 0 1)\n(1 0 1)\n(1 1 1)\n(0 1 1)\n)\n";
    }
    {
        std::ofstream file(directory / "constant" / "polyMesh" / "faces");
        WriteFoamHeader(file, "faceList", "constant/polyMesh", "faces");
        file << "6\n(\n4(0 3 2 1)\n4(4 5 6 7)\n4(0 1 5 4)\n"
             << "4(1 2 6 5)\n4(2 3 7 6)\n4(3 0 4 7)\n)\n";
    }
    {
        std::ofstream file(directory / "constant" / "polyMesh" / "owner");
        WriteFoamHeader(file, "labelList", "constant/polyMesh", "owner");
        file << "6\n(\n0\n0\n0\n0\n0\n0\n)\n";
    }
    {
        std::ofstream file(directory / "constant" / "polyMesh" / "neighbour");
        WriteFoamHeader(file, "labelList", "constant/polyMesh", "neighbour");
        file << "0\n(\n)\n";
    }
    {
        std::ofstream file(directory / "constant" / "polyMesh" / "boundary");
        WriteFoamHeader(file, "polyBoundaryMesh", "constant/polyMesh", "boundary");
        file << "6\n(\n"
             << "bottom { type wall; nFaces 1; startFace 0; }\n"
             << "top { type wall; nFaces 1; startFace 1; }\n"
             << "front { type wall; nFaces 1; startFace 2; }\n"
             << "right { type wall; nFaces 1; startFace 3; }\n"
             << "back { type wall; nFaces 1; startFace 4; }\n"
             << "left { type wall; nFaces 1; startFace 5; }\n)\n";
    }
    const std::string boundaryNames = "bottom top front right back left";
    {
        std::ofstream file(directory / "0" / "U");
        WriteFoamHeader(file, "volVectorField", "0", "U");
        file << "dimensions [0 1 -1 0 0 0 0];\ninternalField uniform (1 2 3);\n"
             << "boundaryField\n{\n";
        std::istringstream names(boundaryNames);
        std::string name;
        while (names >> name) {
            file << name << " { type zeroGradient; }\n";
        }
        file << "}\n";
    }
    {
        std::ofstream file(directory / "0" / "p");
        WriteFoamHeader(file, "volScalarField", "0", "p");
        file << "dimensions [0 2 -2 0 0 0 0];\ninternalField uniform 7;\n"
             << "boundaryField\n{\n";
        std::istringstream names(boundaryNames);
        std::string name;
        while (names >> name) {
            file << name << " { type zeroGradient; }\n";
        }
        file << "}\n";
    }
    return directory / "tiny.foam";
}

void ValidateOpenFoam(const fs::path& marker) {
    const std::vector<io::DataSetFrame> frames =
        io::EnumerateDataSetFrames(marker);
    Require(!frames.empty() &&
                frames.front().kind == io::DataSetFrameKind::OpenFoamTimeStep,
            "OpenFOAM fixture did not expose a native time frame.");
    vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
    Require(data != nullptr && data->GetNumberOfCells() == 1 &&
                data->GetNumberOfPoints() == 8,
            "OpenFOAM internal mesh topology is incorrect.");
    Require(data->GetCellData()->GetArray("U") != nullptr &&
                data->GetCellData()->GetArray("p") != nullptr,
            "OpenFOAM cell fields U/p were not preserved.");
}

}  // namespace

int main() {
    try {
        const fs::path root =
            fs::path(SCP_REPOSITORY_ROOT) / ".tests" / "native-readers" /
            "SCP-ResamplerIO";
        fs::remove_all(root);
        fs::create_directories(root);
        ValidateNek(CreateNekFixture(root));
        ValidateNekFloat(CreateNekFloatFixture(root));
        ValidateMultiFileNek(CreateMultiFileNekFixture(root));
        ValidateNek2D(CreateNek2DFixture(root));
        ValidateCurvedNekMetrics();
        ValidateConservativeCurvedNekBounds();
        ValidateCellFreeElementBvh();
        ValidateEnSight(CreateEnSightFixture(root));
        ValidateOpenFoam(CreateOpenFoamFixture(root));
        std::cout << "Native reader tests passed: Nek5000, EnSight Gold, OpenFOAM.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Native reader test failed: " << error.what() << '\n';
        return 1;
    }
}
