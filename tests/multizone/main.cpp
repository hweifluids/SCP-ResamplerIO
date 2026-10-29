#include <streamcenterplus/io/DataSetReader.h>
#include <streamcenterplus_resample/Core.h>

#include <vtk_hdf5.h>

#include <vtkAppendFilter.h>
#include <vtkDataArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkGenericCell.h>
#include <vtkHexahedron.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkStaticCellLocator.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVoxel.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace io = streamcenterplus::io;
namespace resample = streamcenterplus::resample;

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

bool NearlyEqual(double lhs, double rhs, double tolerance = 1.0e-5) {
    return std::abs(lhs - rhs) <= tolerance *
           std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

vtkIdType GridPointId(int i, int j, int k, int divisions) {
    const vtkIdType width = static_cast<vtkIdType>(divisions + 1);
    return (static_cast<vtkIdType>(k) * width + static_cast<vtkIdType>(j)) * width +
           static_cast<vtkIdType>(i);
}

vtkSmartPointer<vtkUnstructuredGrid> CreateBoxGrid(double minimum,
                                                    double maximum,
                                                    int divisions,
                                                    bool scalarVelocityComponents,
                                                    double pressure,
                                                    double velocityX,
                                                    bool scalarVelocityOnCells = false,
                                                    bool firstCellVoxel = false) {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    const vtkIdType pointCount = static_cast<vtkIdType>(divisions + 1) *
                                 static_cast<vtkIdType>(divisions + 1) *
                                 static_cast<vtkIdType>(divisions + 1);
    const vtkIdType cellCount = static_cast<vtkIdType>(divisions) *
                                static_cast<vtkIdType>(divisions) *
                                static_cast<vtkIdType>(divisions);
    points->SetNumberOfPoints(pointCount);
    vtkNew<vtkDoubleArray> pressureArray;
    pressureArray->SetName("Pressure");
    pressureArray->SetNumberOfComponents(1);
    pressureArray->SetNumberOfTuples(pointCount);
    vtkNew<vtkDoubleArray> velocity;
    vtkNew<vtkDoubleArray> u;
    vtkNew<vtkDoubleArray> v;
    vtkNew<vtkDoubleArray> w;
    if (scalarVelocityComponents) {
        u->SetName("uLocal");
        v->SetName("vLocal");
        w->SetName("wLocal");
        for (vtkDoubleArray* array : {u.GetPointer(), v.GetPointer(), w.GetPointer()}) {
            array->SetNumberOfComponents(1);
            array->SetNumberOfTuples(scalarVelocityOnCells ? cellCount : pointCount);
        }
    } else {
        velocity->SetName("Velocity");
        velocity->SetNumberOfComponents(3);
        velocity->SetNumberOfTuples(pointCount);
    }

    vtkIdType pointId = 0;
    for (int k = 0; k <= divisions; ++k) {
        for (int j = 0; j <= divisions; ++j) {
            for (int i = 0; i <= divisions; ++i, ++pointId) {
                const auto coordinate = [&](int index) {
                    return minimum + (maximum - minimum) *
                                         static_cast<double>(index) /
                                         static_cast<double>(divisions);
                };
                points->SetPoint(pointId, coordinate(i), coordinate(j), coordinate(k));
                pressureArray->SetValue(pointId, pressure);
                if (scalarVelocityComponents && !scalarVelocityOnCells) {
                    u->SetValue(pointId, velocityX);
                    v->SetValue(pointId, 0.0);
                    w->SetValue(pointId, 0.0);
                } else if (!scalarVelocityComponents) {
                    velocity->SetTuple3(pointId, velocityX, 0.0, 0.0);
                }
            }
        }
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->GetPointData()->AddArray(pressureArray);
    if (!scalarVelocityComponents) {
        grid->GetPointData()->AddArray(velocity);
    }

    for (int k = 0; k < divisions; ++k) {
        for (int j = 0; j < divisions; ++j) {
            for (int i = 0; i < divisions; ++i) {
                const std::array<vtkIdType, 8> ids = {
                    GridPointId(i, j, k, divisions),
                    GridPointId(i + 1, j, k, divisions),
                    GridPointId(i + 1, j + 1, k, divisions),
                    GridPointId(i, j + 1, k, divisions),
                    GridPointId(i, j, k + 1, divisions),
                    GridPointId(i + 1, j, k + 1, divisions),
                    GridPointId(i + 1, j + 1, k + 1, divisions),
                    GridPointId(i, j + 1, k + 1, divisions),
                };
                if (firstCellVoxel && i == 0 && j == 0 && k == 0) {
                    vtkNew<vtkVoxel> cell;
                    for (int point = 0; point < 8; ++point) {
                        cell->GetPointIds()->SetId(
                            point, ids[static_cast<std::size_t>(point)]);
                    }
                    grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
                } else {
                    vtkNew<vtkHexahedron> cell;
                    for (int point = 0; point < 8; ++point) {
                        cell->GetPointIds()->SetId(
                            point, ids[static_cast<std::size_t>(point)]);
                    }
                    grid->InsertNextCell(cell->GetCellType(), cell->GetPointIds());
                }
            }
        }
    }
    if (scalarVelocityComponents) {
        if (scalarVelocityOnCells) {
            for (vtkIdType cellId = 0; cellId < cellCount; ++cellId) {
                u->SetValue(cellId, velocityX);
                v->SetValue(cellId, 0.0);
                w->SetValue(cellId, 0.0);
            }
            grid->GetCellData()->AddArray(u);
            grid->GetCellData()->AddArray(v);
            grid->GetCellData()->AddArray(w);
        } else {
            grid->GetPointData()->AddArray(u);
            grid->GetPointData()->AddArray(v);
            grid->GetPointData()->AddArray(w);
        }
    }
    return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> CreateCellVelocityBoxGrid(double minimum,
                                                                double maximum,
                                                                int divisions,
                                                                double pressure,
                                                                double velocityX,
                                                                double cellMarker,
                                                                const char* velocityName = "Velocity",
                                                                bool addStaleVelocityAlias = false) {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        CreateBoxGrid(minimum, maximum, divisions, false, pressure, velocityX);
    grid->GetPointData()->RemoveArray("Velocity");

    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName(velocityName);
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(grid->GetNumberOfCells());
    vtkNew<vtkDoubleArray> marker;
    marker->SetName("CellMarker");
    marker->SetNumberOfComponents(1);
    marker->SetNumberOfTuples(grid->GetNumberOfCells());
    vtkNew<vtkDoubleArray> cellPressure;
    cellPressure->SetName("Pressure");
    cellPressure->SetNumberOfComponents(1);
    cellPressure->SetNumberOfTuples(grid->GetNumberOfCells());
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
        velocity->SetTuple3(cellId, velocityX, 0.0, 0.0);
        marker->SetValue(cellId, cellMarker);
        cellPressure->SetValue(cellId, pressure + 100.0);
    }
    grid->GetCellData()->AddArray(velocity);
    if (addStaleVelocityAlias) {
        vtkNew<vtkDoubleArray> staleVelocity;
        staleVelocity->SetName("Velocity");
        staleVelocity->SetNumberOfComponents(3);
        staleVelocity->SetNumberOfTuples(grid->GetNumberOfCells());
        for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); ++cellId) {
            staleVelocity->SetTuple3(cellId, -100.0, 0.0, 0.0);
        }
        grid->GetCellData()->AddArray(staleVelocity);
    }
    grid->GetCellData()->AddArray(marker);
    grid->GetCellData()->AddArray(cellPressure);
    return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> CreateLargeBoundsTinyVolumeTetra() {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(4);
    points->SetPoint(0, 0.0, 0.0, 0.0);
    points->SetPoint(1, 10.0, 0.0, 10.0);
    points->SetPoint(2, 0.0, 10.0, 10.0);
    // The fourth point is only slightly off the p0/p1/p2 plane. The cell has
    // a huge AABB but a physical volume of 100 * 0.1 / 6.
    points->SetPoint(3, 10.0, 10.0, 20.1);

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    vtkNew<vtkTetra> tetra;
    for (vtkIdType pointId = 0; pointId < 4; ++pointId) {
        tetra->GetPointIds()->SetId(pointId, pointId);
    }
    grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());

    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(4);
    velocity->FillComponent(0, 0.0);
    velocity->FillComponent(1, 0.0);
    velocity->FillComponent(2, 0.0);
    grid->GetPointData()->AddArray(velocity);
    return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> CreatePlanarQuad(double minimum,
                                                       double maximum) {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(4);
    points->SetPoint(0, minimum, minimum, 0.0);
    points->SetPoint(1, maximum, minimum, 0.0);
    points->SetPoint(2, maximum, maximum, 0.0);
    points->SetPoint(3, minimum, maximum, 0.0);

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    vtkNew<vtkQuad> quad;
    for (vtkIdType pointId = 0; pointId < 4; ++pointId) {
        quad->GetPointIds()->SetId(pointId, pointId);
    }
    grid->InsertNextCell(quad->GetCellType(), quad->GetPointIds());

    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(4);
    velocity->FillComponent(0, 0.0);
    velocity->FillComponent(1, 0.0);
    velocity->FillComponent(2, 0.0);
    grid->GetPointData()->AddArray(velocity);
    return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> CreateLargeBoundsTinyAreaTriangle() {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(3);
    points->SetPoint(0, 0.0, 0.0, 0.0);
    points->SetPoint(1, 10.0, 10.0, 0.0);
    // This triangle has a 10 by 10 AABB but only 0.005 physical area.
    points->SetPoint(2, 0.01, 0.011, 0.0);

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    vtkNew<vtkTriangle> triangle;
    for (vtkIdType pointId = 0; pointId < 3; ++pointId) {
        triangle->GetPointIds()->SetId(pointId, pointId);
    }
    grid->InsertNextCell(triangle->GetCellType(), triangle->GetPointIds());

    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(3);
    velocity->FillComponent(0, 0.0);
    velocity->FillComponent(1, 0.0);
    velocity->FillComponent(2, 0.0);
    grid->GetPointData()->AddArray(velocity);
    return grid;
}

void AddIntegerPointArray(vtkUnstructuredGrid* grid, const char* name, int value) {
    Require(grid != nullptr && name != nullptr,
            "Cannot add an integer point array to a null synthetic grid.");
    vtkNew<vtkIntArray> array;
    array->SetName(name);
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(grid->GetNumberOfPoints());
    array->FillValue(value);
    grid->GetPointData()->AddArray(array);
}

void AddFixedReferenceWeightArrays(vtkUnstructuredGrid* grid, double value) {
    Require(grid != nullptr, "Cannot add reference weights to a null synthetic grid.");
    for (const char* name : {io::kNek5000ReferenceWeightArray,
                             io::kNek5000JacobianDeterminantArray,
                             io::kNek5000QuadratureWeightArray,
                             "IntegrationWeight"}) {
        vtkNew<vtkDoubleArray> weight;
        weight->SetName(name);
        weight->SetNumberOfComponents(1);
        weight->SetNumberOfTuples(grid->GetNumberOfPoints());
        weight->FillValue(value);
        grid->GetPointData()->AddArray(weight);
    }
}

void WriteGrid(const fs::path& path, vtkUnstructuredGrid* grid) {
    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(path.string().c_str());
    writer->SetInputData(grid);
    Require(writer->Write() == 1, "Cannot write synthetic VTK grid: " + path.string());
}

void WriteNekCubeSource(const fs::path& descriptor,
                        const std::string& fieldPrefix,
                        double minimum,
                        double maximum) {
    std::ofstream descriptorOutput(descriptor);
    Require(static_cast<bool>(descriptorOutput),
            "Cannot write synthetic Nek5000 descriptor: " + descriptor.string());
    descriptorOutput << "filetemplate: " << fieldPrefix << "%01d.f%05d\n"
                     << "firsttimestep: 0\n"
                     << "numtimesteps: 1\n";
    descriptorOutput.close();

    const fs::path fieldPath = descriptor.parent_path() /
                               (fieldPrefix + "0.f00000");
    std::ofstream fieldOutput(fieldPath, std::ios::binary);
    Require(static_cast<bool>(fieldOutput),
            "Cannot write synthetic Nek5000 field: " + fieldPath.string());
    std::string header = "#std 8 2 2 2 1 1 0 0 0 1 X 0.0 F";
    header.resize(132, ' ');
    fieldOutput.write(header.data(), static_cast<std::streamsize>(header.size()));
    const float endianTag = 6.54321f;
    const std::int32_t elementId = 1;
    fieldOutput.write(reinterpret_cast<const char*>(&endianTag), sizeof(endianTag));
    fieldOutput.write(reinterpret_cast<const char*>(&elementId), sizeof(elementId));

    std::array<std::array<double, 8>, 3> coordinates{};
    std::size_t point = 0;
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i, ++point) {
                coordinates[0][point] = i == 0 ? minimum : maximum;
                coordinates[1][point] = j == 0 ? minimum : maximum;
                coordinates[2][point] = k == 0 ? minimum : maximum;
            }
        }
    }
    for (const std::array<double, 8>& component : coordinates) {
        fieldOutput.write(reinterpret_cast<const char*>(component.data()),
                          static_cast<std::streamsize>(component.size() * sizeof(double)));
    }
    Require(static_cast<bool>(fieldOutput),
            "Failed while writing synthetic Nek5000 field: " + fieldPath.string());
}

void RequireHdf(herr_t status, const std::string& operation) {
    Require(status >= 0, "Synthetic VTKHDF fixture failed to " + operation + ".");
}

void RequireHdfHandle(hid_t handle, const std::string& operation) {
    Require(handle >= 0, "Synthetic VTKHDF fixture failed to " + operation + ".");
}

void WriteVtkHdfTimeMetadata(const fs::path& path,
                             const std::array<double, 3>& times) {
    const hid_t file = H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC,
                                 H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(file, "create the HDF5 file");
    const hid_t vtkGroup = H5Gcreate2(file, "/VTKHDF", H5P_DEFAULT,
                                       H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(vtkGroup, "create the VTKHDF group");

    const hid_t scalarSpace = H5Screate(H5S_SCALAR);
    RequireHdfHandle(scalarSpace, "create scalar metadata space");
    const char typeName[] = "UnstructuredGrid";
    const hid_t typeNameType = H5Tcopy(H5T_C_S1);
    RequireHdfHandle(typeNameType, "create Type metadata type");
    RequireHdf(H5Tset_size(typeNameType, sizeof(typeName) - 1),
               "size Type metadata");
    RequireHdf(H5Tset_cset(typeNameType, H5T_CSET_ASCII),
               "mark Type metadata as ASCII");
    RequireHdf(H5Tset_strpad(typeNameType, H5T_STR_SPACEPAD),
               "configure Type metadata padding");
    const hid_t typeAttribute = H5Acreate2(vtkGroup, "Type", typeNameType,
                                           scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(typeAttribute, "create Type metadata");
    RequireHdf(H5Awrite(typeAttribute, typeNameType, typeName),
               "write Type metadata");
    RequireHdf(H5Aclose(typeAttribute), "close Type metadata");
    RequireHdf(H5Tclose(typeNameType), "close Type metadata type");

    const hsize_t versionDimensions[1] = {2};
    const hid_t versionSpace = H5Screate_simple(1, versionDimensions, nullptr);
    RequireHdfHandle(versionSpace, "create Version metadata space");
    const hid_t versionAttribute = H5Acreate2(vtkGroup, "Version", H5T_NATIVE_INT,
                                              versionSpace, H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(versionAttribute, "create Version metadata");
    const int version[2] = {2, 0};
    RequireHdf(H5Awrite(versionAttribute, H5T_NATIVE_INT, version),
               "write Version metadata");
    RequireHdf(H5Aclose(versionAttribute), "close Version metadata");
    RequireHdf(H5Sclose(versionSpace), "close Version metadata space");

    const hsize_t stepDimensions[1] = {times.size()};
    const hid_t pointsSpace = H5Screate_simple(1, stepDimensions, nullptr);
    RequireHdfHandle(pointsSpace, "create NumberOfPoints space");
    const hid_t points = H5Dcreate2(vtkGroup, "NumberOfPoints", H5T_NATIVE_LLONG,
                                    pointsSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(points, "create NumberOfPoints metadata");
    const long long pointCounts[3] = {0, 0, 0};
    RequireHdf(H5Dwrite(points, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, pointCounts),
               "write NumberOfPoints metadata");
    RequireHdf(H5Dclose(points), "close NumberOfPoints metadata");
    RequireHdf(H5Sclose(pointsSpace), "close NumberOfPoints space");

    const hid_t steps = H5Gcreate2(vtkGroup, "Steps", H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(steps, "create Steps group");
    const hid_t stepCountAttribute = H5Acreate2(steps, "NSteps", H5T_NATIVE_INT,
                                                scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(stepCountAttribute, "create NSteps metadata");
    const int stepCount = static_cast<int>(times.size());
    RequireHdf(H5Awrite(stepCountAttribute, H5T_NATIVE_INT, &stepCount),
               "write NSteps metadata");
    RequireHdf(H5Aclose(stepCountAttribute), "close NSteps metadata");
    const hid_t valuesSpace = H5Screate_simple(1, stepDimensions, nullptr);
    RequireHdfHandle(valuesSpace, "create time Values space");
    const hid_t values = H5Dcreate2(steps, "Values", H5T_NATIVE_DOUBLE,
                                    valuesSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    RequireHdfHandle(values, "create time Values metadata");
    RequireHdf(H5Dwrite(values, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, times.data()),
               "write time Values metadata");
    RequireHdf(H5Dclose(values), "close time Values metadata");
    RequireHdf(H5Sclose(valuesSpace), "close time Values space");
    RequireHdf(H5Gclose(steps), "close Steps group");
    RequireHdf(H5Sclose(scalarSpace), "close scalar metadata space");
    RequireHdf(H5Gclose(vtkGroup), "close VTKHDF group");
    RequireHdf(H5Fclose(file), "close the HDF5 file");
}

void WritePvd(const fs::path& path,
              const std::string& firstFile,
              double firstTime,
              const std::string& secondFile,
              double secondTime) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic PVD: " + path.string());
    output << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
           << "  <Collection>\n"
           << "    <DataSet timestep=\"" << firstTime << "\" file=\"" << firstFile
           << "\"/>\n"
           << "    <DataSet timestep=\"" << secondTime << "\" file=\"" << secondFile
           << "\"/>\n"
           << "  </Collection>\n"
           << "</VTKFile>\n";
}

void WriteManifest(const fs::path& path,
                   const std::string& overlapPolicy,
                   const std::string& foregroundSource,
                   const std::string& backgroundSource = "background.pvd",
                   const std::string& meshKind = "sliding") {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << std::setprecision(17);
    output << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=" << meshKind << "\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n"
           << "zone1_source=" << foregroundSource << "\n"
           << "zone1_motion=rigid\n"
           << "zone1_translation=2,2,2\n"
           << "zone1_rotation_axis=0,0,1\n"
           << "zone1_angular_velocity=" << (kPi / 2.0) << "\n"
           << "zone1_phase=0\n"
           << "zone1_velocity_frame=relative\n"
           << "zone1_velocity_components=uLocal,vLocal,wLocal\n"
           << "overlap_policy=" << overlapPolicy << "\n";
}

void WriteStaticWorldManifest(const fs::path& path,
                              const std::string& meshKind,
                              const std::string& backgroundSource = "background0.vtu",
                              const std::string& foregroundSource = "foreground-static.vtu",
                              int backgroundZone = -1) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=" << meshKind << "\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n"
           << "zone1_source=" << foregroundSource << "\n"
           << "zone1_motion=static\n"
           << "zone1_velocity_frame=world\n"
           << "overlap_policy=cut_background\n";
    if (backgroundZone >= 0) {
        output << "background_zone=" << backgroundZone << "\n";
    }
}

void WriteHugeAxisManifest(const fs::path& path,
                           const std::string& backgroundSource,
                           const std::string& foregroundSource) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << std::setprecision(17)
           << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=sliding\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n"
           << "zone1_source=" << foregroundSource << "\n"
           << "zone1_motion=rigid\n"
           << "zone1_rotation_axis=1e308,1e308,1e308\n"
           << "zone1_angular_velocity=0\n"
           << "zone1_phase=" << kPi << "\n"
           << "zone1_velocity_frame=world\n"
           << "overlap_policy=cut_background\n";
}

void WriteStaticVelocityManifest(const fs::path& path,
                                 const std::string& backgroundSource,
                                 const std::string& foregroundSource,
                                 const std::string& overlapPolicy,
                                 const std::string& velocityArray = {},
                                 const std::string& velocityComponents = {}) {
    Require(velocityArray.empty() || velocityComponents.empty(),
            "Synthetic manifest cannot declare both velocity-array forms.");
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=nested\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n";
    if (!velocityArray.empty()) {
        output << "zone0_velocity_array=" << velocityArray << "\n";
    }
    if (!velocityComponents.empty()) {
        output << "zone0_velocity_components=" << velocityComponents << "\n";
    }
    output << "zone1_source=" << foregroundSource << "\n"
           << "zone1_motion=static\n"
           << "zone1_velocity_frame=world\n";
    if (!velocityArray.empty()) {
        output << "zone1_velocity_array=" << velocityArray << "\n";
    }
    if (!velocityComponents.empty()) {
        output << "zone1_velocity_components=" << velocityComponents << "\n";
    }
    output << "overlap_policy=" << overlapPolicy << "\n";
}

void WriteScalarOnlyManifest(const fs::path& path,
                             const std::string& backgroundSource,
                             bool declareMissingVelocity = false) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=nested\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n";
    if (declareMissingVelocity) {
        output << "zone0_velocity_array=MissingVelocity\n";
    }
    output << "zone1_source=foreground-scalar.vtu\n"
           << "zone1_motion=static\n"
           << "zone1_translation=2,2,2\n"
           << "zone1_velocity_frame=world\n"
           << "overlap_policy=cut_background\n";
}

void WriteCellVelocityManifest(const fs::path& path,
                               const std::string& backgroundSource =
                                   "background-cell-velocity.vtu",
                               const std::string& foregroundSource =
                                   "foreground-cell-velocity.vtu",
                               bool useCustomVelocity = false) {
    std::ofstream output(path);
    Require(static_cast<bool>(output), "Cannot write synthetic .scpm: " + path.string());
    output << "streamcenterplus_meshset=1\n"
           << "zone_count=2\n"
           << "mesh_kind=nested\n"
           << "zone0_source=" << backgroundSource << "\n"
           << "zone0_motion=static\n"
           << "zone0_velocity_frame=world\n";
    if (useCustomVelocity) {
        output << "zone0_velocity_array=MyVelocity\n";
    }
    output << "zone1_source=" << foregroundSource << "\n"
           << "zone1_motion=static\n"
           << "zone1_velocity_frame=world\n";
    if (useCustomVelocity) {
        output << "zone1_velocity_array=MyVelocity\n";
    }
    output << "overlap_policy=cut_background\n";
}

vtkIdType FindPoint(vtkDataSet* data, double x, double y, double z) {
    double point[3] = {};
    for (vtkIdType id = 0; id < data->GetNumberOfPoints(); ++id) {
        data->GetPoint(id, point);
        if (NearlyEqual(point[0], x) && NearlyEqual(point[1], y) &&
            NearlyEqual(point[2], z)) {
            return id;
        }
    }
    throw std::runtime_error("Synthetic reference point was not found.");
}

void RequireCoveredVelocity(vtkDataSet* data, bool blend) {
    vtkDataArray* velocity = data->GetPointData()->GetArray("Velocity");
    vtkDataArray* pressure = data->GetPointData()->GetArray("Pressure");
    vtkDataArray* ownership = data->GetPointData()->GetArray(io::kMultiZoneOwnershipArray);
    vtkDataArray* weight = data->GetPointData()->GetArray(io::kMultiZoneBlendWeightArray);
    Require(velocity != nullptr && velocity->GetNumberOfComponents() == 3 &&
                pressure != nullptr && ownership != nullptr && weight != nullptr,
            "Canonical two-zone reader did not expose world Velocity and provenance.");
    Require(data->GetPointData()->GetArray("uLocal") == nullptr &&
                data->GetPointData()->GetArray("vLocal") == nullptr &&
                data->GetPointData()->GetArray("wLocal") == nullptr,
            "Canonical reader retained a zone-frame velocity candidate.");
    const vtkIdType pointId = FindPoint(data, 3.0, 2.0, 2.0);
    const double foregroundY = 1.0 + kPi / 2.0;
    const double blendFactor = blend ? 0.5 : 1.0;
    Require(NearlyEqual(velocity->GetComponent(pointId, 0), 10.0 * (1.0 - blendFactor)) &&
                NearlyEqual(velocity->GetComponent(pointId, 1), foregroundY * blendFactor) &&
                NearlyEqual(pressure->GetComponent(pointId, 0),
                            (1.0 - blendFactor) * 5.0 + blendFactor * 20.0) &&
                NearlyEqual(ownership->GetComponent(pointId, 0), blend ? 2.0 : 1.0) &&
                NearlyEqual(weight->GetComponent(pointId, 0), blendFactor),
            "Sliding relative velocity was not transformed/merged in the world reference frame.");
    Require(data->GetCellData()->GetArray(io::kMultiZoneBackgroundCutArray) != nullptr &&
                data->GetPointData()->GetArray(io::kMultiZoneInterfaceArray) != nullptr &&
                data->GetCellData()->GetArray(io::kMultiZoneInterfaceArray) != nullptr,
            "Canonical two-zone reader is missing interface/cut diagnostics.");
}

void ValidateStaticMeshKinds(const fs::path& directory,
                             vtkUnstructuredGrid* background) {
    vtkSmartPointer<vtkUnstructuredGrid> foreground =
        CreateBoxGrid(1.0, 3.0, 2, false, 20.0, 3.0);
    WriteGrid(directory / "foreground-static.vtu", foreground);

    for (const std::string& meshKind : {"static", "nested", "multigrid"}) {
        const fs::path manifest = directory / (meshKind + ".scpm");
        WriteStaticWorldManifest(manifest, meshKind);
        const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(manifest);
        Require(frames.size() == 1 &&
                    frames.front().kind == io::DataSetFrameKind::MultiZoneManifestTimeStep &&
                    frames.front().multiZoneBackgroundIndex == 0,
                "Static two-zone " + meshKind +
                    " manifest did not enumerate a canonical frame.");
        vtkSmartPointer<vtkDataSet> data = io::ReadDataSet(frames.front());
        vtkDataArray* velocity = data->GetPointData()->GetArray(
            io::kMultiZoneCanonicalVelocityArray);
        Require(data->GetNumberOfPoints() == background->GetNumberOfPoints() &&
                    velocity != nullptr && velocity->GetNumberOfComponents() == 3 &&
                    NearlyEqual(velocity->GetComponent(FindPoint(data, 2.0, 2.0, 2.0), 0),
                                3.0),
                "Static world-velocity " + meshKind +
                    " manifest did not read as a two-zone mesh set.");
    }
}

void ValidateScalarOnlyManifest(const fs::path& directory,
                                vtkUnstructuredGrid* background) {
    vtkSmartPointer<vtkUnstructuredGrid> scalarBackground =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> scalarForeground =
        CreateBoxGrid(-1.5, 1.5, 3, false, 20.0, 1.0);
    scalarBackground->GetPointData()->RemoveArray("Velocity");
    scalarForeground->GetPointData()->RemoveArray("Velocity");
    WriteGrid(directory / "background-scalar.vtu", scalarBackground);
    WriteGrid(directory / "foreground-scalar.vtu", scalarForeground);

    const fs::path scalarManifest = directory / "scalar-only.scpm";
    WriteScalarOnlyManifest(scalarManifest, "background-scalar.vtu");
    const std::vector<io::DataSetFrame> scalarFrames =
        io::EnumerateDataSetFrames(scalarManifest);
    Require(scalarFrames.size() == 1 &&
                scalarFrames.front().kind == io::DataSetFrameKind::MultiZoneManifestTimeStep,
            "Scalar-only two-zone manifest did not enumerate.");
    vtkSmartPointer<vtkDataSet> scalarData = io::ReadDataSet(scalarFrames.front());
    vtkDataArray* pressure = scalarData->GetPointData()->GetArray("Pressure");
    Require(scalarData->GetNumberOfPoints() == background->GetNumberOfPoints() &&
                pressure != nullptr &&
                NearlyEqual(pressure->GetComponent(FindPoint(scalarData, 2.0, 2.0, 2.0), 0),
                            20.0) &&
                scalarData->GetPointData()->GetArray(io::kMultiZoneCanonicalVelocityArray) ==
                    nullptr,
            "Scalar-only two-zone input did not transform geometry and compose scalar data "
            "without inventing Velocity.");

    const fs::path partialVelocityManifest = directory / "partial-velocity.scpm";
    WriteScalarOnlyManifest(partialVelocityManifest, "background0.vtu");
    vtkSmartPointer<vtkDataSet> partialVelocity = io::ReadDataSet(
        io::EnumerateDataSetFrames(partialVelocityManifest).front());
    Require(partialVelocity->GetPointData()->GetArray(
                io::kMultiZoneCanonicalVelocityArray) == nullptr,
            "A one-zone velocity was exposed as a complete canonical Velocity.");

    const fs::path invalidVelocityManifest = directory / "invalid-velocity.scpm";
    WriteScalarOnlyManifest(invalidVelocityManifest, "background-scalar.vtu", true);
    bool rejected = false;
    try {
        static_cast<void>(io::ReadDataSet(
            io::EnumerateDataSetFrames(invalidVelocityManifest).front()));
    } catch (const std::exception&) {
        rejected = true;
    }
    Require(rejected, "An invalid explicitly declared velocity was accepted.");
}

void ValidateCellVelocityResampler(const fs::path& directory) {
    vtkSmartPointer<vtkUnstructuredGrid> background =
        CreateCellVelocityBoxGrid(0.0, 4.0, 4, 5.0, 10.0, 7.0);
    vtkSmartPointer<vtkUnstructuredGrid> foreground =
        CreateCellVelocityBoxGrid(1.0, 3.0, 2, 20.0, 3.0, 9.0);
    AddFixedReferenceWeightArrays(background, 7.0);
    AddFixedReferenceWeightArrays(foreground, 11.0);
    WriteGrid(directory / "background-cell-velocity.vtu", background);
    WriteGrid(directory / "foreground-cell-velocity.vtu", foreground);
    const fs::path manifest = directory / "cell-velocity.scpm";
    WriteCellVelocityManifest(manifest);

    vtkSmartPointer<vtkDataSet> composed = io::ReadDataSet(
        io::EnumerateDataSetFrames(manifest).front());
    Require(composed->GetCellData()->GetArray("Velocity") == nullptr &&
                composed->GetCellData()->GetArray("CellMarker") != nullptr &&
                composed->GetCellData()->GetArray("Pressure") != nullptr,
            "Mesh-set cell velocity cleanup removed the wrong cell arrays.");

    const vtkIdType coveredPoint = FindPoint(composed, 2.0, 2.0, 2.0);
    for (const char* name : {io::kNek5000ReferenceWeightArray,
                             io::kNek5000JacobianDeterminantArray,
                             io::kNek5000QuadratureWeightArray,
                             "IntegrationWeight"}) {
        vtkDataArray* weight = composed->GetPointData()->GetArray(name);
        Require(weight != nullptr && NearlyEqual(weight->GetComponent(coveredPoint, 0), 7.0),
                "A fixed-background physical weight was composed from the foreground.");
    }

    vtkSmartPointer<vtkDataSet> prepared = resample::PreparePointData(composed);
    vtkDataArray* velocity = prepared->GetPointData()->GetArray(
        io::kMultiZoneCanonicalVelocityArray);
    Require(velocity != nullptr && velocity->GetNumberOfComponents() == 3 &&
                NearlyEqual(velocity->GetComponent(FindPoint(prepared, 2.0, 2.0, 2.0), 0),
                            3.0),
            "PreparePointData replaced canonical point Velocity with stale cell velocity.");
    vtkDataArray* pressure = prepared->GetPointData()->GetArray("Pressure");
    Require(pressure != nullptr &&
                NearlyEqual(pressure->GetComponent(FindPoint(prepared, 2.0, 2.0, 2.0), 0),
                            20.0),
            "PreparePointData replaced a composed point scalar with CellData of the same name.");

    vtkSmartPointer<vtkUnstructuredGrid> customBackground = CreateCellVelocityBoxGrid(
        0.0, 4.0, 4, 5.0, 10.0, 7.0, "MyVelocity", true);
    vtkSmartPointer<vtkUnstructuredGrid> customForeground = CreateCellVelocityBoxGrid(
        1.0, 3.0, 2, 20.0, 3.0, 9.0, "MyVelocity", true);
    WriteGrid(directory / "background-cell-custom-velocity.vtu", customBackground);
    WriteGrid(directory / "foreground-cell-custom-velocity.vtu", customForeground);
    const fs::path customManifest = directory / "cell-custom-velocity.scpm";
    WriteCellVelocityManifest(customManifest, "background-cell-custom-velocity.vtu",
                              "foreground-cell-custom-velocity.vtu", true);
    vtkSmartPointer<vtkDataSet> customComposed = io::ReadDataSet(
        io::EnumerateDataSetFrames(customManifest).front());
    Require(customComposed->GetCellData()->GetArray("MyVelocity") == nullptr &&
                customComposed->GetCellData()->GetArray("Velocity") == nullptr &&
                customComposed->GetCellData()->GetArray("CellMarker") != nullptr,
            "Custom cell velocity cleanup retained a stale standard alias.");
    vtkSmartPointer<vtkDataSet> customPrepared =
        resample::PreparePointData(customComposed);
    vtkDataArray* customVelocity = customPrepared->GetPointData()->GetArray(
        io::kMultiZoneCanonicalVelocityArray);
    Require(customVelocity != nullptr && customVelocity->GetNumberOfComponents() == 3 &&
                NearlyEqual(customVelocity->GetComponent(
                                FindPoint(customPrepared, 2.0, 2.0, 2.0), 0),
                            3.0),
            "A stale CellData Velocity alias replaced custom canonical velocity.");
}

vtkSmartPointer<vtkUnstructuredGrid> AppendDisjointBackground(
    vtkUnstructuredGrid* first,
    vtkUnstructuredGrid* second) {
    Require(first != nullptr && second != nullptr,
            "Cannot append null synthetic background grids.");
    vtkNew<vtkAppendFilter> append;
    append->MergePointsOff();
    append->AddInputData(first);
    append->AddInputData(second);
    append->Update();
    vtkSmartPointer<vtkUnstructuredGrid> output =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    output->DeepCopy(append->GetOutput());
    return output;
}

vtkSmartPointer<vtkUnstructuredGrid> CreateCrossingGapTetra() {
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->InsertNextPoint(0.75, 0.25, 0.25);
    points->InsertNextPoint(0.75, 0.75, 0.25);
    points->InsertNextPoint(0.75, 0.25, 0.75);
    points->InsertNextPoint(1.20, 0.30, 0.30);

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    vtkNew<vtkTetra> tetra;
    for (vtkIdType pointId = 0; pointId < 4; ++pointId) {
        tetra->GetPointIds()->SetId(pointId, pointId);
    }
    grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());

    vtkNew<vtkDoubleArray> pressure;
    pressure->SetName("Pressure");
    pressure->SetNumberOfComponents(1);
    pressure->SetNumberOfTuples(4);
    pressure->FillValue(20.0);
    grid->GetPointData()->AddArray(pressure);
    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(4);
    velocity->FillValue(3.0);
    grid->GetPointData()->AddArray(velocity);
    return grid;
}

void RequireManifestReadRejected(const fs::path& manifest,
                                 const std::string& expectedError,
                                 const std::string& message) {
    bool rejected = false;
    std::string errorText;
    try {
        const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(manifest);
        Require(!frames.empty(), "Synthetic mesh-set manifest unexpectedly had no frames.");
        static_cast<void>(io::ReadDataSet(frames.back()));
    } catch (const std::exception& error) {
        rejected = true;
        errorText = error.what();
    }
    Require(rejected && errorText.find(expectedError) != std::string::npos, message);
}

void ValidatePhysicalVolumeBackgroundSelection(const fs::path& directory) {
    vtkSmartPointer<vtkUnstructuredGrid> largeBoundsTinyVolume =
        CreateLargeBoundsTinyVolumeTetra();
    vtkSmartPointer<vtkUnstructuredGrid> compactPhysicalBackground =
        CreateBoxGrid(0.0, 2.0, 1, false, 5.0, 10.0);
    WriteGrid(directory / "large-bounds-tiny-volume.vtu", largeBoundsTinyVolume);
    WriteGrid(directory / "compact-physical-background.vtu", compactPhysicalBackground);

    const fs::path manifest = directory / "physical-volume-selection.scpm";
    WriteStaticWorldManifest(manifest, "nested", "large-bounds-tiny-volume.vtu",
                             "compact-physical-background.vtu");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(manifest);
    Require(frames.size() == 1 && frames.front().multiZoneBackgroundIndex == 1,
            "Default mesh-set background selection used an AABB instead of physical 3-D volume.");

    vtkSmartPointer<vtkUnstructuredGrid> largeBoundsTinyArea =
        CreateLargeBoundsTinyAreaTriangle();
    vtkSmartPointer<vtkUnstructuredGrid> compactPlanarBackground =
        CreatePlanarQuad(0.0, 2.0);
    WriteGrid(directory / "large-bounds-tiny-area.vtu", largeBoundsTinyArea);
    WriteGrid(directory / "compact-planar-background.vtu", compactPlanarBackground);

    const fs::path planarManifest = directory / "physical-area-selection.scpm";
    WriteStaticWorldManifest(planarManifest, "nested", "large-bounds-tiny-area.vtu",
                             "compact-planar-background.vtu");
    const std::vector<io::DataSetFrame> planarFrames =
        io::EnumerateDataSetFrames(planarManifest);
    Require(planarFrames.size() == 1 && planarFrames.front().multiZoneBackgroundIndex == 1,
            "Automatic mesh-set background selection rejected or miscompared valid 2-D "
            "reference areas.");
}

void ValidateNekAutomaticBackgroundSelection(const fs::path& directory) {
    const fs::path nekDescriptor = directory / "auto-background.nek5000";
    WriteNekCubeSource(nekDescriptor, "auto-background", 0.0, 4.0);
    vtkSmartPointer<vtkUnstructuredGrid> compactVtkZone =
        CreateBoxGrid(0.0, 2.0, 1, false, 5.0, 10.0);
    WriteGrid(directory / "nek-auto-compact-zone.vtu", compactVtkZone);

    const fs::path manifest = directory / "nek-auto-background.scpm";
    WriteStaticWorldManifest(manifest, "nested", nekDescriptor.filename().string(),
                             "nek-auto-compact-zone.vtu");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(manifest);
    Require(frames.size() == 1 && frames.front().multiZoneBackgroundIndex == 0,
            "Automatic mesh-set background selection did not build Nek5000 reference "
            "topology before measuring its physical volume.");
}

void ValidateVtkHdfNativeTimes(const fs::path& directory) {
    const fs::path hdfDirectory = directory / "vtkhdf-native-times";
    fs::create_directories(hdfDirectory);
    const fs::path source = hdfDirectory / "nonuniform.vtkhdf";
    const std::array<double, 3> expectedTimes = {0.0, 0.125, 0.625};
    WriteVtkHdfTimeMetadata(source, expectedTimes);

    const auto requireNativeTimes = [&](const std::vector<io::DataSetFrame>& frames,
                                        const char* description) {
        Require(frames.size() == expectedTimes.size(),
                std::string(description) + " did not enumerate every VTKHDF step.");
        for (std::size_t step = 0; step < expectedTimes.size(); ++step) {
            Require(NearlyEqual(frames[step].time, expectedTimes[step]) &&
                        frames[step].readerTime.has_value() &&
                        NearlyEqual(*frames[step].readerTime, expectedTimes[step]),
                    std::string(description) +
                        " replaced non-uniform VTKHDF time metadata with fallback cadence.");
        }
    };
    requireNativeTimes(io::EnumerateDataSetFrames(source, 7.0),
                       "Direct VTKHDF enumeration");
    requireNativeTimes(io::EnumerateDataSetFrames(hdfDirectory, 7.0),
                       "Folder VTKHDF enumeration");
}

void ValidateVelocityTupleGuards(const fs::path& directory) {
    const auto truncate = [](vtkDataArray* array) {
        Require(array != nullptr && array->GetNumberOfTuples() > 1,
                "Synthetic velocity array cannot be truncated.");
        array->SetNumberOfTuples(array->GetNumberOfTuples() - 1);
    };

    vtkSmartPointer<vtkUnstructuredGrid> autoBackground =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> autoForeground =
        CreateBoxGrid(1.0, 3.0, 2, false, 20.0, 3.0);
    truncate(autoBackground->GetPointData()->GetArray("Velocity"));
    WriteGrid(directory / "truncated-auto-background.vtu", autoBackground);
    WriteGrid(directory / "truncated-auto-foreground.vtu", autoForeground);
    const fs::path autoManifest = directory / "truncated-auto-velocity.scpm";
    WriteStaticWorldManifest(autoManifest, "nested", "truncated-auto-background.vtu",
                             "truncated-auto-foreground.vtu");
    RequireManifestReadRejected(
        autoManifest, "Cannot read VTK XML data",
        "A truncated automatically selected point Velocity array was accepted.");

    vtkSmartPointer<vtkUnstructuredGrid> explicitBackground =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> explicitForeground =
        CreateBoxGrid(1.0, 3.0, 2, false, 20.0, 3.0);
    explicitBackground->GetPointData()->GetArray("Velocity")->SetName("MyVelocity");
    explicitForeground->GetPointData()->GetArray("Velocity")->SetName("MyVelocity");
    truncate(explicitBackground->GetPointData()->GetArray("MyVelocity"));
    WriteGrid(directory / "truncated-explicit-background.vtu", explicitBackground);
    WriteGrid(directory / "truncated-explicit-foreground.vtu", explicitForeground);
    const fs::path explicitManifest = directory / "truncated-explicit-velocity.scpm";
    WriteStaticVelocityManifest(explicitManifest, "truncated-explicit-background.vtu",
                                "truncated-explicit-foreground.vtu", "cut_background",
                                "MyVelocity");
    RequireManifestReadRejected(
        explicitManifest, "Cannot read VTK XML data",
        "A truncated explicitly selected point Velocity array was accepted.");

    vtkSmartPointer<vtkUnstructuredGrid> componentBackground =
        CreateBoxGrid(0.0, 4.0, 4, true, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> componentForeground =
        CreateBoxGrid(1.0, 3.0, 2, true, 20.0, 3.0);
    truncate(componentBackground->GetPointData()->GetArray("uLocal"));
    WriteGrid(directory / "truncated-components-background.vtu", componentBackground);
    WriteGrid(directory / "truncated-components-foreground.vtu", componentForeground);
    const fs::path componentManifest = directory / "truncated-components-velocity.scpm";
    WriteStaticVelocityManifest(componentManifest, "truncated-components-background.vtu",
                                "truncated-components-foreground.vtu", "cut_background", {},
                                "uLocal,vLocal,wLocal");
    RequireManifestReadRejected(
        componentManifest, "Cannot read VTK XML data",
        "A truncated explicitly selected scalar velocity component was accepted.");

    vtkSmartPointer<vtkUnstructuredGrid> cellBackground =
        CreateCellVelocityBoxGrid(0.0, 4.0, 4, 5.0, 10.0, 7.0);
    vtkSmartPointer<vtkUnstructuredGrid> cellForeground =
        CreateCellVelocityBoxGrid(1.0, 3.0, 2, 20.0, 3.0, 9.0);
    truncate(cellBackground->GetCellData()->GetArray("Velocity"));
    WriteGrid(directory / "truncated-cell-background.vtu", cellBackground);
    WriteGrid(directory / "truncated-cell-foreground.vtu", cellForeground);
    const fs::path cellManifest = directory / "truncated-cell-velocity.scpm";
    WriteCellVelocityManifest(cellManifest, "truncated-cell-background.vtu",
                              "truncated-cell-foreground.vtu");
    RequireManifestReadRejected(
        cellManifest, "Cannot read VTK XML data",
        "A truncated cell Velocity array was accepted.");
}

void ValidateBlendFieldGuards(const fs::path& directory) {
    vtkSmartPointer<vtkUnstructuredGrid> background =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> foreground =
        CreateBoxGrid(1.0, 3.0, 2, false, 20.0, 3.0);
    AddIntegerPointArray(background, "MaterialId", 1);
    AddIntegerPointArray(foreground, "MaterialId", 2);
    WriteGrid(directory / "blend-integer-background.vtu", background);
    WriteGrid(directory / "blend-integer-foreground.vtu", foreground);
    const fs::path manifest = directory / "blend-integer.scpm";
    WriteStaticVelocityManifest(manifest, "blend-integer-background.vtu",
                                "blend-integer-foreground.vtu", "blend_50_50");
    vtkSmartPointer<vtkDataSet> composed = io::ReadDataSet(
        io::EnumerateDataSetFrames(manifest).front());
    Require(composed->GetPointData()->GetArray("MaterialId") == nullptr,
            "A non-floating PointData field survived a 50:50 mesh-set blend.");
    Require(composed->GetPointData()->GetArray("Pressure") != nullptr,
            "A floating-point mesh-set field was incorrectly removed from a 50:50 blend.");
}

void ValidateFailClosedMeshSetGuards(const fs::path& directory) {
    vtkSmartPointer<vtkUnstructuredGrid> coarseBackground =
        CreateBoxGrid(0.0, 4.0, 1, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> embeddedForeground =
        CreateBoxGrid(1.75, 2.25, 1, false, 20.0, 3.0);
    WriteGrid(directory / "coarse-background.vtu", coarseBackground);
    WriteGrid(directory / "embedded-foreground.vtu", embeddedForeground);
    const fs::path unrepresentableManifest = directory / "unrepresentable.scpm";
    WriteStaticWorldManifest(unrepresentableManifest, "nested", "coarse-background.vtu",
                             "embedded-foreground.vtu");
    RequireManifestReadRejected(
        unrepresentableManifest,
        "contains no fixed background degree of freedom",
        "A foreground cell without a fixed background degree of freedom was accepted.");

    vtkSmartPointer<vtkUnstructuredGrid> disjointBackground = AppendDisjointBackground(
        CreateBoxGrid(0.0, 1.0, 1, false, 5.0, 10.0),
        CreateBoxGrid(3.0, 4.0, 1, false, 5.0, 10.0));
    WriteGrid(directory / "disjoint-background.vtu", disjointBackground);
    const fs::path physicalContainmentManifest = directory / "physical-containment.scpm";
    WriteStaticWorldManifest(physicalContainmentManifest, "nested", "disjoint-background.vtu",
                             "embedded-foreground.vtu");
    RequireManifestReadRejected(
        physicalContainmentManifest,
        "physical background volume",
        "A foreground in a background AABB hole was accepted as physically contained.");

    vtkSmartPointer<vtkUnstructuredGrid> edgeCrossingRight =
        CreateBoxGrid(1.1, 2.1, 8, false, 5.0, 10.0);
    for (vtkIdType pointId = 0;
         pointId < edgeCrossingRight->GetNumberOfPoints(); ++pointId) {
        double point[3] = {};
        edgeCrossingRight->GetPoint(pointId, point);
        edgeCrossingRight->GetPoints()->SetPoint(pointId, point[0], point[1] - 1.1,
                                                  point[2] - 1.1);
    }
    edgeCrossingRight->GetPoints()->Modified();
    edgeCrossingRight->Modified();
    vtkSmartPointer<vtkUnstructuredGrid> edgeCrossingBackground =
        AppendDisjointBackground(
            CreateBoxGrid(0.0, 1.0, 8, false, 5.0, 10.0), edgeCrossingRight);
    vtkSmartPointer<vtkUnstructuredGrid> edgeCrossingForeground =
        CreateCrossingGapTetra();
    WriteGrid(directory / "edge-crossing-background.vtu", edgeCrossingBackground);
    WriteGrid(directory / "edge-crossing-foreground.vtu", edgeCrossingForeground);
    vtkSmartPointer<vtkDataSet> persistedEdgeCrossingBackground = io::ReadDataSet(
        directory / "edge-crossing-background.vtu");
    vtkNew<vtkStaticCellLocator> edgeCrossingLocator;
    edgeCrossingLocator->SetDataSet(persistedEdgeCrossingBackground);
    edgeCrossingLocator->BuildLocator();
    vtkNew<vtkGenericCell> edgeCrossingCell;
    std::vector<double> edgeCrossingWeights(static_cast<std::size_t>(
        std::max<vtkIdType>(1, persistedEdgeCrossingBackground->GetMaxCellSize())));
    double edgeCrossingGapPoint[3] = {1.05, 0.35, 0.35};
    double edgeCrossingPcoords[3] = {};
    int edgeCrossingSubId = 0;
    Require(edgeCrossingLocator->FindCell(edgeCrossingGapPoint, 0.0, edgeCrossingCell,
                                          edgeCrossingSubId, edgeCrossingPcoords,
                                          edgeCrossingWeights.data()) < 0,
            "Synthetic edge-crossing fixture did not retain the intended background gap.");
    double edgeCrossingCentre[3] = {};
    for (vtkIdType pointId = 0;
         pointId < edgeCrossingForeground->GetNumberOfPoints(); ++pointId) {
        double point[3] = {};
        edgeCrossingForeground->GetPoint(pointId, point);
        edgeCrossingCentre[0] += point[0];
        edgeCrossingCentre[1] += point[1];
        edgeCrossingCentre[2] += point[2];
        edgeCrossingSubId = 0;
        edgeCrossingPcoords[0] = edgeCrossingPcoords[1] = edgeCrossingPcoords[2] = 0.0;
        Require(edgeCrossingLocator->FindCell(point, 0.0, edgeCrossingCell,
                                              edgeCrossingSubId, edgeCrossingPcoords,
                                              edgeCrossingWeights.data()) >= 0,
                "Synthetic edge-crossing foreground vertex was not covered by its background.");
    }
    const double pointCount = static_cast<double>(edgeCrossingForeground->GetNumberOfPoints());
    edgeCrossingCentre[0] /= pointCount;
    edgeCrossingCentre[1] /= pointCount;
    edgeCrossingCentre[2] /= pointCount;
    edgeCrossingSubId = 0;
    edgeCrossingPcoords[0] = edgeCrossingPcoords[1] = edgeCrossingPcoords[2] = 0.0;
    Require(edgeCrossingLocator->FindCell(edgeCrossingCentre, 0.0, edgeCrossingCell,
                                          edgeCrossingSubId, edgeCrossingPcoords,
                                          edgeCrossingWeights.data()) >= 0,
            "Synthetic edge-crossing foreground centre was not covered by its background.");
    const fs::path edgeCrossingManifest = directory / "edge-crossing-gap.scpm";
    WriteStaticWorldManifest(edgeCrossingManifest, "nested", "edge-crossing-background.vtu",
                             "edge-crossing-foreground.vtu", 0);
    const std::vector<io::DataSetFrame> edgeCrossingFrames =
        io::EnumerateDataSetFrames(edgeCrossingManifest);
    Require(edgeCrossingFrames.size() == 1 &&
                edgeCrossingFrames.front().multiZoneBackgroundIndex == 0,
            "Synthetic edge-crossing manifest did not retain zone 0 as its background.");
    RequireManifestReadRejected(
        edgeCrossingManifest,
        "cell edge crosses outside the physical background volume",
        "A foreground cell crossing a background gap was accepted when its vertices and centre were covered.");

    vtkSmartPointer<vtkUnstructuredGrid> hugeAxisBackground =
        CreateBoxGrid(-2.0, 2.0, 8, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> hugeAxisForeground =
        CreateBoxGrid(0.75, 1.25, 1, false, 20.0, 3.0);
    WriteGrid(directory / "huge-axis-background.vtu", hugeAxisBackground);
    WriteGrid(directory / "huge-axis-foreground.vtu", hugeAxisForeground);
    const fs::path hugeAxisManifest = directory / "huge-axis.scpm";
    WriteHugeAxisManifest(hugeAxisManifest, "huge-axis-background.vtu",
                          "huge-axis-foreground.vtu");
    const std::vector<io::DataSetFrame> hugeAxisFrames =
        io::EnumerateDataSetFrames(hugeAxisManifest);
    Require(hugeAxisFrames.size() == 1,
            "Huge finite rotation-axis manifest did not enumerate one frame.");
    vtkSmartPointer<vtkDataSet> hugeAxisComposed =
        io::ReadDataSet(hugeAxisFrames.front());
    vtkDataArray* hugeAxisOwnership = hugeAxisComposed->GetPointData()->GetArray(
        io::kMultiZoneOwnershipArray);
    const vtkIdType hugeAxisCentre = FindPoint(hugeAxisComposed, 1.0, 1.0, 1.0);
    Require(hugeAxisOwnership != nullptr &&
                NearlyEqual(hugeAxisOwnership->GetComponent(hugeAxisCentre, 0), 1.0),
            "A huge finite rotation axis did not preserve its direction during world transform.");

    // A PVD indirection must not turn an already-composed two-zone manifest
    // into a source zone of another manifest.  That would silently exceed the
    // two-zone contract and transform canonical world-frame fields twice.
    const fs::path indirectInnerManifest = directory / "indirect-inner.scpm";
    WriteStaticWorldManifest(indirectInnerManifest, "nested", "coarse-background.vtu",
                             "embedded-foreground.vtu");
    WritePvd(directory / "indirect-inner.pvd", "indirect-inner.scpm", 0.0,
             "indirect-inner.scpm", 1.0);
    const fs::path indirectOuterManifest = directory / "indirect-outer.scpm";
    WriteStaticWorldManifest(indirectOuterManifest, "nested", "indirect-inner.pvd",
                             "embedded-foreground.vtu");
    RequireManifestReadRejected(
        indirectOuterManifest,
        "Nested .scpm sources",
        "A PVD-indirected nested .scpm source was accepted as a two-zone manifest.");

    // A physical field must be available in both zones.  Leaving a
    // background-only field attached would silently use its stale values in
    // the foreground-covered region.
    vtkSmartPointer<vtkUnstructuredGrid> schemaBackground =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> schemaForeground =
        CreateBoxGrid(1.0, 3.0, 2, false, 20.0, 3.0);
    schemaForeground->GetPointData()->RemoveArray("Pressure");
    vtkNew<vtkDoubleArray> staleBackgroundCellPressure;
    staleBackgroundCellPressure->SetName("Pressure");
    staleBackgroundCellPressure->SetNumberOfComponents(1);
    staleBackgroundCellPressure->SetNumberOfTuples(schemaBackground->GetNumberOfCells());
    staleBackgroundCellPressure->FillValue(105.0);
    schemaBackground->GetCellData()->AddArray(staleBackgroundCellPressure);
    WriteGrid(directory / "schema-background.vtu", schemaBackground);
    WriteGrid(directory / "schema-foreground.vtu", schemaForeground);
    const fs::path missingFieldManifest = directory / "missing-field.scpm";
    WriteStaticWorldManifest(missingFieldManifest, "nested", "schema-background.vtu",
                             "schema-foreground.vtu");
    const std::vector<io::DataSetFrame> missingFieldFrames =
        io::EnumerateDataSetFrames(missingFieldManifest);
    vtkSmartPointer<vtkDataSet> missingFieldData =
        io::ReadDataSet(missingFieldFrames.front());
    Require(missingFieldData->GetPointData()->GetArray("Pressure") == nullptr,
            "A background-only PointData field survived multi-zone composition.");
    Require(missingFieldData->GetPointData()->GetArray(
                io::kMultiZoneCanonicalVelocityArray) != nullptr,
            "Removing an unavailable PointData field removed a shared canonical field.");
    Require(missingFieldData->GetCellData()->GetArray("Pressure") != nullptr,
            "The stale CellData regression fixture was not retained by the reader.");
    vtkSmartPointer<vtkDataSet> missingFieldPrepared =
        resample::PreparePointData(missingFieldData);
    Require(missingFieldPrepared->GetPointData()->GetArray("Pressure") == nullptr &&
                missingFieldPrepared->GetCellData()->GetNumberOfArrays() == 0,
            "Resample promoted a background-only multi-zone CellData field to PointData.");

    vtkSmartPointer<vtkUnstructuredGrid> locallyMovedForeground =
        CreateBoxGrid(-1.4, 1.6, 3, true, 20.0, 1.0, true);
    WriteGrid(directory / "foreground-local-moving1.vtu", locallyMovedForeground);
    WritePvd(directory / "foreground-local-moving.pvd", "foreground0.vtu", 0.0,
             "foreground-local-moving1.vtu", 1.0);
    const fs::path movingLocalManifest = directory / "moving-local-foreground.scpm";
    WriteManifest(movingLocalManifest, "cut_background", "foreground-local-moving.pvd");
    RequireManifestReadRejected(
        movingLocalManifest,
        "foreground local geometry",
        "A moving local foreground with zone/relative velocity was accepted.");

    const fs::path unknownKeyManifest = directory / "unknown-key.scpm";
    WriteManifest(unknownKeyManifest, "cut_background", "foreground.pvd");
    std::ofstream unknownKeyOutput(unknownKeyManifest, std::ios::app);
    Require(static_cast<bool>(unknownKeyOutput), "Cannot append synthetic invalid manifest key.");
    unknownKeyOutput << "zone1_velocity_fram=relative\n";
    unknownKeyOutput.close();
    bool rejected = false;
    try {
        static_cast<void>(io::EnumerateDataSetFrames(unknownKeyManifest));
    } catch (const std::exception& error) {
        rejected = true;
        Require(std::string(error.what()).find("Unknown mesh-set manifest key") !=
                    std::string::npos,
                "An unknown mesh-set manifest key failed for the wrong reason.");
    }
    Require(rejected, "An unknown mesh-set manifest key was silently accepted.");
}

void ValidateManifestReaderAndResampler(const fs::path& root) {
    const fs::path directory = root / "two-zone";
    fs::create_directories(directory);
    vtkSmartPointer<vtkUnstructuredGrid> background =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0);
    vtkSmartPointer<vtkUnstructuredGrid> foreground =
        CreateBoxGrid(-1.5, 1.5, 3, true, 20.0, 1.0, true);
    WriteGrid(directory / "background0.vtu", background);
    WriteGrid(directory / "background1.vtu", background);
    WriteGrid(directory / "foreground0.vtu", foreground);
    WriteGrid(directory / "foreground1.vtu", foreground);
    ValidateStaticMeshKinds(directory, background);
    ValidatePhysicalVolumeBackgroundSelection(directory);
    ValidateNekAutomaticBackgroundSelection(directory);
    ValidateVtkHdfNativeTimes(directory);
    ValidateScalarOnlyManifest(directory, background);
    ValidateCellVelocityResampler(directory);
    ValidateVelocityTupleGuards(directory);
    ValidateBlendFieldGuards(directory);
    WritePvd(directory / "background.pvd", "background0.vtu", 0.0,
             "background1.vtu", 1.0);
    WritePvd(directory / "foreground.pvd", "foreground0.vtu", 0.0,
             "foreground1.vtu", 1.0);
    ValidateFailClosedMeshSetGuards(directory);

    const fs::path cutManifest = directory / "cut.scpm";
    WriteManifest(cutManifest, "cut_background", "foreground.pvd");
    Require(io::IsMultiZoneMeshManifest(cutManifest) && io::IsSupportedDataSource(cutManifest),
            "Valid .scpm input was not recognized as a data source.");
    const std::vector<io::DataSetFrame> frames = io::EnumerateDataSetFrames(cutManifest);
    Require(frames.size() == 2 &&
                frames.front().kind == io::DataSetFrameKind::MultiZoneManifestTimeStep &&
                frames.front().multiZoneBackgroundIndex == 0,
            "Two-zone manifest did not enumerate stable canonical frames.");
    vtkSmartPointer<vtkDataSet> cut = io::ReadDataSet(frames[1]);
    Require(cut->GetNumberOfPoints() == background->GetNumberOfPoints(),
            "Canonical reader did not retain the larger background topology.");
    RequireCoveredVelocity(cut, false);

    const fs::path staticBackgroundTimeManifest = directory / "static-background-time.scpm";
    WriteManifest(staticBackgroundTimeManifest, "cut_background", "foreground.pvd",
                  "background0.vtu");
    const std::vector<io::DataSetFrame> staticBackgroundTimeFrames =
        io::EnumerateDataSetFrames(staticBackgroundTimeManifest);
    vtkSmartPointer<vtkDataSet> staticBackgroundTime = io::ReadDataSet(
        staticBackgroundTimeFrames.at(1));
    vtkDataArray* logicalTime = staticBackgroundTime->GetFieldData()->GetArray("TimeValue");
    Require(logicalTime != nullptr && logicalTime->GetNumberOfTuples() == 1 &&
                NearlyEqual(logicalTime->GetComponent(0, 0),
                            staticBackgroundTimeFrames.at(1).time),
            "A composed mesh-set frame did not expose its logical frame time in FieldData.");

    vtkDataArray* sourceInterface = cut->GetPointData()->GetArray(
        io::kMultiZoneInterfaceArray);
    Require(sourceInterface != nullptr,
            "Synthetic two-zone source is missing its point interface diagnostic.");
    const vtkIdType sourceInterfaceTupleCount = sourceInterface->GetNumberOfTuples();
    std::vector<double> sourceInterfaceValues(
        static_cast<std::size_t>(sourceInterfaceTupleCount));
    bool sawInterface = false;
    for (vtkIdType pointId = 0; pointId < sourceInterfaceTupleCount; ++pointId) {
        sourceInterfaceValues[static_cast<std::size_t>(pointId)] =
            sourceInterface->GetComponent(pointId, 0);
        sawInterface = sawInterface || sourceInterfaceValues[static_cast<std::size_t>(pointId)] != 0.0;
    }
    Require(sawInterface, "Synthetic two-zone source did not create an interface stencil.");
    vtkSmartPointer<vtkDataSet> preparedCut = resample::PreparePointData(cut);
    vtkDataArray* preparedInterface = preparedCut->GetPointData()->GetArray(
        io::kMultiZoneInterfaceArray);
    Require(preparedInterface != nullptr &&
                preparedInterface->GetNumberOfTuples() == sourceInterfaceTupleCount,
            "PreparePointData dropped the canonical point interface diagnostic.");
    for (vtkIdType pointId = 0; pointId < preparedInterface->GetNumberOfTuples(); ++pointId) {
        Require(NearlyEqual(preparedInterface->GetComponent(pointId, 0),
                            sourceInterfaceValues[static_cast<std::size_t>(pointId)]),
                "PreparePointData overwrote the canonical point interface diagnostic.");
    }
    vtkDataArray* preparedCutMask = preparedCut->GetPointData()->GetArray(
        io::kMultiZoneBackgroundCutArray);
    Require(preparedCutMask != nullptr,
            "PreparePointData did not conservatively promote the cell cut diagnostic.");
    bool sawCut = false;
    for (vtkIdType pointId = 0; pointId < preparedCutMask->GetNumberOfTuples(); ++pointId) {
        sawCut = sawCut || preparedCutMask->GetComponent(pointId, 0) != 0.0;
    }
    Require(sawCut, "PreparePointData lost all covered-cell cut diagnostics.");

    const std::vector<resample::FrameInfo> resampleFrames =
        resample::EnumerateDataSetFrames(cutManifest);
    Require(resampleFrames.size() == 2 &&
                resampleFrames.front().kind ==
                    resample::DataSetFrameKind::MultiZoneManifestTimeStep &&
                resampleFrames.front().multiZoneBackgroundIndex == 0,
            "Resampler wrapper did not preserve the two-zone frame kind.");
    const std::array<double, 6> resampleBounds =
        resample::ReadDataSetBounds(resampleFrames[1]);
    Require(NearlyEqual(resampleBounds[0], 0.0) && NearlyEqual(resampleBounds[1], 4.0) &&
                NearlyEqual(resampleBounds[2], 0.0) && NearlyEqual(resampleBounds[3], 4.0) &&
                NearlyEqual(resampleBounds[4], 0.0) && NearlyEqual(resampleBounds[5], 4.0),
            "Resampler bounds path did not preserve the canonical background bounds.");
    vtkSmartPointer<vtkDataSet> resampleData = resample::ReadDataSet(resampleFrames[1]);
    RequireCoveredVelocity(resampleData, false);

    const fs::path blendManifest = directory / "blend.scpm";
    WriteManifest(blendManifest, "blend_50_50", "foreground.pvd");
    const std::vector<io::DataSetFrame> blendFrames =
        io::EnumerateDataSetFrames(blendManifest);
    RequireCoveredVelocity(io::ReadDataSet(blendFrames[1]), true);

    WritePvd(directory / "foreground-bad-times.pvd", "foreground0.vtu", 0.0,
             "foreground1.vtu", 2.0);
    const fs::path badManifest = directory / "bad-times.scpm";
    WriteManifest(badManifest, "cut_background", "foreground-bad-times.pvd");
    bool rejected = false;
    try {
        static_cast<void>(io::EnumerateDataSetFrames(badManifest));
    } catch (const std::exception&) {
        rejected = true;
    }
    Require(rejected, "Mismatched dynamic zone times were not rejected.");

    const fs::path wrongKindManifest = directory / "wrong-kind.scpm";
    WriteManifest(wrongKindManifest, "cut_background", "foreground.pvd",
                  "background.pvd", "static");
    rejected = false;
    try {
        static_cast<void>(io::EnumerateDataSetFrames(wrongKindManifest));
    } catch (const std::exception&) {
        rejected = true;
    }
    Require(rejected, "A moving rigid zone was accepted with a static mesh kind.");

    vtkSmartPointer<vtkUnstructuredGrid> movedBackground =
        CreateBoxGrid(0.1, 4.1, 4, false, 5.0, 10.0);
    WriteGrid(directory / "background-moving1.vtu", movedBackground);
    WritePvd(directory / "background-moving.pvd", "background0.vtu", 0.0,
             "background-moving1.vtu", 1.0);
    const fs::path movingBackgroundManifest = directory / "moving-background.scpm";
    WriteManifest(movingBackgroundManifest, "cut_background", "foreground.pvd",
                  "background-moving.pvd");
    rejected = false;
    try {
        const std::vector<io::DataSetFrame> movingFrames =
            io::EnumerateDataSetFrames(movingBackgroundManifest);
        static_cast<void>(io::ReadDataSet(movingFrames[1]));
    } catch (const std::exception&) {
        rejected = true;
    }
    Require(rejected, "A changing selected background geometry was not rejected.");

    vtkSmartPointer<vtkUnstructuredGrid> topologyChangedBackground =
        CreateBoxGrid(0.0, 4.0, 4, false, 5.0, 10.0, false, true);
    WriteGrid(directory / "background-topology1.vtu", topologyChangedBackground);
    WritePvd(directory / "background-topology.pvd", "background0.vtu", 0.0,
             "background-topology1.vtu", 1.0);
    const fs::path topologyBackgroundManifest = directory / "topology-background.scpm";
    WriteManifest(topologyBackgroundManifest, "cut_background", "foreground.pvd",
                  "background-topology.pvd");
    rejected = false;
    try {
        const std::vector<io::DataSetFrame> topologyFrames =
            io::EnumerateDataSetFrames(topologyBackgroundManifest);
        static_cast<void>(io::ReadDataSet(topologyFrames[1]));
    } catch (const std::exception&) {
        rejected = true;
    }
    Require(rejected, "A selected background with changed topology was not rejected.");

    vtkSmartPointer<vtkDataSet> ordinary = io::ReadDataSet(directory / "background0.vtu");
    Require(ordinary->GetPointData()->GetArray(io::kMultiZoneOwnershipArray) == nullptr &&
                ordinary->GetPointData()->GetArray("Velocity") != nullptr,
            "Ordinary single-grid reader incorrectly took the two-zone path.");
}

}  // namespace

int main() {
    try {
        const fs::path root = fs::path(SCP_RESAMPLER_IO_TEST_OUTPUT_ROOT) /
                              "multizone";
        fs::remove_all(root);
        fs::create_directories(root);
        ValidateManifestReaderAndResampler(root);
        std::cout << "Two-zone manifest reader/resampler test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Two-zone manifest reader/resampler test failed: " << error.what()
                  << '\n';
        return 1;
    }
}
