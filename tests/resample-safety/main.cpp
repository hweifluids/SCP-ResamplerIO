#include "streamcenterplus_resample/Core.h"
#include "Integrity.h"
#include "../EnSightGoldFixture.h"

#include <vtkDataArray.h>
#include <vtkCellType.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace resample = streamcenterplus::resample;

namespace {

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::string ReadText(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    Require(static_cast<bool>(input), "Cannot read test file: " + path.string());
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
}

void WriteText(const fs::path& path, const std::string& text) {
    fs::create_directories(path.parent_path());
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    Require(static_cast<bool>(output), "Cannot create test file: " + path.string());
    output.write(text.data(), static_cast<std::streamsize>(text.size()));
    output.close();
    Require(!output.fail(), "Cannot close test file: " + path.string());
}

std::string ExpectThrows(const std::function<void()>& operation,
                         const std::string& expectedText) {
    try {
        operation();
    } catch (const std::exception& exception) {
        const std::string message = exception.what();
        Require(message.find(expectedText) != std::string::npos,
                "Exception did not contain '" + expectedText + "': " + message);
        return message;
    }
    throw std::runtime_error("Expected operation to throw: " + expectedText);
}

void WriteInputImage(const fs::path& path, double value) {
    vtkNew<vtkImageData> image;
    image->SetDimensions(2, 2, 2);
    image->SetOrigin(0.0, 0.0, 0.0);
    image->SetSpacing(1.0, 1.0, 1.0);
    vtkNew<vtkDoubleArray> values;
    values->SetName("Signal");
    values->SetNumberOfComponents(1);
    values->SetNumberOfTuples(image->GetNumberOfPoints());
    values->FillValue(value);
    image->GetPointData()->AddArray(values);
    image->GetPointData()->SetActiveScalars("Signal");
    resample::WriteVti(path, image);
}

double ReadFirstSignal(const fs::path& path) {
    vtkSmartPointer<vtkDataSet> data = resample::ReadDataSet(path);
    vtkDataArray* values = data->GetPointData()->GetArray("Signal");
    Require(values != nullptr && values->GetNumberOfTuples() > 0,
            "Resampled output has no Signal array: " + path.string());
    return values->GetComponent(0, 0);
}

resample::Parameters BasicParameters(const fs::path& input,
                                     const fs::path& output,
                                     const std::string& basename) {
    resample::Parameters params;
    params.deckFile = output / "test.resample.scpo";
    params.inputPath = input;
    params.outputDir = output;
    params.outputBasename = basename;
    params.dimensions = {3, 3, 3};
    params.boundsMode = resample::BoundsMode::Global;
    params.method = resample::InterpolationMethod::Linear;
    params.linear.computeTolerance = false;
    params.linear.tolerance = 1.0e-8;
    return params;
}

void TestStrictParameters(const fs::path& root) {
    const fs::path decks = root / "parameters";
    fs::create_directories(decks);
    auto deck = [&](const std::string& name, const std::string& body) {
        const fs::path path = decks / (name + ".resample.scpo");
        WriteText(path, body);
        return path;
    };
    const std::string input = "input_path = dummy.vti\n";

    ExpectThrows([&] { (void)resample::LoadParameters(deck("missing", "dimensions=[2,2,2]\n")); },
                 "Missing required input path");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("bool", input + "enabled=perhaps\n")); },
                 "Invalid boolean");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("integer", input + "dimensions=[2x,2,2]\n")); },
                 "complete in-range integer");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("double", input + "physical_dt=1junk\n")); },
                 "complete finite numeric");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("nan", input + "idw.null_value=nan\n")); },
                 "complete finite numeric");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck("empty-number", input + "linear.tolerance=\n"));
    }, "finite number, not an empty value");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck("empty-output", input + "output_dir=\n"));
    }, "output_dir was provided but is empty");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("linear", input + "linear.tolerance=-1\n")); },
                 "Linear tolerance");
    ExpectThrows([&] { (void)resample::LoadParameters(deck("idw", input + "idw.neighbors=0\n")); },
                 "IDW neighbors");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck(
            "kriging", input + "kriging.neighbors=2\nkriging.min_neighbors=3\n"));
    }, "min_neighbors cannot exceed");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck(
            "kriging-radius", input + "kriging.neighbor_mode=radius\n"
                                      "kriging.neighbors=2\nkriging.min_neighbors=3\n"));
    }, "min_neighbors cannot exceed");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck(
            "kriging-neighbor-limit", input + "kriging.neighbors=65535\n"));
    }, "Kriging neighbors cannot exceed");
    vtkNew<vtkImageData> directKrigingInput;
    resample::Parameters directKriging;
    directKriging.method = resample::InterpolationMethod::Kriging;
    directKriging.kriging.neighbors = 65535;
    ExpectThrows([&] {
        (void)resample::ResampleFrame(directKrigingInput, directKriging,
                                      {0.0, 1.0, 0.0, 1.0, 0.0, 1.0});
    }, "Kriging neighbors must be between");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck("traversal", input + "output_basename=../escape\n"));
    }, "output_basename");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck("separator", input + "output_basename=sub/name\n"));
    }, "output_basename");
    ExpectThrows([&] {
        (void)resample::LoadParameters(deck("absolute", input + "output_basename=/escape\n"));
    }, "output_basename");

    const resample::Parameters fields = resample::LoadParameters(deck(
        "interpolation-fields",
        input +
            "linear.compute_tolerance=false\n"
            "linear.tolerance=0.00025\n"
            "linear.snap_to_cell_with_closest_point=true\n"
            "idw.neighbor_mode=radius\n"
            "idw.neighbors=7\n"
            "idw.radius=auto\n"
            "idw.power=3.5\n"
            "idw.null_strategy=null_value\n"
            "idw.null_value=-4.25\n"
            "kriging.variogram=gaussian\n"
            "kriging.neighbor_mode=radius\n"
            "kriging.neighbors=9\n"
            "kriging.min_neighbors=3\n"
            "kriging.radius=2.5\n"
            "kriging.range=auto\n"
            "kriging.sill=4.5\n"
            "kriging.nugget=0.25\n"
            "kriging.regularization=0.000001\n"
            "kriging.null_value=-8.5\n"
            "kriging.fallback=closest\n"
            "kriging.fallback_idw_power=1.75\n"));
    Require(!fields.linear.computeTolerance &&
                std::abs(fields.linear.tolerance - 0.00025) < 1.0e-15 &&
                fields.linear.snapToCellWithClosestPoint &&
                fields.idw.neighborMode == resample::NeighborMode::Radius &&
                fields.idw.neighbors == 7 && fields.idw.radius == 0.0 &&
                std::abs(fields.idw.power - 3.5) < 1.0e-15 &&
                fields.idw.nullStrategy == resample::NullStrategy::NullValue &&
                std::abs(fields.idw.nullValue + 4.25) < 1.0e-15 &&
                fields.kriging.variogram == resample::KrigingVariogram::Gaussian &&
                fields.kriging.neighborMode == resample::NeighborMode::Radius &&
                fields.kriging.neighbors == 9 && fields.kriging.minNeighbors == 3 &&
                std::abs(fields.kriging.radius - 2.5) < 1.0e-15 &&
                fields.kriging.range == 0.0 &&
                std::abs(fields.kriging.sill - 4.5) < 1.0e-15 &&
                std::abs(fields.kriging.nugget - 0.25) < 1.0e-15 &&
                std::abs(fields.kriging.regularization - 0.000001) < 1.0e-18 &&
                std::abs(fields.kriging.nullValue + 8.5) < 1.0e-15 &&
                fields.kriging.fallback == resample::KrigingFallback::Closest &&
                std::abs(fields.kriging.fallbackIdwPower - 1.75) < 1.0e-15,
            "Linear/IDW/Kriging parameter fields were not loaded exactly.");
}

void TestPvdAttributeHandling(const fs::path& root) {
    const fs::path directory = root / "pvd-attributes";
    const fs::path frame = directory / "frame&name.vti";
    WriteInputImage(frame, 2.0);

    const fs::path generatedPvd = directory / "generated.pvd";
    resample::WritePvd(generatedPvd, {{2.5, frame.filename()}});
    Require(ReadText(generatedPvd).find("frame&amp;name.vti") != std::string::npos,
            "PVD writer did not XML-escape the generated frame name.");
    const std::vector<resample::FrameInfo> generatedFrames =
        resample::EnumerateDataSetFrames(generatedPvd);
    Require(generatedFrames.size() == 1 &&
                generatedFrames.front().path.lexically_normal() == frame.lexically_normal() &&
                std::abs(generatedFrames.front().time - 2.5) < 1.0e-12,
            "PVD reader did not round-trip an XML-escaped frame name and timestep.");

    for (const std::string& timestep : {std::string("nan"), std::string("1.5ms")}) {
        const fs::path invalidPvd = directory / ("invalid-time-" + timestep + ".pvd");
        WriteText(invalidPvd,
                  "<?xml version=\"1.0\"?>\n"
                  "<VTKFile type=\"Collection\"><Collection>\n"
                  "<DataSet timestep=\"" + timestep + "\" file=\"frame&amp;name.vti\"/>\n"
                  "</Collection></VTKFile>\n");
        ExpectThrows([&] { (void)resample::EnumerateDataSetFrames(invalidPvd); },
                     "Invalid PVD timestep");
    }

    const fs::path commentedPvd = directory / "commented.pvd";
    WriteText(commentedPvd,
              "<?xml version=\"1.0\"?>\n"
              "<VTKFile type=\"Collection\"><Collection>\n"
              "<!-- <DataSet timestep=\"99\" file=\"missing.vti\"/> -->\n"
              "<![CDATA[<DataSet timestep=\"98\" file=\"also-missing.vti\"/>]]>\n"
              "<DataSet timestep=\"3.5\" file=\"frame&amp;name.vti\"/>\n"
              "</Collection></VTKFile>\n");
    const std::vector<resample::FrameInfo> commentedFrames =
        resample::EnumerateDataSetFrames(commentedPvd);
    Require(commentedFrames.size() == 1 &&
                commentedFrames.front().path.lexically_normal() == frame.lexically_normal() &&
                std::abs(commentedFrames.front().time - 3.5) < 1.0e-12,
            "PVD parser treated a comment or CDATA payload as a DataSet entry.");

    const fs::path collectionDirectory = directory / "directory-collection";
    const fs::path collectionFrame = collectionDirectory / "frame0.vti";
    WriteInputImage(collectionFrame, 4.0);
    resample::WritePvd(collectionDirectory / "series.pvd", {{7.5, collectionFrame.filename()}});
    const std::vector<resample::FrameInfo> directoryFrames =
        resample::EnumerateDataSetFrames(collectionDirectory);
    Require(directoryFrames.size() == 1 &&
                directoryFrames.front().path.lexically_normal() == collectionFrame.lexically_normal() &&
                std::abs(directoryFrames.front().time - 7.5) < 1.0e-12,
            "Directory enumeration did not use its explicit PVD collection.");

    const fs::path cyclicPvd = directory / "cyclic.pvd";
    WriteText(cyclicPvd,
              "<?xml version=\"1.0\"?>\n"
              "<VTKFile type=\"Collection\"><Collection>\n"
              "<DataSet timestep=\"0\" file=\"cyclic.pvd\"/>\n"
              "</Collection></VTKFile>\n");
    ExpectThrows([&] { (void)resample::EnumerateDataSetFrames(cyclicPvd); },
                 "Recursive PVD collection");

    const fs::path longTagDirectory = directory / "long-numeric-tag";
    const fs::path longTagFrame = longTagDirectory /
        "frame999999999999999999999999999999999999999999999999999999.vti";
    WriteInputImage(longTagFrame, 5.0);
    const std::vector<resample::FrameInfo> longTagFrames =
        resample::EnumerateDataSetFrames(longTagDirectory);
    Require(longTagFrames.size() == 1 &&
                longTagFrames.front().path.lexically_normal() == longTagFrame.lexically_normal(),
            "A long numeric filename prevented deterministic folder enumeration.");
}

void WriteNekHeader(const fs::path& fieldPath, const std::string& header) {
    Require(header.size() <= 132, "Synthetic Nek header exceeds 132 bytes.");
    fs::create_directories(fieldPath.parent_path());
    std::ofstream output(fieldPath, std::ios::binary | std::ios::trunc);
    Require(static_cast<bool>(output), "Cannot create synthetic Nek field.");
    std::array<char, 132> bytes{};
    bytes.fill(' ');
    std::copy(header.begin(), header.end(), bytes.begin());
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    const float endianTag = 6.54321f;
    output.write(reinterpret_cast<const char*>(&endianTag), sizeof(endianTag));
    output.close();
    Require(!output.fail(), "Cannot close synthetic Nek field.");
}

void ExpectRejectedNekHeader(const fs::path& root,
                             const std::string& label,
                             const std::string& header,
                             const std::string& expected) {
    const fs::path directory = root / "nek" / label;
    const fs::path descriptor = directory / "case.nek5000";
    WriteText(descriptor,
              "filetemplate: field%05d.fld\nfirsttimestep: 0\nnumtimesteps: 1\n");
    WriteNekHeader(directory / "field00000.fld", header);
    ExpectThrows([&] { (void)resample::EnumerateDataSetFrames(descriptor); }, expected);
}

void TestMaliciousNekHeaders(const fs::path& root) {
    ExpectRejectedNekHeader(root, "huge-dimensions",
                            "#std 8 2147483647 2147483647 2 1 1 0 0 0 1 X 0 F",
                            "resource limits");
    ExpectRejectedNekHeader(root, "huge-file-count",
                            "#std 8 2 2 2 1 1 0 0 0 2147483647 X 0 F",
                            "resource limits");
    ExpectRejectedNekHeader(root, "huge-stored-values",
                            "#std 8 2 2 2 10000000 10000000 0 0 0 1 S20 0 F",
                            "stored-value resource limit");
    ExpectRejectedNekHeader(root, "trailing-token",
                            "#std 8 2junk 2 2 1 1 0 0 0 1 X 0 F",
                            "Invalid Nek5000 nx");
    ExpectRejectedNekHeader(root, "truncated",
                            "#std 8 2 2 2 1 1 0 0 0 1 XUP 0 F",
                            "truncated");
}

void TestInvalidInputPreservesOutputs(const fs::path& root) {
    const fs::path output = root / "invalid-input-output";
    fs::create_directories(output);
    WriteText(output / "preserve.pvd", "pvd-sentinel");
    WriteText(output / "unrelated.vti", "vti-sentinel");
    resample::Parameters params = BasicParameters(root / "missing-input", output, "preserve");
    ExpectThrows([&] { (void)resample::RunSolver(params); }, "does not exist");
    Require(ReadText(output / "preserve.pvd") == "pvd-sentinel",
            "Invalid input changed the pre-existing namespace PVD.");
    Require(ReadText(output / "unrelated.vti") == "vti-sentinel",
            "Invalid input removed an unrelated VTI.");

    const fs::path validInput = root / "valid-input.vti";
    WriteInputImage(validInput, 1.0);
    params = BasicParameters(validInput, output, "preserve");
    params.dimensions = {0, 3, 3};
    ExpectThrows([&] { (void)resample::RunSolver(params); }, "dimensions");
    Require(ReadText(output / "preserve.pvd") == "pvd-sentinel" &&
                ReadText(output / "unrelated.vti") == "vti-sentinel",
            "Invalid executable parameters changed pre-existing outputs.");

    const fs::path collisionOutput = root / "namespace-collision";
    const fs::path collisionEntry = collisionOutput / "collision_f000000.vti";
    fs::create_directories(collisionEntry);
    WriteText(collisionEntry / "sentinel.txt", "keep-collision");
    params = BasicParameters(validInput, collisionOutput, "collision");
    ExpectThrows([&] { (void)resample::RunSolver(params); },
                 "namespace already contains files");
    Require(ReadText(collisionEntry / "sentinel.txt") == "keep-collision",
            "An unowned non-file namespace collision was removed or overwritten.");

    const fs::path caseCollisionOutput = root / "namespace-case-collision";
    fs::create_directories(caseCollisionOutput);
    WriteText(caseCollisionOutput / "COLLISION.PVD", "keep-case-collision");
    params = BasicParameters(validInput, caseCollisionOutput, "collision");
    ExpectThrows([&] { (void)resample::RunSolver(params); },
                 "aliases an existing filesystem entry");
    Require(ReadText(caseCollisionOutput / "COLLISION.PVD") == "keep-case-collision",
            "A case-insensitive unowned namespace collision was overwritten.");
}

void TestContentBoundResumeAndOwnedCleanup(const fs::path& root) {
    const fs::path input = root / "resume-input";
    const fs::path output = root / "resume-output";
    fs::create_directories(input);
    fs::create_directories(output);
    WriteInputImage(input / "frame0.vti", 1.0);
    WriteInputImage(input / "frame1.vti", 2.0);
    WriteText(output / "other.pvd", "unrelated-pvd");
    WriteText(output / "other.vti", "unrelated-vti");

    resample::Parameters params = BasicParameters(input, output, "cache");
    Require(resample::RunSolver(params) == 0, "Initial resample run failed.");
    Require(std::abs(ReadFirstSignal(output / "cache_f000000.vti") - 1.0) < 1.0e-12,
            "Initial frame did not contain the first source content.");
    Require(ReadText(output / "other.pvd") == "unrelated-pvd" &&
                ReadText(output / "other.vti") == "unrelated-vti",
            "Fresh run changed unrelated VTI/PVD outputs.");

    // The grid, dimensions, bounds, timestamps and array schema stay the same;
    // only the ordered source content changes. A geometry-only resume check
    // would incorrectly keep frame zero.
    WriteInputImage(input / "frame0.vti", 9.0);
    params.solverStatus = "error";
    Require(resample::RunSolver(params) == 0, "Content-change resume run failed.");
    Require(std::abs(ReadFirstSignal(output / "cache_f000000.vti") - 9.0) < 1.0e-12,
            "Resume reused stale same-grid input content.");
    Require(ReadText(output / "other.pvd") == "unrelated-pvd" &&
                ReadText(output / "other.vti") == "unrelated-vti",
            "Resume cleanup changed unrelated VTI/PVD outputs.");
    const std::string manifest = ReadText(output / "cache.scp-resample-manifest");
    Require(manifest.find("state=complete") != std::string::npos &&
                manifest.find("content_sha256=") != std::string::npos &&
                manifest.find("schema_sha256=") != std::string::npos &&
                manifest.find("grid_sha256=") != std::string::npos &&
                manifest.find("output.0.sha256=") != std::string::npos,
            "Resume manifest does not bind required input/output identities.");

    // Keep the cached image type, dimensions and bounds valid while changing
    // only its payload. The manifest checksum must reject this otherwise
    // geometry-compatible output.
    vtkSmartPointer<vtkDataSet> tamperedData =
        resample::ReadDataSet(output / "cache_f000000.vti");
    vtkImageData* tamperedImage = vtkImageData::SafeDownCast(tamperedData);
    Require(tamperedImage != nullptr, "Cached output is not vtkImageData before tamper test.");
    vtkDataArray* tamperedSignal = tamperedImage->GetPointData()->GetArray("Signal");
    Require(tamperedSignal != nullptr, "Cached output has no Signal array before tamper test.");
    for (vtkIdType tuple = 0; tuple < tamperedSignal->GetNumberOfTuples(); ++tuple) {
        tamperedSignal->SetComponent(tuple, 0, 123.0);
    }
    resample::WriteVti(output / "cache_f000000.vti", tamperedImage);
    Require(resample::RunSolver(params) == 0, "Output-checksum resume run failed.");
    Require(std::abs(ReadFirstSignal(output / "cache_f000000.vti") - 9.0) < 1.0e-12,
            "Resume reused a payload-tampered output with compatible geometry.");
}

void TestAtomicWriteFailure(const fs::path& root) {
    const fs::path directoryTarget = root / "atomic-failure" / "blocked.pvd";
    fs::create_directories(directoryTarget);
    WriteText(directoryTarget / "sentinel.txt", "keep");
    ExpectThrows([&] { resample::WritePvd(directoryTarget, {}); }, "Atomic file replacement failed");
    Require(ReadText(directoryTarget / "sentinel.txt") == "keep",
            "Failed atomic PVD publish damaged the existing destination directory.");
    for (const fs::directory_entry& entry : fs::directory_iterator(directoryTarget.parent_path())) {
        Require(entry.path().filename().string().find("blocked.pvd.tmp.") != 0,
                "Failed atomic publish left a temporary PVD file behind.");
    }

    const fs::path summaryTarget = root / "atomic-failure" / "blocked-summary.txt";
    fs::create_directories(summaryTarget);
    WriteText(summaryTarget / "sentinel.txt", "keep-summary");
    ExpectThrows([&] { resample::WriteTextFile(summaryTarget, "replacement"); },
                 "Atomic file replacement failed");
    Require(ReadText(summaryTarget / "sentinel.txt") == "keep-summary",
            "Failed atomic summary publish damaged the existing destination directory.");
    for (const fs::directory_entry& entry : fs::directory_iterator(summaryTarget.parent_path())) {
        Require(entry.path().filename().string().find("blocked-summary.txt.tmp.") != 0,
                "Failed atomic publish left a temporary summary file behind.");
    }
}

void TestEnSightGoldInputPaths(const fs::path& root) {
    streamcenterplus::test::EnSightGoldFixtureOptions fixtureOptions;
    fixtureOptions.includeBoundaryPart = true;
    const fs::path descriptor = streamcenterplus::test::CreateTransientEnSightGoldFixture(
        root / "encas-input", "region.EnCaS", fixtureOptions);
    const std::vector<resample::FrameInfo> frames =
        resample::EnumerateDataSetFrames(descriptor, 7.0);
    const std::vector<resample::FrameInfo> directoryFrames =
        resample::EnumerateDataSetFrames(descriptor.parent_path(), 7.0);
    Require(frames.size() == 3 && directoryFrames.size() == 3,
            "Resample file/directory inputs did not enumerate ENCAS native frames.");
    for (std::size_t index = 0; index < frames.size(); ++index) {
        Require(frames[index].kind == resample::DataSetFrameKind::EnSightGoldTimeStep &&
                    frames[index].readerTime.has_value() &&
                    std::abs(frames[index].time -
                             streamcenterplus::test::kEnSightGoldTimes[index]) < 1.0e-6 &&
                    directoryFrames[index].time == frames[index].time,
                "Resample replaced the EnSight native frame metadata.");
    }

    const std::array<double, 6> expectedBounds = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    const resample::FrameInfo plainFrame{descriptor, 0.0};
    Require(resample::ReadDataSetBounds(descriptor) == expectedBounds &&
                resample::ReadDataSetBounds(plainFrame) == expectedBounds &&
                resample::ReadDataSetBounds(frames.back()) == expectedBounds,
            "Resample bounds reads did not dispatch ENCAS/static geometry correctly.");
    const auto requirePressure = [&](vtkDataSet* data, double expected) {
        Require(data != nullptr, "Resample did not produce an EnSight data set.");
        vtkDataArray* pressure = data->GetPointData()->GetArray("pressure");
        Require(pressure != nullptr && pressure->GetNumberOfTuples() > 0 &&
                    std::abs(pressure->GetComponent(0, 0) - expected) < 1.0e-10,
                "Resample read the wrong EnSight variable time step.");
    };
    vtkSmartPointer<vtkDataSet> first = resample::ReadDataSet(descriptor);
    Require(first->GetNumberOfPoints() == 8 && first->GetNumberOfCells() == 1,
            "Resample retained separate EnSight boundary parts beside volume parts.");
    requirePressure(first, 0.0);
    requirePressure(resample::ReadDataSet(plainFrame), 0.0);
    requirePressure(resample::ReadDataSet(frames.back()), 20.0);

    const fs::path output = root / "encas-output";
    const resample::Parameters params = BasicParameters(descriptor, output, "encas");
    Require(resample::RunSolver(params) == 0,
            "Transient EnSight Gold resampling failed.");
    const std::vector<resample::FrameInfo> resultFrames =
        resample::EnumerateDataSetFrames(output / "encas.pvd");
    Require(resultFrames.size() == 3,
            "EnSight resampling did not publish all variable frames.");
    for (std::size_t index = 0; index < resultFrames.size(); ++index) {
        Require(std::abs(resultFrames[index].time - frames[index].time) < 1.0e-10,
                "Resampled PVD changed the native EnSight timestamps.");
        requirePressure(resample::ReadDataSet(resultFrames[index]), 10.0 * index);
    }
}

void TestEnSightTimeProvenance(const fs::path& root) {
    const fs::path descriptor =
        streamcenterplus::test::CreateQuantizedEnSightGoldFixture(root / "quantized-encas");
    resample::Parameters params = BasicParameters(descriptor, root / "quantized-output", "quantized");
    Require(resample::RunSolver(params) == 0, "Quantized EnSight resampling failed.");
    const fs::path pvd = params.outputDir / "quantized.pvd";
    const auto frames = resample::EnumerateDataSetFrames(pvd);
    Require(frames.size() == 8 && frames.front().time == 4.25700 &&
                frames[2].time == 4.25901 && frames.back().time == 4.26401 &&
                frames.front().ensightTime.has_value(),
            "Resample output dropped original EnSight times or validated time provenance.");
    Require(ReadText(pvd).find(streamcenterplus::ensight_time::kTimeQuantumAttribute) !=
                std::string::npos,
            "Quantized EnSight Resample PVD did not publish its decimal precision.");
    params.solverStatus = "error";
    Require(resample::RunSolver(params) == 0 &&
                resample::EnumerateDataSetFrames(pvd).front().ensightTime.has_value(),
            "Resample resume dropped EnSight time provenance.");
    params.outputDir = root / "quantized-cropped-output";
    params.outputBasename = "cropped";
    params.solverStatus = "unsolved";
    params.frameRange.start = 2;
    params.frameRange.end = 3;
    Require(resample::RunSolver(params) == 0, "Two-frame EnSight selection failed to resample.");
    const fs::path cropped = params.outputDir / "cropped.pvd";
    Require(ReadText(cropped).find(streamcenterplus::ensight_time::kTimeQuantumAttribute) ==
                std::string::npos &&
                resample::EnumerateDataSetFrames(cropped).front().time == 4.25901,
            "Strict two-frame output retained stale EnSight approximation metadata or changed times.");
}

void WriteCompositeVtuBlock(const fs::path& path, int frame, bool boundary) {
    const std::size_t pointCount = boundary ? 4 : 8;
    std::ostringstream xml;
    xml << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
        << "<UnstructuredGrid><Piece NumberOfPoints=\"" << pointCount
        << "\" NumberOfCells=\"1\"><PointData>\n";
    const auto scalar = [&](const char* name, int component) {
        xml << "<DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">";
        for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
            const auto& point = streamcenterplus::test::kEnSightGoldPoints[pointIndex];
            xml << (component < 0 ? 10.0 * frame + point[0] + 2.0 * point[1] + 3.0 * point[2]
                                 : 100.0 * frame + 10.0 * component + point[component]) << ' ';
        }
        xml << "</DataArray>\n";
    };
    scalar("pressure", -1);
    scalar("x_velocity", 0);
    scalar("y_velocity", 1);
    scalar("z_velocity", 2);
    xml << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">";
    for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
        for (int component = 0; component < 3; ++component) {
            xml << 100.0 * frame + 10.0 * component +
                       streamcenterplus::test::kEnSightGoldPoints[pointIndex][component] << ' ';
        }
    }
    xml << "</DataArray></PointData><CellData/><Points>"
        << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
    for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
        for (double coordinate : streamcenterplus::test::kEnSightGoldPoints[pointIndex]) {
            xml << coordinate << ' ';
        }
    }
    xml << "</DataArray></Points><Cells>"
        << "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">";
    for (std::size_t pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
        xml << pointIndex << ' ';
    }
    xml << "</DataArray><DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"
        << pointCount << "</DataArray>"
        << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
        << (boundary ? VTK_QUAD : VTK_HEXAHEDRON)
        << "</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>\n";
    WriteText(path, xml.str());
}

void TestPvdVtmCompositeInput(const fs::path& root) {
    const fs::path input = root / "pvd-vtm-input";
    const fs::path collection = input / "series.pvd";
    const std::array<double, 3> times = {4.257, 4.258, 4.25901};
    std::vector<std::pair<double, fs::path>> entries;
    for (int frame = 0; frame < 3; ++frame) {
        const std::string prefix = "frame_" + std::to_string(frame);
        const fs::path volumeName = prefix + "_volume.vtu";
        const fs::path boundaryName = prefix + "_boundary.vtu";
        const fs::path compositeName = prefix + ".VtM";
        WriteCompositeVtuBlock(input / volumeName, frame, false);
        WriteCompositeVtuBlock(input / boundaryName, frame, true);
        WriteText(input / compositeName,
            "<?xml version=\"1.0\"?>\n"
            "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">"
            "<vtkMultiBlockDataSet><DataSet index=\"0\" name=\"volume\" file=\"" +
            volumeName.generic_string() + "\"/><DataSet index=\"1\" name=\"boundary\" file=\"" +
            boundaryName.generic_string() + "\"/></vtkMultiBlockDataSet></VTKFile>\n");
        entries.push_back({times[static_cast<std::size_t>(frame)], compositeName});
    }
    resample::WritePvd(collection, entries, 1.0e-5);
    const auto frames = resample::EnumerateDataSetFrames(collection);
    Require(frames.size() == 3 && frames.front().ensightTime.has_value(),
            "PVD/VTM discovery lost its frame times or EnSight time provenance.");
    const std::array<double, 6> expectedBounds = {0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    Require(resample::ReadDataSetBounds(frames.front()) == expectedBounds &&
                resample::ReadDataSetBounds(frames.front().path) == expectedBounds,
            "PVD/VTM bounds did not use the shared composite reader.");
    Require(resample::ReadDataSet(frames.front().path)->GetNumberOfCells() == 2,
            "Direct VTM input did not dispatch through the shared composite reader.");
    const std::array<const char*, 3> componentNames = {"x_velocity", "y_velocity", "z_velocity"};
    for (std::size_t frame = 0; frame < frames.size(); ++frame) {
        Require(frames[frame].time == times[frame], "PVD/VTM input changed an original time.");
        vtkSmartPointer<vtkDataSet> data = resample::ReadDataSet(frames[frame]);
        Require(data->GetNumberOfPoints() == 12 && data->GetNumberOfCells() == 2 &&
                    data->GetCellType(0) == VTK_HEXAHEDRON && data->GetCellType(1) == VTK_QUAD,
                "Generic VTM must preserve volume and boundary blocks without EnSight filtering.");
        for (const char* name : {"pressure", "x_velocity", "y_velocity", "z_velocity", "velocity"}) {
            vtkDataArray* array = data->GetPointData()->GetArray(name);
            Require(array != nullptr && array->GetNumberOfTuples() == 12 &&
                        array->GetNumberOfComponents() == (std::string(name) == "velocity" ? 3 : 1),
                    "PVD/VTM assembly dropped or truncated an original field.");
        }
        vtkDataArray* pressure = data->GetPointData()->GetArray("pressure");
        vtkDataArray* velocity = data->GetPointData()->GetArray("velocity");
        for (vtkIdType pointId = 0; pointId < data->GetNumberOfPoints(); ++pointId) {
            double point[3] = {};
            data->GetPoint(pointId, point);
            Require(std::abs(pressure->GetComponent(pointId, 0) -
                             (10.0 * frame + point[0] + 2.0 * point[1] + 3.0 * point[2])) < 1.0e-12,
                    "PVD/VTM pressure values were changed while appending blocks.");
            for (int component = 0; component < 3; ++component) {
                vtkDataArray* scalar = data->GetPointData()->GetArray(componentNames[component]);
                const double expected = 100.0 * frame + 10.0 * component + point[component];
                Require(std::abs(velocity->GetComponent(pointId, component) -
                                 expected) < 1.0e-12 &&
                            std::abs(scalar->GetComponent(pointId, 0) - expected) < 1.0e-12,
                        "PVD/VTM velocity values were changed while appending blocks.");
            }
        }
    }
    const fs::path output = root / "pvd-vtm-output";
    const auto params = BasicParameters(collection, output, "composite");
    Require(resample::RunSolver(params) == 0, "PVD referencing multi-block VTM frames failed to resample.");
    const auto resultFrames = resample::EnumerateDataSetFrames(output / "composite.pvd");
    Require(resultFrames.size() == 3 && resultFrames.front().ensightTime.has_value(),
            "PVD/VTM resampling lost frame count or validated time provenance.");
    for (std::size_t frame = 0; frame < resultFrames.size(); ++frame) {
        vtkSmartPointer<vtkDataSet> image = resample::ReadDataSet(resultFrames[frame]);
        vtkDataArray* pressure = image->GetPointData()->GetArray("pressure");
        vtkDataArray* velocity = image->GetPointData()->GetArray("velocity");
        Require(resultFrames[frame].time == times[frame] && image->GetNumberOfPoints() == 27 &&
                    pressure != nullptr && velocity != nullptr &&
                    std::abs(pressure->GetComponent(0, 0) - 10.0 * frame) < 1.0e-12 &&
                    std::abs(velocity->GetComponent(0, 0) - 100.0 * frame) < 1.0e-12,
                "PVD/VTM output changed original times or resampled field values.");
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2 || argv[1] == nullptr || std::string(argv[1]).empty()) {
            throw std::runtime_error("Usage: scp_resampler_io_safety <fixture-root>");
        }
        const fs::path root = fs::absolute(argv[1]).lexically_normal();
        const bool insideTests = std::any_of(root.begin(), root.end(), [](const fs::path& part) {
            return part == ".tests";
        });
        if (!insideTests || root.filename() != "safety-regression") {
            throw std::runtime_error(
                "Refusing to reset safety fixture outside .tests/.../safety-regression: " +
                root.string());
        }
        fs::remove_all(root);
        fs::create_directories(root);
        Require(streamcenterplus::resample::integrity::Sha256String("abc") ==
                    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
                "SHA-256 implementation failed its known-answer test.");
        const std::string largeKnownAnswer(1000000, 'a');
        Require(streamcenterplus::resample::integrity::Sha256String(largeKnownAnswer) ==
                    "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
                "SHA-256 implementation failed its million-byte known-answer test.");
        TestStrictParameters(root);
        TestPvdAttributeHandling(root);
        TestMaliciousNekHeaders(root);
        TestInvalidInputPreservesOutputs(root);
        TestContentBoundResumeAndOwnedCleanup(root);
        TestAtomicWriteFailure(root);
        TestEnSightGoldInputPaths(root);
        TestEnSightTimeProvenance(root);
        TestPvdVtmCompositeInput(root);
        std::cout << "Resampler safety regression tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Resampler safety regression test failed: " << exception.what() << "\n";
        return 1;
    }
}
