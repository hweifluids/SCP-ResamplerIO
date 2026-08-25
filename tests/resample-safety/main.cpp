#include "streamcenterplus_resample/Core.h"
#include "Integrity.h"

#include <vtkDataArray.h>
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
        std::cout << "Resampler safety regression tests passed.\n";
        return 0;
    } catch (const std::exception& exception) {
        std::cerr << "Resampler safety regression test failed: " << exception.what() << "\n";
        return 1;
    }
}
