#include "streamcenterplus_resample/Core.h"

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPointData.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using streamcenterplus::resample::DataSetFrameKind;
using streamcenterplus::resample::FrameInfo;
using streamcenterplus::resample::Parameters;

constexpr std::array<int, 3> kDimensions{16, 8, 6};
constexpr int kFirstFrame = 490;
constexpr int kLastFrame = 491;

void Require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::filesystem::path TestRoot() {
    return std::filesystem::path(SCP_RESAMPLER_IO_TEST_OUTPUT_ROOT).lexically_normal();
}

std::filesystem::path DescriptorPath() {
#ifdef _WIN32
    return std::filesystem::path(
        LR"(Z:\Hemisphere\5_Results\2_Re1370\260614_HM1370_500_550_N1000\hemi.nek5000)");
#else
    return TestRoot() / "smoke" / "SCP-ResamplerIO" / "source" / "hemi.nek5000";
#endif
}

std::string FrameFileName(int index) {
    std::ostringstream name;
    name << "nek5000_adapter_f";
    name.width(6);
    name.fill('0');
    name << index << ".vti";
    return name.str();
}

std::string ReadText(const std::filesystem::path& path) {
    std::ifstream input(path);
    Require(static_cast<bool>(input), "Cannot read integration-test output: " + path.string());
    std::ostringstream text;
    text << input.rdbuf();
    return text.str();
}

std::size_t CountOccurrences(const std::string& text, const std::string& token) {
    std::size_t count = 0;
    std::size_t position = 0;
    while ((position = text.find(token, position)) != std::string::npos) {
        ++count;
        position += token.size();
    }
    return count;
}

void ValidateSourceFrame(const FrameInfo& frame) {
    Require(frame.kind == DataSetFrameKind::Nek5000TimeStep,
            "The adapter did not mark the source as a Nek5000 time step.");
    Require(frame.readerTime.has_value() && std::isfinite(*frame.readerTime),
            "The adapter did not retain the Nek5000 reader time.");
    Require(std::abs(frame.time - *frame.readerTime) <= 1.0e-12,
            "The public frame time differs from the Nek5000 reader time.");

    vtkSmartPointer<vtkDataSet> data = streamcenterplus::resample::ReadDataSet(frame);
    Require(data != nullptr && data->GetNumberOfPoints() > 0,
            "The Nek5000 adapter returned an empty source frame.");
    vtkDataArray* velocity = data->GetPointData()->GetArray("Velocity");
    Require(velocity != nullptr && velocity->GetNumberOfComponents() == 3,
            "The Nek5000 adapter did not expose point-associated Velocity(float3).");

    const std::array<double, 6> bounds = streamcenterplus::resample::ReadDataSetBounds(frame);
    for (double value : bounds) {
        Require(std::isfinite(value), "The Nek5000 adapter returned non-finite bounds.");
    }
}

void ValidateOutput(const std::filesystem::path& outputDirectory, int frameIndex) {
    const std::filesystem::path path = outputDirectory / FrameFileName(frameIndex);
    Require(std::filesystem::is_regular_file(path),
            "The generic solver did not write the expected VTI frame: " + path.string());

    vtkSmartPointer<vtkDataSet> data = streamcenterplus::resample::ReadDataSet(path);
    vtkImageData* image = vtkImageData::SafeDownCast(data);
    Require(image != nullptr, "The generic solver output is not vtkImageData.");
    int dimensions[3] = {};
    image->GetDimensions(dimensions);
    for (std::size_t axis = 0; axis < kDimensions.size(); ++axis) {
        Require(dimensions[axis] == kDimensions[axis],
                "The generic solver output dimensions are incorrect.");
    }
    vtkDataArray* velocity = image->GetPointData()->GetArray("Velocity");
    Require(velocity != nullptr && velocity->GetNumberOfComponents() == 3,
            "The generic solver output does not retain Velocity(float3).");
    Require(image->GetPointData()->GetArray("vtkValidPointMask") != nullptr,
            "The generic solver output is missing vtkValidPointMask.");
}

}  // namespace

int main() {
    try {
        const std::filesystem::path descriptor = DescriptorPath();
        if (!std::filesystem::is_regular_file(descriptor)) {
            std::cout << "SKIP: Nek5000 integration source is unavailable at the fixed location: "
                      << descriptor.string() << "\n";
            return 77;
        }

        const std::vector<FrameInfo> advertised =
            streamcenterplus::resample::EnumerateDataSetFrames(descriptor);
        Require(advertised.size() > static_cast<std::size_t>(kLastFrame),
                "The Nek5000 descriptor does not advertise frames 490 and 491.");

        const std::filesystem::path outputDirectory =
            TestRoot() / "nek5000-integration" / "nek5000-general-adapter";
        std::filesystem::remove_all(outputDirectory);
        std::filesystem::create_directories(outputDirectory);

        Parameters params;
        params.inputPath = descriptor;
        params.outputDir = outputDirectory;
        params.outputBasename = "nek5000_adapter";
        params.dimensions = kDimensions;
        params.method = streamcenterplus::resample::InterpolationMethod::VtkResampleToImage;
        params.boundsMode = streamcenterplus::resample::BoundsMode::Global;
        params.frameRange.start = kFirstFrame;
        params.frameRange.end = kLastFrame;

        const streamcenterplus::resample::SequenceInfo selected =
            streamcenterplus::resample::CollectInputSequence(params);
        Require(selected.frames.size() == 2,
                "Frame-range selection did not retain exactly two Nek5000 frames.");
        for (std::size_t index = 0; index < selected.frames.size(); ++index) {
            const FrameInfo& frame = selected.frames[index];
            const FrameInfo& expected = advertised[static_cast<std::size_t>(kFirstFrame) + index];
            Require(std::abs(frame.time - expected.time) <= 1.0e-12,
                    "Frame-range selection changed an advertised Nek5000 time.");
            ValidateSourceFrame(frame);
        }

        Require(streamcenterplus::resample::RunSolver(params) == 0,
                "The generic Resample solver returned a failure code.");
        ValidateOutput(outputDirectory, 0);
        ValidateOutput(outputDirectory, 1);

        const std::filesystem::path pvdPath = outputDirectory / "nek5000_adapter.pvd";
        const std::string pvd = ReadText(pvdPath);
        Require(CountOccurrences(pvd, "<DataSet ") == 2,
                "The PVD collection does not contain exactly two frames.");
        for (const FrameInfo& frame : selected.frames) {
            std::ostringstream time;
            time.precision(17);
            time << frame.time;
            Require(pvd.find("timestep=\"" + time.str() + "\"") != std::string::npos,
                    "The PVD collection does not preserve an advertised Nek5000 time.");
        }

        std::cout << "SCP-ResamplerIO Nek5000 integration passed: frames "
                  << kFirstFrame << "-" << kLastFrame << ", dimensions "
                  << kDimensions[0] << "x" << kDimensions[1] << "x" << kDimensions[2]
                  << ", output=" << outputDirectory.string() << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "SCP-ResamplerIO Nek5000 integration failed: " << ex.what() << "\n";
        return 1;
    }
}
