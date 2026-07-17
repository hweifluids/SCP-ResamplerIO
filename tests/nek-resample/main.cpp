#include <streamcenterplus/io/DataSetReader.h>
#include <streamcenterplus_resample/Core.h>

#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace io = streamcenterplus::io;
namespace resample = streamcenterplus::resample;

namespace {

bool NearlyEqual(double lhs, double rhs, double tolerance = 2.0e-10) {
    return std::abs(lhs - rhs) <=
           tolerance * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

vtkSmartPointer<vtkUnstructuredGrid> CreateSpectralElement() {
    const double inner = 1.0 / std::sqrt(5.0);
    const std::array<double, 4> nodes = {-1.0, -inner, inner, 1.0};
    vtkNew<vtkPoints> points;
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(64);
    vtkNew<vtkDoubleArray> velocity;
    velocity->SetName("Velocity");
    velocity->SetNumberOfComponents(3);
    velocity->SetNumberOfTuples(64);
    vtkIdType pointId = 0;
    for (int k = 0; k < 4; ++k) {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i, ++pointId) {
                const double x = 1.0 + nodes[static_cast<std::size_t>(i)];
                const double y = 2.0 + nodes[static_cast<std::size_t>(j)];
                const double z = 3.0 + nodes[static_cast<std::size_t>(k)];
                points->SetPoint(pointId, x, y, z);
                const double tuple[3] = {
                    x * x + 2.0 * y + 3.0 * z,
                    x * y,
                    z * z * z,
                };
                velocity->SetTuple(pointId, tuple);
            }
        }
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);
    grid->GetPointData()->AddArray(velocity);

    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(io::kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(4, 4, 4);
    grid->GetFieldData()->AddArray(dimensions);
    vtkNew<vtkIntArray> elementCount;
    elementCount->SetName(io::kNek5000ElementCountField);
    elementCount->InsertNextValue(1);
    grid->GetFieldData()->AddArray(elementCount);
    return grid;
}

}  // namespace

int main() {
    try {
        vtkSmartPointer<vtkUnstructuredGrid> input = CreateSpectralElement();
        resample::Parameters parameters;
        parameters.method = resample::InterpolationMethod::Linear;
        parameters.dimensions = {5, 5, 5};
        const std::array<double, 6> bounds = {0.0, 2.0, 1.0, 3.0, 2.0, 4.0};
        vtkSmartPointer<vtkImageData> output =
            resample::ResampleFrame(input, parameters, bounds);
        if (output == nullptr || output->GetNumberOfPoints() != 125) {
            throw std::runtime_error("Spectral resampler returned an invalid image.");
        }
        vtkDataArray* velocity = output->GetPointData()->GetArray("Velocity");
        if (velocity == nullptr || velocity->GetNumberOfComponents() != 3) {
            throw std::runtime_error("Spectral resampler dropped the velocity array.");
        }
        const vtkIdType pointId = 1 + 5 * 1 + 25 * 1;
        const double x = 0.5;
        const double y = 1.5;
        const double z = 2.5;
        if (!NearlyEqual(velocity->GetComponent(pointId, 0),
                         x * x + 2.0 * y + 3.0 * z) ||
            !NearlyEqual(velocity->GetComponent(pointId, 1), x * y) ||
            !NearlyEqual(velocity->GetComponent(pointId, 2), z * z * z)) {
            throw std::runtime_error(
                "Resampler did not preserve the represented GLL polynomial.");
        }
        if (output->GetFieldData()->GetArray(io::kNek5000GllDimensionsField) !=
            nullptr) {
            throw std::runtime_error(
                "Regular-grid output retained invalid spectral-element metadata.");
        }
        std::cout << "Nek5000 spectral resampling test passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Nek5000 spectral resampling test failed: " << error.what()
                  << '\n';
        return 1;
    }
}
