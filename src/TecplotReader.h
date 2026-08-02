#pragma once

#include <filesystem>

#include <vtkSmartPointer.h>

class vtkDataObject;

namespace streamcenterplus::io::detail {

// Classic Tecplot binary files begin with #!TDV112 or #!TDV191.  ASCII
// Tecplot .dat files (and ASCII data misnamed .plt) are handled by VTK's
// vtkTecplotReader instead.
bool IsClassicTecplotBinaryFile(const std::filesystem::path& path);

// Reads the subset of the documented classic binary PLT format that carries
// ORDERED, fixed finite-element, and face-based finite-element flow meshes.
// The caller owns composite flattening and source-block provenance so it
// behaves the same as the VTK readers.
vtkSmartPointer<vtkDataObject> ReadClassicTecplotBinary(
    const std::filesystem::path& path);

}  // namespace streamcenterplus::io::detail
