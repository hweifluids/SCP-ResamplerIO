#include "streamcenterplus/io/DataSetReader.h"
#include "TecplotReader.h"

#include <streamcenterplus/VtkHdfIO.h>

#include <vtkAppendFilter.h>
#include <vtkBox.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkCellSizeFilter.h>
#include <vtkCellType.h>
#include <vtkCallbackCommand.h>
#include <vtkClipDataSet.h>
#include <vtkContourFilter.h>
#include <vtkCutter.h>
#include <vtkDataSetAttributes.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataArray.h>
#include <vtkDataArrayCollection.h>
#include <vtkDataArraySelection.h>
#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkGenericCell.h>
#include <vtkGenericEnSightReader.h>
#include <vtkInformation.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkNonMergingPointLocator.h>
#include <vtkOpenFOAMReader.h>
#include <vtkPOpenFOAMReader.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkRemoveUnusedPoints.h>
#include <vtkSmartPointer.h>
#include <vtkSMPTools.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkStringArray.h>
#include <vtkStaticCellLocator.h>
#include <vtkStaticPointLocator.h>
#include <vtkTecplotReader.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLGenericDataObjectReader.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace streamcenterplus::io {

namespace {

constexpr const char* kSourceBlockIdArray = "SCP_SourceBlockId";
constexpr const char* kTimeValueArray = "TimeValue";

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string LowerExtension(const std::filesystem::path& path) {
    return ToLower(streamcenterplus::vtkhdf::Extension(path));
}

std::string TrimAscii(const std::string& value) {
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string StripComment(const std::string& line) {
    const std::size_t hash = line.find('#');
    return TrimAscii(hash == std::string::npos ? line : line.substr(0, hash));
}

bool IsFiniteTime(double value) {
    return std::isfinite(value);
}

bool IsRegularFileWithExtension(const std::filesystem::path& path,
                                const std::string& extension) {
    return std::filesystem::is_regular_file(path) && LowerExtension(path) == extension;
}

bool IsOrdinaryVtkExtension(const std::string& extension) {
    return extension == ".vtk" || extension == ".vti" || extension == ".vtr" ||
           extension == ".vts" || extension == ".vtu" || extension == ".vtp" ||
           extension == ".vtm" || extension == ".vtkhdf";
}

bool IsTecplotExtension(const std::string& extension) {
    return extension == ".dat" || extension == ".plt";
}

bool IsDirectDataExtension(const std::string& extension) {
    return IsOrdinaryVtkExtension(extension) || IsTecplotExtension(extension);
}

constexpr const char* kCanonicalVelocityArray = kMultiZoneCanonicalVelocityArray;

enum class MultiZoneOverlapPolicy {
    CutBackground,
    Blend50_50,
};

enum class MultiZoneMotion {
    Static,
    Rigid,
};

enum class MultiZoneVelocityFrame {
    World,
    Zone,
    Relative,
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct MultiZoneZone {
    std::filesystem::path source;
    MultiZoneMotion motion = MultiZoneMotion::Static;
    Vec3 translation;
    Vec3 rotationAxis{0.0, 0.0, 1.0};
    Vec3 rotationOrigin;
    Vec3 linearVelocity;
    double angularVelocity = 0.0;
    double phase = 0.0;
    MultiZoneVelocityFrame velocityFrame = MultiZoneVelocityFrame::World;
    std::string velocityArray;
    std::array<std::string, 3> velocityComponents;
};

struct MultiZoneManifest {
    std::filesystem::path path;
    std::string meshKind;
    MultiZoneOverlapPolicy overlapPolicy = MultiZoneOverlapPolicy::CutBackground;
    int backgroundZone = -1;
    std::array<MultiZoneZone, 2> zones;
};

bool IsZero(const Vec3& value) {
    return value.x == 0.0 && value.y == 0.0 && value.z == 0.0;
}

double Dot(const Vec3& lhs, const Vec3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vec3 Cross(const Vec3& lhs, const Vec3& rhs) {
    return {lhs.y * rhs.z - lhs.z * rhs.y,
            lhs.z * rhs.x - lhs.x * rhs.z,
            lhs.x * rhs.y - lhs.y * rhs.x};
}

Vec3 Add(const Vec3& lhs, const Vec3& rhs) {
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Vec3 Subtract(const Vec3& lhs, const Vec3& rhs) {
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Vec3 Scale(const Vec3& value, double factor) {
    return {value.x * factor, value.y * factor, value.z * factor};
}

double Norm(const Vec3& value) {
    return std::sqrt(Dot(value, value));
}

Vec3 ParseVec3(const std::string& value, const std::string& key) {
    std::string normalized = value;
    std::replace(normalized.begin(), normalized.end(), ',', ' ');
    std::replace(normalized.begin(), normalized.end(), ';', ' ');
    std::istringstream input(normalized);
    std::array<double, 3> values{};
    for (double& component : values) {
        if (!(input >> component) || !std::isfinite(component)) {
            throw std::runtime_error("Manifest key '" + key +
                                     "' must contain three finite numbers.");
        }
    }
    std::string trailing;
    if (input >> trailing) {
        throw std::runtime_error("Manifest key '" + key +
                                 "' must contain exactly three numbers.");
    }
    return {values[0], values[1], values[2]};
}

std::array<std::string, 3> ParseNameTriplet(const std::string& value,
                                             const std::string& key) {
    std::string normalized = value;
    std::replace(normalized.begin(), normalized.end(), ',', ' ');
    std::replace(normalized.begin(), normalized.end(), ';', ' ');
    std::istringstream input(normalized);
    std::array<std::string, 3> names{};
    for (std::string& name : names) {
        if (!(input >> name) || name.empty()) {
            throw std::runtime_error("Manifest key '" + key +
                                     "' must contain exactly three array names.");
        }
    }
    std::string trailing;
    if (input >> trailing) {
        throw std::runtime_error("Manifest key '" + key +
                                 "' must contain exactly three array names.");
    }
    return names;
}

double ParseFiniteDouble(const std::string& value, const std::string& key) {
    std::size_t used = 0;
    double parsed = 0.0;
    try {
        parsed = std::stod(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("Manifest key '" + key +
                                 "' must be a finite number.");
    }
    if (used != value.size() || !std::isfinite(parsed)) {
        throw std::runtime_error("Manifest key '" + key +
                                 "' must be a finite number.");
    }
    return parsed;
}

int ParsePositiveInt(const std::string& value, const std::string& key) {
    std::size_t used = 0;
    int parsed = 0;
    try {
        parsed = std::stoi(value, &used);
    } catch (const std::exception&) {
        throw std::runtime_error("Manifest key '" + key + "' must be an integer.");
    }
    if (used != value.size() || parsed <= 0) {
        throw std::runtime_error("Manifest key '" + key +
                                 "' must be a positive integer.");
    }
    return parsed;
}

const std::string* FindManifestValue(
    const std::map<std::string, std::string>& values,
    const std::string& key) {
    const auto found = values.find(key);
    return found == values.end() ? nullptr : &found->second;
}

std::string RequiredManifestValue(const std::map<std::string, std::string>& values,
                                  const std::string& key,
                                  const std::filesystem::path& manifestPath) {
    const std::string* value = FindManifestValue(values, key);
    if (value == nullptr || value->empty()) {
        throw std::runtime_error("Manifest '" + manifestPath.string() +
                                 "' is missing required key '" + key + "'.");
    }
    return *value;
}

std::string OptionalManifestValue(const std::map<std::string, std::string>& values,
                                  const std::string& key,
                                  const std::string& fallback) {
    const std::string* value = FindManifestValue(values, key);
    return value == nullptr ? fallback : *value;
}

std::string ManifestMeshKind(const std::map<std::string, std::string>& values,
                             const std::filesystem::path& manifestPath) {
    for (const char* key : {"mesh_kind", "grid_kind", "layout"}) {
        if (const std::string* value = FindManifestValue(values, key)) {
            const std::string kind = ToLower(*value);
            if (kind == "two_zone_static" || kind == "multizone") {
                return "static";
            }
            if (kind == "two_zone_nested") {
                return "nested";
            }
            if (kind == "two_zone_multigrid" || kind == "multi_grid") {
                return "multigrid";
            }
            if (kind == "two_zone_sliding") {
                return "sliding";
            }
            return kind;
        }
    }
    throw std::runtime_error("Manifest '" + manifestPath.string() +
                             "' is missing required key 'mesh_kind'.");
}

bool IsKnownMultiZoneManifestKey(const std::string& key) {
    if (key == "streamcenterplus_meshset" || key == "zone_count" ||
        key == "mesh_kind" || key == "grid_kind" || key == "layout" ||
        key == "overlap_policy" || key == "background_zone") {
        return true;
    }
    for (int index = 0; index < 2; ++index) {
        const std::string prefix = "zone" + std::to_string(index) + "_";
        for (const char* suffix : {"source", "motion", "translation", "rotation_axis",
                                   "rotation_origin", "linear_velocity", "angular_velocity",
                                   "phase", "velocity_frame", "velocity_array",
                                   "velocity_components"}) {
            if (key == prefix + suffix) {
                return true;
            }
        }
    }
    return false;
}

void ValidateMultiZoneManifestKeys(const std::map<std::string, std::string>& values,
                                   const std::filesystem::path& manifestPath) {
    int meshKindKeys = 0;
    for (const auto& [key, value] : values) {
        static_cast<void>(value);
        if (!IsKnownMultiZoneManifestKey(key)) {
            throw std::runtime_error("Unknown mesh-set manifest key '" + key +
                                     "' in '" + manifestPath.string() + "'.");
        }
        if (key == "mesh_kind" || key == "grid_kind" || key == "layout") {
            ++meshKindKeys;
        }
    }
    if (meshKindKeys != 1) {
        throw std::runtime_error("Manifest '" + manifestPath.string() +
                                 "' must declare exactly one of mesh_kind, grid_kind or layout.");
    }
}

MultiZoneManifest ParseMultiZoneManifest(const std::filesystem::path& sourcePath) {
    const std::filesystem::path realPath =
        streamcenterplus::vtkhdf::RealPath(sourcePath);
    if (!std::filesystem::is_regular_file(realPath)) {
        throw std::runtime_error("Mesh-set manifest does not exist: " +
                                 sourcePath.string());
    }
    if (LowerExtension(sourcePath) != ".scpm") {
        throw std::runtime_error("Expected a .scpm mesh-set manifest: " +
                                 sourcePath.string());
    }

    std::ifstream input(realPath);
    if (!input) {
        throw std::runtime_error("Cannot open mesh-set manifest: " +
                                 sourcePath.string());
    }
    std::map<std::string, std::string> values;
    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        line = StripComment(line);
        if (line.empty()) {
            continue;
        }
        const std::size_t separator = line.find('=');
        if (separator == std::string::npos) {
            throw std::runtime_error("Malformed mesh-set manifest line " +
                                     std::to_string(lineNumber) + " in '" +
                                     sourcePath.string() + "': expected key=value.");
        }
        const std::string key = ToLower(TrimAscii(line.substr(0, separator)));
        const std::string value = TrimAscii(line.substr(separator + 1));
        if (key.empty() || value.empty()) {
            throw std::runtime_error("Malformed mesh-set manifest line " +
                                     std::to_string(lineNumber) + " in '" +
                                     sourcePath.string() + ".");
        }
        if (!values.emplace(key, value).second) {
            throw std::runtime_error("Duplicate mesh-set manifest key '" + key +
                                     "' in '" + sourcePath.string() + "'.");
        }
    }

    ValidateMultiZoneManifestKeys(values, sourcePath);

    if (RequiredManifestValue(values, "streamcenterplus_meshset", sourcePath) != "1") {
        throw std::runtime_error("Manifest '" + sourcePath.string() +
                                 "' must declare streamcenterplus_meshset=1.");
    }
    if (ParsePositiveInt(RequiredManifestValue(values, "zone_count", sourcePath),
                         "zone_count") != 2) {
        throw std::runtime_error("Manifest '" + sourcePath.string() +
                                 "' must declare exactly zone_count=2.");
    }

    MultiZoneManifest manifest;
    manifest.path = std::filesystem::absolute(realPath).lexically_normal();
    manifest.meshKind = ManifestMeshKind(values, sourcePath);
    if (manifest.meshKind != "static" && manifest.meshKind != "nested" &&
        manifest.meshKind != "multigrid" && manifest.meshKind != "sliding") {
        throw std::runtime_error("Manifest '" + sourcePath.string() +
                                 "' has unsupported mesh_kind '" + manifest.meshKind +
                                 "'.");
    }

    const std::string overlap = ToLower(OptionalManifestValue(
        values, "overlap_policy", "cut_background"));
    if (overlap == "cut_background") {
        manifest.overlapPolicy = MultiZoneOverlapPolicy::CutBackground;
    } else if (overlap == "blend_50_50") {
        manifest.overlapPolicy = MultiZoneOverlapPolicy::Blend50_50;
    } else {
        throw std::runtime_error("Manifest '" + sourcePath.string() +
                                 "' has unsupported overlap_policy '" + overlap + "'.");
    }

    for (int index = 0; index < 2; ++index) {
        const std::string prefix = "zone" + std::to_string(index) + "_";
        MultiZoneZone& zone = manifest.zones[static_cast<std::size_t>(index)];
        const std::filesystem::path configuredSource = RequiredManifestValue(
            values, prefix + "source", sourcePath);
        zone.source = configuredSource.is_absolute()
                          ? configuredSource
                          : (manifest.path.parent_path() / configuredSource);
        zone.source = zone.source.lexically_normal();
        if (LowerExtension(zone.source) == ".scpm") {
            throw std::runtime_error("Nested .scpm sources are not supported by '" +
                                     sourcePath.string() + "'.");
        }

        const std::string motion = ToLower(OptionalManifestValue(
            values, prefix + "motion", "static"));
        if (motion == "static") {
            zone.motion = MultiZoneMotion::Static;
        } else if (motion == "rigid") {
            zone.motion = MultiZoneMotion::Rigid;
        } else {
            throw std::runtime_error("Manifest key '" + prefix +
                                     "motion' must be static or rigid.");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "translation")) {
            zone.translation = ParseVec3(*value, prefix + "translation");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "rotation_axis")) {
            zone.rotationAxis = ParseVec3(*value, prefix + "rotation_axis");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "rotation_origin")) {
            zone.rotationOrigin = ParseVec3(*value, prefix + "rotation_origin");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "linear_velocity")) {
            zone.linearVelocity = ParseVec3(*value, prefix + "linear_velocity");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "angular_velocity")) {
            zone.angularVelocity = ParseFiniteDouble(*value, prefix + "angular_velocity");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "phase")) {
            zone.phase = ParseFiniteDouble(*value, prefix + "phase");
        }
        const double axisNorm = Norm(zone.rotationAxis);
        if (axisNorm == 0.0) {
            throw std::runtime_error("Manifest key '" + prefix +
                                     "rotation_axis' must be non-zero.");
        }
        zone.rotationAxis = Scale(zone.rotationAxis, 1.0 / axisNorm);
        if (zone.motion == MultiZoneMotion::Static &&
            (zone.angularVelocity != 0.0 || !IsZero(zone.linearVelocity))) {
            throw std::runtime_error("Static zone '" + std::to_string(index) +
                                     "' cannot declare angular_velocity or linear_velocity.");
        }

        const std::string velocityFrame = ToLower(OptionalManifestValue(
            values, prefix + "velocity_frame", "world"));
        if (velocityFrame == "world") {
            zone.velocityFrame = MultiZoneVelocityFrame::World;
        } else if (velocityFrame == "zone") {
            zone.velocityFrame = MultiZoneVelocityFrame::Zone;
        } else if (velocityFrame == "relative") {
            zone.velocityFrame = MultiZoneVelocityFrame::Relative;
        } else {
            throw std::runtime_error("Manifest key '" + prefix +
                                     "velocity_frame' must be world, zone or relative.");
        }
        if (const std::string* value = FindManifestValue(values, prefix + "velocity_array")) {
            zone.velocityArray = *value;
        }
        if (const std::string* value = FindManifestValue(values, prefix + "velocity_components")) {
            if (!zone.velocityArray.empty()) {
                throw std::runtime_error("Manifest zone " + std::to_string(index) +
                                         " cannot declare both velocity_array and velocity_components.");
            }
            zone.velocityComponents = ParseNameTriplet(*value, prefix + "velocity_components");
        }
    }

    if (const std::string* value = FindManifestValue(values, "background_zone")) {
        std::size_t used = 0;
        int parsed = -1;
        try {
            parsed = std::stoi(*value, &used);
        } catch (const std::exception&) {
            throw std::runtime_error("Manifest key 'background_zone' must be 0 or 1.");
        }
        if (used != value->size() || (parsed != 0 && parsed != 1)) {
            throw std::runtime_error("Manifest key 'background_zone' must be 0 or 1.");
        }
        manifest.backgroundZone = parsed;
    }
    if (manifest.meshKind == "sliding" &&
        manifest.zones[0].motion != MultiZoneMotion::Rigid &&
        manifest.zones[1].motion != MultiZoneMotion::Rigid) {
        throw std::runtime_error("A sliding mesh-set manifest must declare at least one rigid zone.");
    }
    const bool hasDynamicRigidZone = std::any_of(
        manifest.zones.begin(), manifest.zones.end(), [](const MultiZoneZone& zone) {
            return zone.motion == MultiZoneMotion::Rigid &&
                   (zone.angularVelocity != 0.0 || !IsZero(zone.linearVelocity));
        });
    if (hasDynamicRigidZone && manifest.meshKind != "sliding") {
        throw std::runtime_error("A manifest with moving rigid zones must use mesh_kind=sliding.");
    }

    return manifest;
}

bool NearlyEqualTime(double lhs, double rhs) {
    return std::abs(lhs - rhs) <=
           1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

const DataSetFrame& SelectZoneFrame(const std::vector<DataSetFrame>& frames,
                                    int manifestFrameIndex,
                                    const std::filesystem::path& manifestPath) {
    if (frames.empty()) {
        throw std::runtime_error("Mesh-set source has no frames in '" +
                                 manifestPath.string() + "'.");
    }
    if (frames.size() == 1) {
        return frames.front();
    }
    if (manifestFrameIndex < 0 ||
        static_cast<std::size_t>(manifestFrameIndex) >= frames.size()) {
        throw std::runtime_error("Mesh-set frame index is out of range in '" +
                                 manifestPath.string() + "'.");
    }
    return frames[static_cast<std::size_t>(manifestFrameIndex)];
}

struct MeshReferenceMeasure {
    int dimension = -1;
    double value = 0.0;
};

std::optional<MeshReferenceMeasure> ComputeMeshReferenceMeasure(vtkDataSet* dataSet);

int ResolveMultiZoneBackground(const MultiZoneManifest& manifest,
                               const std::vector<DataSetFrame>& first,
                               const std::vector<DataSetFrame>& second) {
    int backgroundZone = manifest.backgroundZone;
    if (backgroundZone < 0) {
        // The choice is made from the two reference frames once, during
        // manifest enumeration.  It must never be re-evaluated from a moving
        // foreground frame because modal solvers require invariant DOF.
        ReadOptions geometryOnly;
        geometryOnly.readAllArrays = false;
        geometryOnly.attachSourceBlockIds = false;
        const bool nativeNekSource =
            first.front().kind == DataSetFrameKind::Nek5000TimeStep ||
            second.front().kind == DataSetFrameKind::Nek5000TimeStep;
        // Native Nek5000 reference measures come from the canonical spectral
        // quadrature arrays, never from invented GLL-to-GLL VTK cells.
        geometryOnly.attachNek5000SpectralMetadata = nativeNekSource;
        geometryOnly.attachNek5000QuadratureArrays = nativeNekSource;
        geometryOnly.failOnVtkXmlReadErrors = true;
        vtkSmartPointer<vtkDataSet> firstReference =
            ReadDataSet(first.front(), geometryOnly);
        vtkSmartPointer<vtkDataSet> secondReference =
            ReadDataSet(second.front(), geometryOnly);
        const std::optional<MeshReferenceMeasure> firstMeasure =
            ComputeMeshReferenceMeasure(firstReference);
        const std::optional<MeshReferenceMeasure> secondMeasure =
            ComputeMeshReferenceMeasure(secondReference);
        if (!firstMeasure.has_value() || !secondMeasure.has_value() ||
            firstMeasure->dimension != secondMeasure->dimension) {
            throw std::runtime_error(
                "Cannot determine comparable physical reference measures for both "
                "mesh-set zones in '" + manifest.path.string() +
                "'. Automatic background selection requires non-degenerate cells of "
                "the same physical dimension; declare background_zone explicitly "
                "after verifying the fixed background geometry.");
        }
        backgroundZone = firstMeasure->value >= secondMeasure->value ? 0 : 1;
    }
    const MultiZoneZone& zone = manifest.zones[static_cast<std::size_t>(backgroundZone)];
    if (zone.motion == MultiZoneMotion::Rigid &&
        (zone.angularVelocity != 0.0 || !IsZero(zone.linearVelocity))) {
        throw std::runtime_error("The selected background zone must be static so the "
                                 "canonical world reference geometry remains fixed: " +
                                 manifest.path.string());
    }
    return backgroundZone;
}

std::vector<DataSetFrame> EnumerateMultiZoneManifestFrames(
    const std::filesystem::path& sourcePath,
    double fallbackPhysicalDt) {
    const MultiZoneManifest manifest = ParseMultiZoneManifest(sourcePath);
    // A zone source may be a PVD collection.  It must still resolve to an
    // ordinary dataset sequence, not another mesh-set manifest: recursively
    // composing a two-zone result would exceed the public two-zone contract
    // and apply a second reference-frame transform to already canonical data.
    // Track the active manifests as well so a PVD-mediated self-reference
    // fails deterministically instead of recursing through enumeration.
    static thread_local std::vector<std::filesystem::path> activeManifests;
    const auto active = std::find(activeManifests.begin(), activeManifests.end(), manifest.path);
    if (active != activeManifests.end()) {
        throw std::runtime_error("Recursive mesh-set manifest source is not supported: " +
                                 manifest.path.string());
    }
    activeManifests.push_back(manifest.path);
    struct ActiveManifestScope {
        std::vector<std::filesystem::path>* paths = nullptr;
        ~ActiveManifestScope() {
            if (paths != nullptr) {
                paths->pop_back();
            }
        }
    } activeManifestScope{&activeManifests};

    const std::vector<DataSetFrame> first =
        EnumerateDataSetFrames(manifest.zones[0].source, fallbackPhysicalDt);
    const std::vector<DataSetFrame> second =
        EnumerateDataSetFrames(manifest.zones[1].source, fallbackPhysicalDt);
    const auto rejectNestedManifestFrames =
        [&](const std::vector<DataSetFrame>& zoneFrames, int zoneIndex) {
            if (std::any_of(zoneFrames.begin(), zoneFrames.end(),
                            [](const DataSetFrame& frame) {
                                return frame.kind ==
                                       DataSetFrameKind::MultiZoneManifestTimeStep;
                            })) {
                throw std::runtime_error(
                    "Nested .scpm sources, including PVD-indirected mesh-set manifests, "
                    "are not supported by '" + manifest.path.string() +
                    "' (zone" + std::to_string(zoneIndex) + ").");
            }
        };
    rejectNestedManifestFrames(first, 0);
    rejectNestedManifestFrames(second, 1);
    if (first.empty() || second.empty()) {
        throw std::runtime_error("Every mesh-set zone must provide at least one frame: " +
                                 sourcePath.string());
    }
    if (first.size() > 1 && second.size() > 1 && first.size() != second.size()) {
        throw std::runtime_error("Two dynamic mesh-set sources must have equal frame counts: " +
                                 sourcePath.string());
    }
    const int backgroundZone = ResolveMultiZoneBackground(manifest, first, second);
    const std::size_t count = std::max(first.size(), second.size());
    std::vector<DataSetFrame> result;
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        const DataSetFrame& firstFrame = SelectZoneFrame(
            first, static_cast<int>(index), sourcePath);
        const DataSetFrame& secondFrame = SelectZoneFrame(
            second, static_cast<int>(index), sourcePath);
        if (first.size() > 1 && second.size() > 1 &&
            !NearlyEqualTime(firstFrame.time, secondFrame.time)) {
            throw std::runtime_error("Two dynamic mesh-set sources must have matching frame times: " +
                                     sourcePath.string());
        }
        DataSetFrame frame;
        frame.path = manifest.path;
        frame.time = first.size() > 1 ? firstFrame.time : secondFrame.time;
        frame.kind = DataSetFrameKind::MultiZoneManifestTimeStep;
        frame.sourceIndex = static_cast<int>(index);
        frame.multiZoneBackgroundIndex = backgroundZone;
        result.push_back(std::move(frame));
    }
    return result;
}

vtkSmartPointer<vtkDataSet> CloneDataSet(vtkDataSet* input) {
    if (input == nullptr) {
        return nullptr;
    }
    vtkSmartPointer<vtkDataSet> output;
    output.TakeReference(vtkDataSet::SafeDownCast(input->NewInstance()));
    if (output == nullptr) {
        throw std::runtime_error("Cannot clone a VTK dataset instance.");
    }
    output->DeepCopy(input);
    return output;
}

void CollectDataSets(vtkDataObject* object,
                     std::vector<vtkSmartPointer<vtkDataSet>>* dataSets) {
    if (object == nullptr || dataSets == nullptr) {
        return;
    }
    if (vtkDataSet* dataSet = vtkDataSet::SafeDownCast(object)) {
        if (dataSet->GetNumberOfPoints() > 0) {
            dataSets->push_back(CloneDataSet(dataSet));
        }
        return;
    }
    vtkCompositeDataSet* composite = vtkCompositeDataSet::SafeDownCast(object);
    if (composite == nullptr) {
        return;
    }
    vtkSmartPointer<vtkCompositeDataIterator> iterator;
    iterator.TakeReference(composite->NewIterator());
    if (iterator == nullptr) {
        return;
    }
    iterator->SkipEmptyNodesOn();
    for (iterator->InitTraversal(); !iterator->IsDoneWithTraversal();
         iterator->GoToNextItem()) {
        CollectDataSets(iterator->GetCurrentDataObject(), dataSets);
    }
}

void AttachSourceBlockId(vtkDataSet* dataSet, int blockId) {
    if (dataSet == nullptr) {
        return;
    }
    vtkCellData* cellData = dataSet->GetCellData();
    if (cellData == nullptr) {
        return;
    }
    vtkDataArray* existing =
        cellData->GetArray(kSourceBlockIdArray);
    if (existing != nullptr &&
        existing->GetNumberOfComponents() == 1 &&
        existing->GetNumberOfTuples() ==
            dataSet->GetNumberOfCells()) {
        // A reader or upstream conversion may already carry stable,
        // application-level block provenance. Preserve it instead of
        // replacing it with the transient composite traversal index.
        return;
    }
    if (existing != nullptr) {
        cellData->RemoveArray(kSourceBlockIdArray);
    }
    vtkNew<vtkIntArray> cells;
    cells->SetName(kSourceBlockIdArray);
    cells->SetNumberOfTuples(dataSet->GetNumberOfCells());
    cells->FillValue(blockId);
    cellData->AddArray(cells);
}

vtkSmartPointer<vtkDataSet> CollapseDataObject(vtkDataObject* object,
                                               const ReadOptions& options,
                                               const std::filesystem::path& sourcePath) {
    std::vector<vtkSmartPointer<vtkDataSet>> dataSets;
    CollectDataSets(object, &dataSets);
    if (dataSets.empty()) {
        throw std::runtime_error("The source produced no non-empty VTK dataset: " +
                                 sourcePath.string());
    }
    if (dataSets.size() == 1) {
        if (options.attachSourceBlockIds) {
            AttachSourceBlockId(dataSets.front(), 0);
        }
        return dataSets.front();
    }

    vtkNew<vtkAppendFilter> append;
    if (options.preserveCompositeBlockBoundaries) {
        append->MergePointsOff();
    } else {
        append->MergePointsOn();
    }
    for (std::size_t i = 0; i < dataSets.size(); ++i) {
        if (options.attachSourceBlockIds) {
            AttachSourceBlockId(dataSets[i], static_cast<int>(i));
        }
        append->AddInputData(dataSets[i]);
    }
    append->Update();
    vtkSmartPointer<vtkDataSet> output = CloneDataSet(append->GetOutput());
    if (output == nullptr || output->GetNumberOfPoints() <= 0) {
        throw std::runtime_error("Failed to assemble the volume blocks from: " +
                                 sourcePath.string());
    }
    return output;
}

Vec3 RotateZoneVector(const MultiZoneZone& zone, const Vec3& value, double time) {
    const double angle = zone.phase +
                         (zone.motion == MultiZoneMotion::Rigid
                              ? zone.angularVelocity * time
                              : 0.0);
    if (angle == 0.0) {
        return value;
    }
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    return Add(Add(Scale(value, cosine), Scale(Cross(zone.rotationAxis, value), sine)),
               Scale(zone.rotationAxis, Dot(zone.rotationAxis, value) * (1.0 - cosine)));
}

Vec3 ZoneTranslationAt(const MultiZoneZone& zone, double time) {
    if (zone.motion != MultiZoneMotion::Rigid) {
        return zone.translation;
    }
    return Add(zone.translation, Scale(zone.linearVelocity, time));
}

Vec3 TransformZonePointToWorld(const MultiZoneZone& zone,
                               const Vec3& localPoint,
                               double time) {
    const Vec3 translation = ZoneTranslationAt(zone, time);
    return Add(Add(translation, zone.rotationOrigin),
               RotateZoneVector(zone, Subtract(localPoint, zone.rotationOrigin), time));
}

Vec3 MeshVelocityAtWorldPoint(const MultiZoneZone& zone,
                              const Vec3& worldPoint,
                              double time) {
    if (zone.motion != MultiZoneMotion::Rigid) {
        return {};
    }
    const Vec3 centre = Add(ZoneTranslationAt(zone, time), zone.rotationOrigin);
    return Add(zone.linearVelocity,
               Scale(Cross(zone.rotationAxis, Subtract(worldPoint, centre)),
                     zone.angularVelocity));
}

vtkDataArray* FindDataArray(vtkDataSetAttributes* attributes, const std::string& name) {
    if (attributes == nullptr || name.empty()) {
        return nullptr;
    }
    if (vtkDataArray* exact = attributes->GetArray(name.c_str())) {
        return exact;
    }
    const std::string expected = ToLower(name);
    for (int index = 0; index < attributes->GetNumberOfArrays(); ++index) {
        vtkDataArray* candidate = attributes->GetArray(index);
        if (candidate != nullptr && candidate->GetName() != nullptr &&
            ToLower(candidate->GetName()) == expected) {
            return candidate;
        }
    }
    return nullptr;
}

vtkDataArray* FindExactDataArray(vtkDataSetAttributes* attributes, const char* name) {
    return attributes == nullptr ? nullptr : attributes->GetArray(name);
}

struct VelocitySelection {
    vtkDataArray* vector = nullptr;
    std::array<vtkDataArray*, 3> components = {nullptr, nullptr, nullptr};
    std::vector<std::string> explicitlyNamedArrays;
};

VelocitySelection ResolveVelocitySelection(vtkDataSetAttributes* attributes,
                                           const MultiZoneZone& zone,
                                           int zoneIndex) {
    if (attributes == nullptr) {
        throw std::runtime_error("Mesh-set zone " + std::to_string(zoneIndex) +
                                 " has no point data for Velocity.");
    }
    VelocitySelection result;
    if (!zone.velocityArray.empty()) {
        result.vector = FindDataArray(attributes, zone.velocityArray);
        if (result.vector == nullptr || result.vector->GetNumberOfComponents() != 3) {
            throw std::runtime_error("Manifest velocity_array for zone " +
                                     std::to_string(zoneIndex) +
                                     " must name a three-component point array.");
        }
        result.explicitlyNamedArrays.push_back(result.vector->GetName());
        return result;
    }

    if (!zone.velocityComponents[0].empty()) {
        for (int component = 0; component < 3; ++component) {
            result.components[static_cast<std::size_t>(component)] = FindDataArray(
                attributes, zone.velocityComponents[static_cast<std::size_t>(component)]);
            vtkDataArray* array = result.components[static_cast<std::size_t>(component)];
            if (array == nullptr || array->GetNumberOfComponents() != 1) {
                throw std::runtime_error("Manifest velocity_components for zone " +
                                         std::to_string(zoneIndex) +
                                         " must name three scalar point arrays.");
            }
            result.explicitlyNamedArrays.push_back(array->GetName());
        }
        return result;
    }

    for (const char* candidate : {"Velocity", "U", "V"}) {
        vtkDataArray* array = FindDataArray(attributes, candidate);
        if (array != nullptr && array->GetNumberOfComponents() == 3) {
            result.vector = array;
            result.explicitlyNamedArrays.push_back(array->GetName());
            return result;
        }
    }

    const auto selectComponents = [&](const std::array<const char*, 3>& names) {
        std::array<vtkDataArray*, 3> arrays = {
            FindExactDataArray(attributes, names[0]),
            FindExactDataArray(attributes, names[1]),
            FindExactDataArray(attributes, names[2]),
        };
        for (vtkDataArray* array : arrays) {
            if (array == nullptr || array->GetNumberOfComponents() != 1) {
                return false;
            }
        }
        result.components = arrays;
        for (vtkDataArray* array : arrays) {
            result.explicitlyNamedArrays.push_back(array->GetName());
        }
        return true;
    };
    if (selectComponents({"u", "v", "w"}) ||
        selectComponents({"U", "V", "W"})) {
        return result;
    }

    throw std::runtime_error(
        "Mesh-set zone " + std::to_string(zoneIndex) +
        " must expose a three-component Velocity/U/V point array, or scalar u/v/w or U/V/W arrays.");
}

Vec3 ReadVelocityTuple(const VelocitySelection& selection, vtkIdType pointId) {
    if (selection.vector != nullptr) {
        return {selection.vector->GetComponent(pointId, 0),
                selection.vector->GetComponent(pointId, 1),
                selection.vector->GetComponent(pointId, 2)};
    }
    return {selection.components[0]->GetComponent(pointId, 0),
            selection.components[1]->GetComponent(pointId, 0),
            selection.components[2]->GetComponent(pointId, 0)};
}

void RequireCompleteVelocityTuples(const VelocitySelection& selection,
                                   vtkIdType expectedTupleCount,
                                   int zoneIndex,
                                   const char* association) {
    if (expectedTupleCount < 0) {
        throw std::runtime_error("Invalid expected tuple count for mesh-set zone " +
                                 std::to_string(zoneIndex) + ".");
    }
    const auto requireArray = [&](vtkDataArray* array) {
        if (array == nullptr || array->GetName() == nullptr ||
            array->GetNumberOfTuples() != expectedTupleCount) {
            const std::string name = array != nullptr && array->GetName() != nullptr
                                         ? array->GetName()
                                         : std::string("<unnamed>");
            throw std::runtime_error(
                std::string(association) + " velocity array '" + name +
                "' in mesh-set zone " + std::to_string(zoneIndex) +
                " must provide exactly " + std::to_string(expectedTupleCount) +
                " tuples.");
        }
    };
    if (selection.vector != nullptr) {
        requireArray(selection.vector);
        return;
    }
    for (vtkDataArray* component : selection.components) {
        requireArray(component);
    }
}

bool HasExplicitVelocitySelection(const MultiZoneZone& zone) {
    return !zone.velocityArray.empty() || !zone.velocityComponents[0].empty();
}

bool IsStandardVelocityAlias(const char* name) {
    if (name == nullptr) {
        return false;
    }
    const std::string normalized = ToLower(name);
    return normalized == "velocity" || normalized == "u" || normalized == "v" ||
           normalized == "w";
}

void RemoveRawVelocityArrays(vtkDataSetAttributes* attributes,
                             const VelocitySelection& selection,
                             bool removeStandardAliases = true) {
    if (attributes == nullptr) {
        return;
    }
    std::vector<std::string> names;
    names.reserve(static_cast<std::size_t>(attributes->GetNumberOfArrays()));
    for (int index = 0; index < attributes->GetNumberOfArrays(); ++index) {
        vtkDataArray* array = attributes->GetArray(index);
        if (array == nullptr || array->GetName() == nullptr) {
            continue;
        }
        bool remove = removeStandardAliases && IsStandardVelocityAlias(array->GetName());
        for (const std::string& explicitName : selection.explicitlyNamedArrays) {
            if (array->GetName() == explicitName) {
                remove = true;
                break;
            }
        }
        if (remove) {
            names.emplace_back(array->GetName());
        }
    }
    for (const std::string& name : names) {
        attributes->RemoveArray(name.c_str());
    }
}

void PromoteCellVelocityToPointData(vtkUnstructuredGrid* dataSet,
                                    const VelocitySelection& selection,
                                    int zoneIndex) {
    if (dataSet == nullptr || dataSet->GetCellData() == nullptr ||
        dataSet->GetPointData() == nullptr) {
        throw std::runtime_error("Cannot promote cell velocity for mesh-set zone " +
                                 std::to_string(zoneIndex) + ".");
    }
    std::vector<vtkDataArray*> sourceArrays;
    if (selection.vector != nullptr) {
        sourceArrays.push_back(selection.vector);
    } else {
        sourceArrays.assign(selection.components.begin(), selection.components.end());
    }
    vtkNew<vtkIdList> pointIds;
    const vtkIdType pointCount = dataSet->GetNumberOfPoints();
    for (vtkDataArray* source : sourceArrays) {
        if (source == nullptr || source->GetName() == nullptr ||
            source->GetNumberOfTuples() != dataSet->GetNumberOfCells()) {
            const std::string name = source != nullptr && source->GetName() != nullptr
                                         ? source->GetName()
                                         : std::string("<unnamed>");
            throw std::runtime_error("Cell velocity array '" + name +
                                     "' in mesh-set zone " +
                                     std::to_string(zoneIndex) +
                                     " must provide exactly " +
                                     std::to_string(dataSet->GetNumberOfCells()) +
                                     " tuples.");
        }
        const int componentCount = source->GetNumberOfComponents();
        std::vector<double> sums(static_cast<std::size_t>(pointCount) *
                                     static_cast<std::size_t>(componentCount),
                                 0.0);
        std::vector<vtkIdType> counts(static_cast<std::size_t>(pointCount), 0);
        for (vtkIdType cellId = 0; cellId < dataSet->GetNumberOfCells(); ++cellId) {
            dataSet->GetCellPoints(cellId, pointIds);
            for (vtkIdType localId = 0; localId < pointIds->GetNumberOfIds(); ++localId) {
                const vtkIdType pointId = pointIds->GetId(localId);
                ++counts[static_cast<std::size_t>(pointId)];
                for (int component = 0; component < componentCount; ++component) {
                    sums[static_cast<std::size_t>(pointId) *
                             static_cast<std::size_t>(componentCount) +
                         static_cast<std::size_t>(component)] +=
                        source->GetComponent(cellId, component);
                }
            }
        }
        vtkSmartPointer<vtkDoubleArray> promoted =
            vtkSmartPointer<vtkDoubleArray>::New();
        promoted->SetName(source->GetName());
        promoted->SetNumberOfComponents(componentCount);
        promoted->SetNumberOfTuples(pointCount);
        for (vtkIdType pointId = 0; pointId < pointCount; ++pointId) {
            const vtkIdType count = counts[static_cast<std::size_t>(pointId)];
            if (count <= 0) {
                throw std::runtime_error("Cell velocity cannot be projected to an orphan point "
                                         "in mesh-set zone " + std::to_string(zoneIndex) + ".");
            }
            for (int component = 0; component < componentCount; ++component) {
                promoted->SetComponent(
                    pointId, component,
                    sums[static_cast<std::size_t>(pointId) *
                             static_cast<std::size_t>(componentCount) +
                         static_cast<std::size_t>(component)] /
                        static_cast<double>(count));
            }
        }
        dataSet->GetPointData()->RemoveArray(source->GetName());
        dataSet->GetPointData()->AddArray(promoted);
    }
    // Keep a single world-frame velocity representation. Leaving the source
    // cell arrays would let Resample's cell-to-point preparation replace the
    // canonical point Velocity with stale zone-frame values.
    RemoveRawVelocityArrays(dataSet->GetCellData(), selection);
}

bool TryResolveVelocitySelection(vtkUnstructuredGrid* dataSet,
                                 const MultiZoneZone& zone,
                                 int zoneIndex,
                                 VelocitySelection* selection) {
    if (dataSet == nullptr || selection == nullptr) {
        throw std::runtime_error("Cannot resolve mesh-set velocity from a null zone.");
    }
    const bool explicitSelection = HasExplicitVelocitySelection(zone);
    try {
        *selection = ResolveVelocitySelection(dataSet->GetPointData(), zone, zoneIndex);
        return true;
    } catch (const std::runtime_error&) {
        vtkCellData* cellData = dataSet->GetCellData();
        if (cellData == nullptr || cellData->GetNumberOfArrays() == 0) {
            if (explicitSelection) {
                throw;
            }
            return false;
        }

        VelocitySelection cellSelection;
        try {
            cellSelection = ResolveVelocitySelection(cellData, zone, zoneIndex);
        } catch (const std::runtime_error&) {
            if (explicitSelection) {
                throw;
            }
            return false;
        }
        PromoteCellVelocityToPointData(dataSet, cellSelection, zoneIndex);
        *selection = ResolveVelocitySelection(dataSet->GetPointData(), zone, zoneIndex);
        return true;
    }
}

void ValidateSourceVelocityTupleCounts(vtkDataSet* source,
                                       const MultiZoneZone& zone,
                                       int zoneIndex) {
    if (source == nullptr) {
        throw std::runtime_error("Cannot validate velocity arrays for a null mesh-set zone.");
    }

    // Validate on the source before vtkAppendFilter materializes an
    // unstructured grid. The append operation may allocate a full point-array
    // length for a truncated input array, masking corrupt tail tuples behind
    // default-initialized values.
    VelocitySelection pointSelection;
    bool hasPointSelection = false;
    try {
        pointSelection = ResolveVelocitySelection(source->GetPointData(), zone, zoneIndex);
        hasPointSelection = true;
    } catch (const std::runtime_error&) {
    }
    if (hasPointSelection) {
        RequireCompleteVelocityTuples(pointSelection, source->GetNumberOfPoints(), zoneIndex,
                                      "Point");
        return;
    }

    VelocitySelection cellSelection;
    bool hasCellSelection = false;
    try {
        cellSelection = ResolveVelocitySelection(source->GetCellData(), zone, zoneIndex);
        hasCellSelection = true;
    } catch (const std::runtime_error&) {
    }
    if (hasCellSelection) {
        RequireCompleteVelocityTuples(cellSelection, source->GetNumberOfCells(), zoneIndex,
                                      "Cell");
    }
}

vtkSmartPointer<vtkUnstructuredGrid> MakeWorldZone(vtkDataSet* source,
                                                    const MultiZoneZone& zone,
                                                    int zoneIndex,
                                                    double time,
                                                    bool includeVelocity,
                                                    bool* hasWorldVelocity = nullptr) {
    if (source == nullptr || source->GetNumberOfPoints() <= 0 ||
        source->GetNumberOfCells() <= 0) {
        throw std::runtime_error("Mesh-set zone " + std::to_string(zoneIndex) +
                                 " has no volume cells or points.");
    }
    if (includeVelocity) {
        ValidateSourceVelocityTupleCounts(source, zone, zoneIndex);
    }
    vtkNew<vtkAppendFilter> append;
    append->MergePointsOff();
    append->AddInputData(source);
    append->Update();
    vtkSmartPointer<vtkUnstructuredGrid> result =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    result->DeepCopy(append->GetOutput());
    if (result->GetPoints() == nullptr || result->GetPointData() == nullptr) {
        throw std::runtime_error("Cannot materialize mesh-set zone " +
                                 std::to_string(zoneIndex) + " as an unstructured grid.");
    }

    vtkPointData* pointData = result->GetPointData();
    VelocitySelection velocity;
    vtkSmartPointer<vtkDoubleArray> worldVelocity;
    bool resolvedVelocity = false;
    if (includeVelocity) {
        resolvedVelocity = TryResolveVelocitySelection(result, zone, zoneIndex, &velocity);
        if (resolvedVelocity) {
            RequireCompleteVelocityTuples(velocity, result->GetNumberOfPoints(), zoneIndex,
                                          "Point");
            worldVelocity = vtkSmartPointer<vtkDoubleArray>::New();
            worldVelocity->SetName(kCanonicalVelocityArray);
            worldVelocity->SetNumberOfComponents(3);
            worldVelocity->SetNumberOfTuples(result->GetNumberOfPoints());
        }
    }

    double point[3] = {};
    for (vtkIdType pointId = 0; pointId < result->GetNumberOfPoints(); ++pointId) {
        result->GetPoint(pointId, point);
        const Vec3 localPoint{point[0], point[1], point[2]};
        const Vec3 worldPoint = TransformZonePointToWorld(zone, localPoint, time);
        if (resolvedVelocity) {
            const Vec3 inputVelocity = ReadVelocityTuple(velocity, pointId);
            Vec3 transformedVelocity = inputVelocity;
            if (zone.velocityFrame != MultiZoneVelocityFrame::World) {
                transformedVelocity = RotateZoneVector(zone, inputVelocity, time);
            }
            if (zone.velocityFrame == MultiZoneVelocityFrame::Relative) {
                transformedVelocity = Add(transformedVelocity,
                                          MeshVelocityAtWorldPoint(zone, worldPoint, time));
            }
            const double tuple[3] = {
                transformedVelocity.x, transformedVelocity.y, transformedVelocity.z};
            worldVelocity->SetTuple(pointId, tuple);
        }
        result->GetPoints()->SetPoint(pointId, worldPoint.x, worldPoint.y, worldPoint.z);
    }
    result->GetPoints()->Modified();
    result->Modified();
    if (resolvedVelocity) {
        RemoveRawVelocityArrays(pointData, velocity);
        RemoveRawVelocityArrays(result->GetCellData(), velocity);
        pointData->RemoveArray(kCanonicalVelocityArray);
        pointData->AddArray(worldVelocity);
        pointData->SetActiveVectors(kCanonicalVelocityArray);
    }
    if (hasWorldVelocity != nullptr) {
        *hasWorldVelocity = resolvedVelocity;
    }
    return result;
}

void RequireSameReferenceGeometry(vtkDataSet* reference,
                                  vtkDataSet* current,
                                  const std::filesystem::path& manifestPath,
                                  const char* subject) {
    if (reference == nullptr || current == nullptr ||
        reference->GetNumberOfPoints() != current->GetNumberOfPoints() ||
        reference->GetNumberOfCells() != current->GetNumberOfCells()) {
        throw std::runtime_error(std::string(subject) + " geometry changes across frames in '" +
                                 manifestPath.string() + "'.");
    }
    double referencePoint[3] = {};
    double currentPoint[3] = {};
    for (vtkIdType pointId = 0; pointId < reference->GetNumberOfPoints(); ++pointId) {
        reference->GetPoint(pointId, referencePoint);
        current->GetPoint(pointId, currentPoint);
        for (int axis = 0; axis < 3; ++axis) {
            const double tolerance = 1.0e-11 * std::max(
                {1.0, std::abs(referencePoint[axis]), std::abs(currentPoint[axis])});
            if (std::abs(referencePoint[axis] - currentPoint[axis]) > tolerance) {
                throw std::runtime_error(
                    std::string(subject) + " geometry changes across frames in '" +
                    manifestPath.string() + "'.");
            }
        }
    }
    vtkNew<vtkIdList> referencePointIds;
    vtkNew<vtkIdList> currentPointIds;
    for (vtkIdType cellId = 0; cellId < reference->GetNumberOfCells(); ++cellId) {
        if (reference->GetCellType(cellId) != current->GetCellType(cellId)) {
            throw std::runtime_error(std::string(subject) + " topology changes across frames in '" +
                                     manifestPath.string() + "'.");
        }
        reference->GetCellPoints(cellId, referencePointIds);
        current->GetCellPoints(cellId, currentPointIds);
        if (referencePointIds->GetNumberOfIds() != currentPointIds->GetNumberOfIds()) {
            throw std::runtime_error(std::string(subject) + " topology changes across frames in '" +
                                     manifestPath.string() + "'.");
        }
        for (vtkIdType pointIndex = 0;
             pointIndex < referencePointIds->GetNumberOfIds(); ++pointIndex) {
            if (referencePointIds->GetId(pointIndex) != currentPointIds->GetId(pointIndex)) {
                throw std::runtime_error(
                    std::string(subject) + " topology changes across frames in '" +
                    manifestPath.string() + "'.");
            }
        }
    }
}

std::optional<MeshReferenceMeasure> ComputeMeshReferenceMeasure(vtkDataSet* dataSet) {
    if (dataSet == nullptr) {
        return std::nullopt;
    }

    if (HasNek5000SpectralMetadata(dataSet) &&
        !IsNek5000SpectralDisplayProxy(dataSet)) {
        const Nek5000SpectralMetadata metadata =
            GetNek5000SpectralMetadata(dataSet);
        vtkDataArray* weights = dataSet->GetPointData() == nullptr
            ? nullptr
            : dataSet->GetPointData()->GetArray(kNek5000QuadratureWeightArray);
        if (weights == nullptr || weights->GetNumberOfComponents() != 1 ||
            weights->GetNumberOfTuples() != dataSet->GetNumberOfPoints()) {
            return std::nullopt;
        }
        double measure = 0.0;
        for (vtkIdType pointId = 0; pointId < weights->GetNumberOfTuples(); ++pointId) {
            const double weight = weights->GetComponent(pointId, 0);
            if (!std::isfinite(weight) || weight < 0.0) {
                return std::nullopt;
            }
            measure += weight;
        }
        if (!std::isfinite(measure) || measure <= 0.0) {
            return std::nullopt;
        }
        return MeshReferenceMeasure{metadata.is3D() ? 3 : 2, measure};
    }

    if (dataSet->GetNumberOfCells() <= 0) {
        return std::nullopt;
    }

    vtkDataArray* ghosts = dataSet->GetCellData() == nullptr
                                ? nullptr
                                : dataSet->GetCellData()->GetArray(
                                      vtkDataSetAttributes::GhostArrayName());
    const vtkIdType cellCount = dataSet->GetNumberOfCells();
    if (ghosts != nullptr && ghosts->GetNumberOfTuples() != cellCount) {
        return std::nullopt;
    }

    // A reference can legitimately be a surface or a line mesh.  Select the
    // highest physical dimension shared by the zones (3-D volume, 2-D area,
    // or 1-D length), rather than comparing AABBs or rejecting valid 2-D
    // mesh sets.  Mixed-dimensional references are rejected by the caller.
    int dimension = -1;
    for (vtkIdType cellId = 0; cellId < cellCount; ++cellId) {
        if (ghosts != nullptr && ghosts->GetComponent(cellId, 0) != 0.0) {
            continue;
        }
        vtkCell* cell = dataSet->GetCell(cellId);
        if (cell == nullptr) {
            return std::nullopt;
        }
        const int cellDimension = cell->GetCellDimension();
        if (cellDimension >= 1 && cellDimension <= 3) {
            dimension = std::max(dimension, cellDimension);
        }
    }
    if (dimension < 1) {
        return std::nullopt;
    }

    constexpr const char* kReferenceLengthArray = "__SCP_MeshReferenceLength";
    constexpr const char* kReferenceAreaArray = "__SCP_MeshReferenceArea";
    constexpr const char* kReferenceVolumeArray = "__SCP_MeshReferenceVolume";
    vtkNew<vtkCellSizeFilter> cellSizes;
    cellSizes->SetInputData(dataSet);
    cellSizes->SetComputeVertexCount(false);
    cellSizes->SetComputeLength(dimension == 1);
    cellSizes->SetComputeArea(dimension == 2);
    cellSizes->SetComputeVolume(dimension == 3);
    cellSizes->SetComputeSum(false);
    cellSizes->SetLengthArrayName(kReferenceLengthArray);
    cellSizes->SetAreaArrayName(kReferenceAreaArray);
    cellSizes->SetVolumeArrayName(kReferenceVolumeArray);
    cellSizes->Update();

    vtkDataSet* measured = vtkDataSet::SafeDownCast(cellSizes->GetOutputDataObject(0));
    const char* measureArrayName = dimension == 3
                                       ? kReferenceVolumeArray
                                       : (dimension == 2 ? kReferenceAreaArray
                                                         : kReferenceLengthArray);
    vtkDataArray* measures = measured == nullptr || measured->GetCellData() == nullptr
                                 ? nullptr
                                 : measured->GetCellData()->GetArray(measureArrayName);
    if (measures == nullptr || measures->GetNumberOfTuples() != cellCount) {
        return std::nullopt;
    }

    bool hasPhysicalCell = false;
    double totalMeasure = 0.0;
    for (vtkIdType cellId = 0; cellId < cellCount; ++cellId) {
        if (ghosts != nullptr && ghosts->GetComponent(cellId, 0) != 0.0) {
            continue;
        }
        vtkCell* cell = dataSet->GetCell(cellId);
        if (cell == nullptr) {
            return std::nullopt;
        }
        if (cell->GetCellDimension() != dimension) {
            continue;
        }
        const double measure = measures->GetComponent(cellId, 0);
        if (!std::isfinite(measure) || measure <= 0.0) {
            return std::nullopt;
        }
        hasPhysicalCell = true;
        totalMeasure += measure;
        if (!std::isfinite(totalMeasure)) {
            return std::nullopt;
        }
    }
    if (!hasPhysicalCell || totalMeasure <= 0.0) {
        return std::nullopt;
    }
    return MeshReferenceMeasure{dimension, totalMeasure};
}

bool FindContainingCell(vtkStaticCellLocator* locator,
                        vtkGenericCell* cell,
                        std::vector<double>* weights,
                        const double point[3]);

void RequireForegroundWithinBackground(vtkDataSet* background,
                                       vtkDataSet* foreground,
                                       const std::filesystem::path& manifestPath) {
    if (background == nullptr || foreground == nullptr) {
        throw std::runtime_error("Cannot validate mesh-set foreground containment in '" +
                                 manifestPath.string() + "'.");
    }
    double backgroundBounds[6] = {};
    double foregroundBounds[6] = {};
    background->GetBounds(backgroundBounds);
    foreground->GetBounds(foregroundBounds);
    for (int axis = 0; axis < 3; ++axis) {
        const double scale = std::max(
            {1.0, std::abs(backgroundBounds[2 * axis]),
             std::abs(backgroundBounds[2 * axis + 1])});
        const double tolerance = 1.0e-9 * scale;
        if (foregroundBounds[2 * axis] < backgroundBounds[2 * axis] - tolerance ||
            foregroundBounds[2 * axis + 1] > backgroundBounds[2 * axis + 1] + tolerance) {
            throw std::runtime_error(
                "The foreground zone extends beyond the fixed background reference mesh in '" +
                manifestPath.string() + "'. Use a larger background zone or a pre-resampled input.");
        }
    }

    vtkNew<vtkStaticCellLocator> backgroundLocator;
    backgroundLocator->SetDataSet(background);
    backgroundLocator->BuildLocator();
    vtkNew<vtkGenericCell> locatorCell;
    std::vector<double> weights(static_cast<std::size_t>(
        std::max<vtkIdType>(1, background->GetMaxCellSize())));
    const auto requireCoveredPoint = [&](const double point[3], const char* description) {
        if (!FindContainingCell(backgroundLocator, locatorCell, &weights, point)) {
            throw std::runtime_error("The foreground " + std::string(description) +
                                     " lies outside the physical background volume in '" +
                                     manifestPath.string() +
                                     "'. Bounds containment alone is insufficient; use a "
                                     "covering background or pre-resample the input.");
        }
    };

    double point[3] = {};
    for (vtkIdType pointId = 0; pointId < foreground->GetNumberOfPoints(); ++pointId) {
        foreground->GetPoint(pointId, point);
        requireCoveredPoint(point, "point");
    }
    vtkNew<vtkIdList> pointIds;
    for (vtkIdType cellId = 0; cellId < foreground->GetNumberOfCells(); ++cellId) {
        foreground->GetCellPoints(cellId, pointIds);
        if (pointIds->GetNumberOfIds() <= 0) {
            continue;
        }
        Vec3 centre;
        for (vtkIdType pointIndex = 0; pointIndex < pointIds->GetNumberOfIds(); ++pointIndex) {
            foreground->GetPoint(pointIds->GetId(pointIndex), point);
            centre = Add(centre, {point[0], point[1], point[2]});
        }
        centre = Scale(centre, 1.0 / static_cast<double>(pointIds->GetNumberOfIds()));
        const double centrePoint[3] = {centre.x, centre.y, centre.z};
        requireCoveredPoint(centrePoint, "cell centre");
    }
}

bool FindContainingCell(vtkStaticCellLocator* locator,
                        vtkGenericCell* cell,
                        std::vector<double>* weights,
                        const double point[3]) {
    if (locator == nullptr || cell == nullptr || weights == nullptr || point == nullptr) {
        return false;
    }
    double query[3] = {point[0], point[1], point[2]};
    double parametricCoordinates[3] = {};
    int subId = 0;
    const vtkIdType cellId = locator->FindCell(
        query, 0.0, cell, subId, parametricCoordinates, weights->data());
    return cellId >= 0;
}

void RequireForegroundRepresentableOnBackground(
    vtkDataSet* background,
    vtkDataSet* foreground,
    const std::filesystem::path& manifestPath) {
    if (background == nullptr || foreground == nullptr) {
        throw std::runtime_error("Cannot validate fixed-reference mesh-set representation in '" +
                                 manifestPath.string() + "'.");
    }
    if (foreground->GetNumberOfCells() <= 0) {
        throw std::runtime_error(
            "The foreground zone has no cells and cannot be sampled on the fixed reference "
            "mesh in '" +
            manifestPath.string() + "'.");
    }

    // The canonical output owns only the background point degrees of freedom.
    // Requiring one of those DOFs for *each* foreground cell prevents a fine
    // foreground from being partly (or wholly) invisible inside a coarse
    // background cell.  A point locator makes this conservative validation
    // proportional to nearby background points rather than all point/cell
    // pairs.
    vtkNew<vtkStaticPointLocator> backgroundPointLocator;
    backgroundPointLocator->SetDataSet(background);
    backgroundPointLocator->BuildLocator();
    vtkNew<vtkIdList> candidatePointIds;
    double point[3] = {};
    for (vtkIdType cellId = 0; cellId < foreground->GetNumberOfCells(); ++cellId) {
        vtkCell* foregroundCell = foreground->GetCell(cellId);
        if (foregroundCell == nullptr || foregroundCell->GetNumberOfPoints() <= 0) {
            throw std::runtime_error(
                "The foreground zone contains an empty cell that cannot be represented on "
                "the fixed reference mesh in '" +
                manifestPath.string() + "'.");
        }

        double bounds[6] = {};
        foregroundCell->GetBounds(bounds);
        const double centre[3] = {
            0.5 * (bounds[0] + bounds[1]),
            0.5 * (bounds[2] + bounds[3]),
            0.5 * (bounds[4] + bounds[5]),
        };
        const double dx = bounds[1] - bounds[0];
        const double dy = bounds[3] - bounds[2];
        const double dz = bounds[5] - bounds[4];
        const double scale = std::max({1.0, std::abs(centre[0]), std::abs(centre[1]),
                                       std::abs(centre[2]), dx, dy, dz});
        const double radius = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz) +
                              1.0e-10 * scale;
        candidatePointIds->Reset();
        backgroundPointLocator->FindPointsWithinRadius(radius, centre, candidatePointIds);

        std::vector<double> weights(static_cast<std::size_t>(
            std::max<vtkIdType>(1, foregroundCell->GetNumberOfPoints())));
        bool represented = false;
        for (vtkIdType candidateIndex = 0;
             candidateIndex < candidatePointIds->GetNumberOfIds(); ++candidateIndex) {
            background->GetPoint(candidatePointIds->GetId(candidateIndex), point);
            double closestPoint[3] = {};
            double parametricCoordinates[3] = {};
            double distanceSquared = 0.0;
            int subId = 0;
            if (foregroundCell->EvaluatePosition(point, closestPoint, subId,
                                                 parametricCoordinates, distanceSquared,
                                                 weights.data()) == 1) {
                represented = true;
                break;
            }
        }
        if (!represented) {
            throw std::runtime_error(
                "Foreground cell " + std::to_string(cellId) +
                " contains no fixed background degree of freedom in '" +
                manifestPath.string() +
                "'. A fine foreground cell fully embedded inside a coarse background cell "
                "cannot be represented on the fixed reference mesh; refine/pre-resample the "
                "background input.");
        }
    }
}

bool IsFloatingPointArray(const vtkDataArray* array) {
    if (array == nullptr) {
        return false;
    }
    return array->GetDataType() == VTK_FLOAT || array->GetDataType() == VTK_DOUBLE;
}

struct PointArrayPair {
    vtkDataArray* target = nullptr;
    vtkDataArray* source = nullptr;
};

bool IsFixedReferencePointArray(vtkDataArray* array) {
    if (array == nullptr || array->GetName() == nullptr) {
        return false;
    }
    const std::string lowerName = ToLower(array->GetName());
    // Nek spectral integration metadata belongs to the fixed background
    // geometry.  Copying foreground values into it would make physical
    // quadrature weights inconsistent with the canonical reference mesh.
    if (lowerName.rfind("scp_nek5000", 0) == 0) {
        return true;
    }

    std::string normalizedName;
    normalizedName.reserve(lowerName.size());
    for (unsigned char character : lowerName) {
        if (std::isalnum(character)) {
            normalizedName.push_back(static_cast<char>(character));
        }
    }
    // These names are commonly attached by CFD exporters for integration or
    // spatial measures.  They are not solution fields and must retain their
    // fixed-background values even when analysis fields are cut/blended.
    return normalizedName == "weight" || normalizedName == "weights" ||
           normalizedName == "quadratureweight" ||
           normalizedName == "quadratureweights" ||
           normalizedName == "referenceweight" ||
           normalizedName == "referenceweights" ||
           normalizedName == "integrationweight" ||
           normalizedName == "integrationweights" ||
           normalizedName == "jacobian" ||
           normalizedName == "jacobiandeterminant" || normalizedName == "detj" ||
           normalizedName == "volume" || normalizedName == "cellvolume" ||
           normalizedName == "area" || normalizedName == "cellarea" ||
           normalizedName == "measure" || normalizedName == "cellmeasure";
}

std::vector<PointArrayPair> SharedPointArrays(vtkPointData* target,
                                              vtkIdType targetPointCount,
                                              vtkPointData* source,
                                              vtkIdType sourcePointCount,
                                              bool requireFloatingPoint) {
    std::vector<PointArrayPair> arrays;
    if (target == nullptr || source == nullptr) {
        return arrays;
    }
    std::vector<std::string> unavailableFieldNames;
    for (int index = 0; index < target->GetNumberOfArrays(); ++index) {
        vtkDataArray* targetArray = target->GetArray(index);
        if (targetArray == nullptr || targetArray->GetName() == nullptr) {
            continue;
        }
        if (IsFixedReferencePointArray(targetArray)) {
            continue;
        }
        vtkDataArray* sourceArray = source->GetArray(targetArray->GetName());
        if (targetArray->GetNumberOfTuples() != targetPointCount ||
            sourceArray == nullptr ||
            sourceArray->GetNumberOfComponents() != targetArray->GetNumberOfComponents() ||
            sourceArray->GetNumberOfTuples() != sourcePointCount ||
            (requireFloatingPoint &&
             (!IsFloatingPointArray(targetArray) || !IsFloatingPointArray(sourceArray)))) {
            // A background-only analysis field would otherwise remain present
            // at covered points after the foreground replaces the flow. The
            // same is true for a truncated VTK point array: sampling skips
            // its absent foreground tuples. Do not expose either mixed-
            // reference field to RAW solvers; retain only fields complete in
            // both zones. A 50:50 blend also has no correct categorical or
            // integer result, so it retains only floating-point solution
            // fields. Fixed reference weights above are the deliberate
            // exception.
            unavailableFieldNames.emplace_back(targetArray->GetName());
            continue;
        }
        arrays.push_back({targetArray, sourceArray});
    }
    for (const std::string& name : unavailableFieldNames) {
        target->RemoveArray(name.c_str());
    }
    return arrays;
}

void AttachMultiZoneProvenance(vtkUnstructuredGrid* background,
                               vtkStaticCellLocator* foregroundLocator,
                               int backgroundZone,
                               int foregroundZone,
                               MultiZoneOverlapPolicy policy) {
    if (background == nullptr || foregroundLocator == nullptr ||
        background->GetPointData() == nullptr || background->GetCellData() == nullptr) {
        throw std::runtime_error("Cannot attach mesh-set provenance to an incomplete dataset.");
    }
    vtkPointData* pointData = background->GetPointData();
    vtkCellData* cellData = background->GetCellData();
    for (const char* name : {kMultiZoneOwnershipArray, kMultiZoneBlendWeightArray,
                             kMultiZoneBackgroundZoneArray,
                             kMultiZoneForegroundZoneArray}) {
        pointData->RemoveArray(name);
    }
    cellData->RemoveArray(kMultiZoneBackgroundCutArray);
    pointData->RemoveArray(kMultiZoneInterfaceArray);
    cellData->RemoveArray(kMultiZoneInterfaceArray);

    vtkNew<vtkIntArray> ownership;
    ownership->SetName(kMultiZoneOwnershipArray);
    ownership->SetNumberOfTuples(background->GetNumberOfPoints());
    ownership->FillValue(0);
    vtkNew<vtkDoubleArray> blendWeight;
    blendWeight->SetName(kMultiZoneBlendWeightArray);
    blendWeight->SetNumberOfTuples(background->GetNumberOfPoints());
    blendWeight->FillValue(0.0);
    vtkNew<vtkIntArray> backgroundIds;
    backgroundIds->SetName(kMultiZoneBackgroundZoneArray);
    backgroundIds->SetNumberOfTuples(background->GetNumberOfPoints());
    backgroundIds->FillValue(backgroundZone);
    vtkNew<vtkIntArray> foregroundIds;
    foregroundIds->SetName(kMultiZoneForegroundZoneArray);
    foregroundIds->SetNumberOfTuples(background->GetNumberOfPoints());
    foregroundIds->FillValue(foregroundZone);

    vtkNew<vtkGenericCell> locatorCell;
    std::vector<double> weights(static_cast<std::size_t>(
        std::max<vtkIdType>(1, foregroundLocator->GetDataSet()->GetMaxCellSize())));
    double point[3] = {};
    for (vtkIdType pointId = 0; pointId < background->GetNumberOfPoints(); ++pointId) {
        background->GetPoint(pointId, point);
        if (!FindContainingCell(foregroundLocator, locatorCell, &weights, point)) {
            continue;
        }
        ownership->SetValue(pointId, policy == MultiZoneOverlapPolicy::Blend50_50 ? 2 : 1);
        blendWeight->SetValue(pointId,
                              policy == MultiZoneOverlapPolicy::Blend50_50 ? 0.5 : 1.0);
    }
    pointData->AddArray(ownership);
    pointData->AddArray(blendWeight);
    pointData->AddArray(backgroundIds);
    pointData->AddArray(foregroundIds);

    vtkNew<vtkUnsignedCharArray> interfacePoints;
    interfacePoints->SetName(kMultiZoneInterfaceArray);
    interfacePoints->SetNumberOfTuples(background->GetNumberOfPoints());
    interfacePoints->FillValue(0);
    vtkNew<vtkUnsignedCharArray> interfaceCells;
    interfaceCells->SetName(kMultiZoneInterfaceArray);
    interfaceCells->SetNumberOfTuples(background->GetNumberOfCells());
    interfaceCells->FillValue(0);
    vtkNew<vtkIdList> pointIds;
    for (vtkIdType cellId = 0; cellId < background->GetNumberOfCells(); ++cellId) {
        background->GetCellPoints(cellId, pointIds);
        if (pointIds->GetNumberOfIds() <= 1) {
            continue;
        }
        const int firstOwnership = ownership->GetValue(pointIds->GetId(0));
        bool crossesOwnership = false;
        for (vtkIdType index = 1; index < pointIds->GetNumberOfIds(); ++index) {
            if (ownership->GetValue(pointIds->GetId(index)) != firstOwnership) {
                crossesOwnership = true;
                break;
            }
        }
        if (!crossesOwnership) {
            continue;
        }
        interfaceCells->SetValue(cellId, 1);
        for (vtkIdType index = 0; index < pointIds->GetNumberOfIds(); ++index) {
            interfacePoints->SetValue(pointIds->GetId(index), 1);
        }
    }
    pointData->AddArray(interfacePoints);
    cellData->AddArray(interfaceCells);

    vtkNew<vtkUnsignedCharArray> backgroundCut;
    backgroundCut->SetName(kMultiZoneBackgroundCutArray);
    backgroundCut->SetNumberOfTuples(background->GetNumberOfCells());
    backgroundCut->FillValue(0);
    if (policy == MultiZoneOverlapPolicy::CutBackground) {
        for (vtkIdType cellId = 0; cellId < background->GetNumberOfCells(); ++cellId) {
            background->GetCellPoints(cellId, pointIds);
            if (pointIds->GetNumberOfIds() <= 0) {
                continue;
            }
            Vec3 centre;
            for (vtkIdType index = 0; index < pointIds->GetNumberOfIds(); ++index) {
                background->GetPoint(pointIds->GetId(index), point);
                centre = Add(centre, {point[0], point[1], point[2]});
            }
            centre = Scale(centre, 1.0 / static_cast<double>(pointIds->GetNumberOfIds()));
            const double centreArray[3] = {centre.x, centre.y, centre.z};
            if (FindContainingCell(foregroundLocator, locatorCell, &weights, centreArray)) {
                backgroundCut->SetValue(cellId, 1);
            }
        }
    }
    cellData->AddArray(backgroundCut);
}

void SetMultiZoneLogicalTime(vtkDataSet* dataSet, double time) {
    if (dataSet == nullptr || !IsFiniteTime(time)) {
        throw std::runtime_error("Cannot attach a non-finite logical time to a mesh-set frame.");
    }
    vtkFieldData* fieldData = dataSet->GetFieldData();
    if (fieldData == nullptr) {
        vtkNew<vtkFieldData> createdFieldData;
        dataSet->SetFieldData(createdFieldData);
        fieldData = dataSet->GetFieldData();
    }
    if (fieldData == nullptr) {
        throw std::runtime_error("Cannot attach logical time to a mesh-set frame.");
    }
    fieldData->RemoveArray(kTimeValueArray);
    vtkNew<vtkDoubleArray> logicalTime;
    logicalTime->SetName(kTimeValueArray);
    logicalTime->SetNumberOfComponents(1);
    logicalTime->SetNumberOfTuples(1);
    logicalTime->SetValue(0, time);
    fieldData->AddArray(logicalTime);
}

vtkSmartPointer<vtkDataSet> ReadMultiZoneManifestFrame(const DataSetFrame& frame,
                                                        const ReadOptions& options) {
    const MultiZoneManifest manifest = ParseMultiZoneManifest(frame.path);
    const std::vector<DataSetFrame> firstFrames =
        EnumerateDataSetFrames(manifest.zones[0].source);
    const std::vector<DataSetFrame> secondFrames =
        EnumerateDataSetFrames(manifest.zones[1].source);
    const DataSetFrame& firstFrame =
        SelectZoneFrame(firstFrames, frame.sourceIndex, manifest.path);
    const DataSetFrame& secondFrame =
        SelectZoneFrame(secondFrames, frame.sourceIndex, manifest.path);
    ReadOptions sourceOptions = options;
    sourceOptions.failOnVtkXmlReadErrors = true;
    vtkSmartPointer<vtkDataSet> firstSource = ReadDataSet(firstFrame, sourceOptions);
    vtkSmartPointer<vtkDataSet> secondSource = ReadDataSet(secondFrame, sourceOptions);
    const int backgroundZone = frame.multiZoneBackgroundIndex == 0 ||
                                       frame.multiZoneBackgroundIndex == 1
                                   ? frame.multiZoneBackgroundIndex
                                   : ResolveMultiZoneBackground(manifest, firstFrames, secondFrames);
    const int foregroundZone = backgroundZone == 0 ? 1 : 0;
    const MultiZoneZone& selectedBackground =
        manifest.zones[static_cast<std::size_t>(backgroundZone)];
    if (selectedBackground.motion == MultiZoneMotion::Rigid &&
        (selectedBackground.angularVelocity != 0.0 ||
         !IsZero(selectedBackground.linearVelocity))) {
        throw std::runtime_error("The selected background zone must be static so the "
                                 "canonical world reference geometry remains fixed: " +
                                 manifest.path.string());
    }

    const auto requireStableForegroundLocalGeometry =
        [&](vtkDataSet* currentSource,
            const std::vector<DataSetFrame>& sourceFrames,
            const MultiZoneZone& zone,
            int zoneIndex) {
            if (frame.sourceIndex <= 0 || sourceFrames.size() <= 1 ||
                zone.velocityFrame == MultiZoneVelocityFrame::World) {
                return;
            }
            ReadOptions geometryOnly = sourceOptions;
            geometryOnly.readAllArrays = false;
            geometryOnly.attachSourceBlockIds = false;
            // Keep the source's topology-affecting options identical to the
            // current frame.  In particular, Nek5000's optional linearized
            // topology must not make an otherwise unchanged local geometry
            // look different merely because this comparison reads no fields.
            vtkSmartPointer<vtkDataSet> referenceSource = ReadDataSet(
                sourceFrames.front(), geometryOnly);
            vtkSmartPointer<vtkUnstructuredGrid> referenceLocal = MakeWorldZone(
                referenceSource, zone, zoneIndex, 0.0, false);
            vtkSmartPointer<vtkUnstructuredGrid> currentLocal = MakeWorldZone(
                currentSource, zone, zoneIndex, 0.0, false);
            RequireSameReferenceGeometry(
                referenceLocal, currentLocal, manifest.path,
                "The foreground local geometry required for zone/relative velocity");
        };
    if (foregroundZone == 0) {
        requireStableForegroundLocalGeometry(firstSource, firstFrames, manifest.zones[0], 0);
    } else {
        requireStableForegroundLocalGeometry(secondSource, secondFrames, manifest.zones[1], 1);
    }

    bool firstHasWorldVelocity = false;
    bool secondHasWorldVelocity = false;
    vtkSmartPointer<vtkUnstructuredGrid> firstWorld = MakeWorldZone(
        firstSource, manifest.zones[0], 0, frame.time, options.readAllArrays,
        &firstHasWorldVelocity);
    vtkSmartPointer<vtkUnstructuredGrid> secondWorld = MakeWorldZone(
        secondSource, manifest.zones[1], 1, frame.time, options.readAllArrays,
        &secondHasWorldVelocity);
    if (options.readAllArrays && !(firstHasWorldVelocity && secondHasWorldVelocity)) {
        // A canonical Velocity is valid only when both zones supplied a velocity
        // that was transformed into the same world reference frame. Never leave
        // a one-zone Velocity behind for a velocity solver to mistake as complete.
        firstWorld->GetPointData()->RemoveArray(kCanonicalVelocityArray);
        secondWorld->GetPointData()->RemoveArray(kCanonicalVelocityArray);
    }
    vtkSmartPointer<vtkUnstructuredGrid> background =
        backgroundZone == 0 ? firstWorld : secondWorld;
    vtkSmartPointer<vtkUnstructuredGrid> foreground =
        foregroundZone == 0 ? firstWorld : secondWorld;
    if (frame.sourceIndex > 0) {
        const std::vector<DataSetFrame>& backgroundFrames =
            backgroundZone == 0 ? firstFrames : secondFrames;
        vtkSmartPointer<vtkDataSet> referenceSource = ReadDataSet(
            backgroundFrames.front(), sourceOptions);
        vtkSmartPointer<vtkUnstructuredGrid> referenceWorld = MakeWorldZone(
            referenceSource, manifest.zones[static_cast<std::size_t>(backgroundZone)],
            backgroundZone, frame.time, options.readAllArrays);
        RequireSameReferenceGeometry(referenceWorld, background, manifest.path,
                                     "The selected background");
    }
    RequireForegroundWithinBackground(background, foreground, manifest.path);

    vtkNew<vtkStaticCellLocator> foregroundLocator;
    foregroundLocator->SetDataSet(foreground);
    foregroundLocator->BuildLocator();
    RequireForegroundRepresentableOnBackground(background, foreground, manifest.path);
    const std::vector<PointArrayPair> arrays = SharedPointArrays(
        background->GetPointData(), background->GetNumberOfPoints(),
        foreground->GetPointData(), foreground->GetNumberOfPoints(),
        manifest.overlapPolicy == MultiZoneOverlapPolicy::Blend50_50);
    vtkNew<vtkGenericCell> locatorCell;
    std::vector<double> weights(static_cast<std::size_t>(
        std::max<vtkIdType>(1, foreground->GetMaxCellSize())));
    double point[3] = {};
    for (vtkIdType targetPointId = 0; targetPointId < background->GetNumberOfPoints();
         ++targetPointId) {
        background->GetPoint(targetPointId, point);
        if (!FindContainingCell(foregroundLocator, locatorCell, &weights, point)) {
            continue;
        }
        vtkIdList* sourcePointIds = locatorCell->GetPointIds();
        if (sourcePointIds == nullptr || sourcePointIds->GetNumberOfIds() <= 0) {
            continue;
        }
        const double blend = manifest.overlapPolicy == MultiZoneOverlapPolicy::Blend50_50
                                 ? 0.5
                                 : 1.0;
        for (const PointArrayPair& pair : arrays) {
            if (pair.target == nullptr || pair.source == nullptr ||
                pair.target->GetNumberOfTuples() != background->GetNumberOfPoints() ||
                pair.source->GetNumberOfTuples() != foreground->GetNumberOfPoints() ||
                (blend < 1.0 && !IsFloatingPointArray(pair.target))) {
                continue;
            }
            for (int component = 0;
                 component < pair.target->GetNumberOfComponents(); ++component) {
                double sampled = 0.0;
                for (vtkIdType localId = 0; localId < sourcePointIds->GetNumberOfIds();
                     ++localId) {
                    sampled += weights[static_cast<std::size_t>(localId)] *
                               pair.source->GetComponent(sourcePointIds->GetId(localId), component);
                }
                const double original = pair.target->GetComponent(targetPointId, component);
                pair.target->SetComponent(targetPointId, component,
                                          (1.0 - blend) * original + blend * sampled);
            }
        }
    }
    AttachMultiZoneProvenance(background, foregroundLocator, backgroundZone,
                              foregroundZone, manifest.overlapPolicy);
    SetMultiZoneLogicalTime(background, frame.time);
    return background;
}

vtkSmartPointer<vtkDataSet> ReadOrdinaryVtk(const std::filesystem::path& path,
                                            const ReadOptions& options) {
    const std::string extension = LowerExtension(path);
    vtkSmartPointer<vtkDataObject> object;
    if (extension == ".vtkhdf") {
        object = streamcenterplus::vtkhdf::ReadDataObject(path);
    } else if (extension == ".vtk") {
        vtkNew<vtkGenericDataObjectReader> reader;
        reader->SetFileName(path.string().c_str());
        if (options.readAllArrays) {
            reader->ReadAllScalarsOn();
            reader->ReadAllVectorsOn();
            reader->ReadAllNormalsOn();
            reader->ReadAllTensorsOn();
            reader->ReadAllFieldsOn();
        }
        reader->Update();
        object = reader->GetOutputDataObject(0);
    } else if (IsTecplotExtension(extension)) {
        if (detail::IsClassicTecplotBinaryFile(path)) {
            object = detail::ReadClassicTecplotBinary(path);
        } else {
            // Tecplot users occasionally retain a .plt extension for an
            // ASCII data file.  VTK's reader handles both that convention and
            // regular .dat files, while classic binary PLT is handled above.
            vtkNew<vtkTecplotReader> reader;
            reader->SetFileName(path.string().c_str());
            reader->Update();
            object = reader->GetOutputDataObject(0);
        }
    } else {
        vtkNew<vtkXMLGenericDataObjectReader> reader;
        bool reportedReadError = false;
        vtkSmartPointer<vtkCallbackCommand> errorObserver;
        if (options.failOnVtkXmlReadErrors) {
            errorObserver = vtkSmartPointer<vtkCallbackCommand>::New();
            errorObserver->SetClientData(&reportedReadError);
            errorObserver->SetCallback(
                [](vtkObject*, unsigned long, void* clientData, void*) {
                    *static_cast<bool*>(clientData) = true;
                });
            // vtkXMLGenericDataObjectReader delegates actual array reads to a
            // concrete child reader. Its output may otherwise look usable after
            // that child rejects a short PointData array, silently turning an
            // invalid Velocity into a missing one. Mesh-set composition opts
            // into propagating either parser or child-reader error.
            reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
            reader->SetReaderErrorObserver(errorObserver);
            reader->SetParserErrorObserver(errorObserver);
        }
        reader->SetFileName(path.string().c_str());
        reader->UpdateInformation();
        if (!options.readAllArrays) {
            if (vtkDataArraySelection* selection = reader->GetPointDataArraySelection()) {
                selection->DisableAllArrays();
            }
            if (vtkDataArraySelection* selection = reader->GetCellDataArraySelection()) {
                selection->DisableAllArrays();
            }
            if (vtkDataArraySelection* selection = reader->GetColumnArraySelection()) {
                selection->DisableAllArrays();
            }
        }
        reader->Update();
        if (options.failOnVtkXmlReadErrors && reportedReadError) {
            throw std::runtime_error("Cannot read VTK XML data from '" +
                                     path.string() + "'.");
        }
        object = reader->GetOutputDataObject(0);
    }
    return CollapseDataObject(object, options, path);
}

std::vector<double> PipelineTimeSteps(vtkInformation* information) {
    std::vector<double> values;
    if (information == nullptr) {
        return values;
    }
    auto* key = vtkStreamingDemandDrivenPipeline::TIME_STEPS();
    if (!information->Has(key)) {
        return values;
    }
    const int count = information->Length(key);
    values.reserve(static_cast<std::size_t>(std::max(0, count)));
    for (int i = 0; i < count; ++i) {
        const double value = information->Get(key, i);
        if (IsFiniteTime(value)) {
            values.push_back(value);
        }
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

std::optional<std::filesystem::path> ResolveOpenFoamCase(
    const std::filesystem::path& sourcePath) {
    if (IsRegularFileWithExtension(sourcePath, ".foam")) {
        return std::filesystem::absolute(sourcePath);
    }
    if (!std::filesystem::is_directory(sourcePath)) {
        return std::nullopt;
    }
    std::vector<std::filesystem::path> markers;
    for (const auto& entry : std::filesystem::directory_iterator(sourcePath)) {
        if (entry.is_regular_file() && LowerExtension(entry.path()) == ".foam") {
            markers.push_back(entry.path());
        }
    }
    std::sort(markers.begin(), markers.end());
    if (!markers.empty()) {
        return std::filesystem::absolute(markers.front());
    }

    const bool recognizable =
        std::filesystem::is_regular_file(sourcePath / "system" / "controlDict") ||
        std::filesystem::is_directory(sourcePath / "constant" / "polyMesh") ||
        std::filesystem::is_directory(sourcePath / "processor0");
    if (!recognizable) {
        return std::nullopt;
    }
    // VTK treats the .foam path as a marker and reads the case from its parent;
    // no marker is created in the user's case directory.
    return std::filesystem::absolute(sourcePath /
                                     (sourcePath.filename().string() + ".foam"));
}

bool IsDecomposedOpenFoamCase(const std::filesystem::path& marker) {
    const std::filesystem::path root = marker.parent_path();
    return std::filesystem::is_directory(root / "processor0");
}

vtkSmartPointer<vtkPOpenFOAMReader> CreateOpenFoamReader(
    const std::filesystem::path& marker,
    bool readAllArrays) {
    vtkSmartPointer<vtkPOpenFOAMReader> reader =
        vtkSmartPointer<vtkPOpenFOAMReader>::New();
    reader->SetCaseType(IsDecomposedOpenFoamCase(marker)
                            ? vtkPOpenFOAMReader::DECOMPOSED_CASE
                            : vtkPOpenFOAMReader::RECONSTRUCTED_CASE);
    reader->SetFileName(marker.string().c_str());
    reader->SetSkipZeroTime(false);
    reader->SetListTimeStepsByControlDict(false);
    reader->SetCreateCellToPoint(false);
    reader->SetCacheMesh(true);
    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->DisableAllPatchArrays();
    if (readAllArrays) {
        reader->EnableAllCellArrays();
        reader->EnableAllPointArrays();
    } else {
        reader->DisableAllCellArrays();
        reader->DisableAllPointArrays();
    }
    reader->UpdateInformation();
    reader->SetReadZones(false);
    reader->DisableAllLagrangianArrays();
    reader->DisableAllPatchArrays();
    for (int i = 0; i < reader->GetNumberOfPatchArrays(); ++i) {
        const char* rawName = reader->GetPatchArrayName(i);
        const std::string name = rawName == nullptr ? std::string() : ToLower(rawName);
        if (name == "internalmesh" || name.find("internal") != std::string::npos) {
            reader->SetPatchArrayStatus(rawName, 1);
        }
    }
    return reader;
}

std::vector<double> OpenFoamTimeSteps(vtkOpenFOAMReader* reader) {
    std::vector<double> values;
    if (reader == nullptr) {
        return values;
    }
    if (vtkDoubleArray* times = reader->GetTimeValues()) {
        values.reserve(static_cast<std::size_t>(times->GetNumberOfTuples()));
        for (vtkIdType i = 0; i < times->GetNumberOfTuples(); ++i) {
            const double value = times->GetValue(i);
            if (IsFiniteTime(value)) {
                values.push_back(value);
            }
        }
    }
    if (values.empty()) {
        values = PipelineTimeSteps(reader->GetOutputInformation(0));
    }
    if (values.empty()) {
        values.push_back(0.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

vtkSmartPointer<vtkDataSet> ReadOpenFoam(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    vtkSmartPointer<vtkPOpenFOAMReader> reader =
        CreateOpenFoamReader(frame.path, options.readAllArrays);
    const double time = frame.readerTime.value_or(frame.time);
    reader->SetTimeValue(time);
    if (reader->UpdateTimeStep(time) == 0) {
        throw std::runtime_error("VTK failed to read OpenFOAM time " +
                                 std::to_string(time) + " from " + frame.path.string());
    }
    return CollapseDataObject(reader->GetOutputDataObject(0), options, frame.path);
}

vtkSmartPointer<vtkGenericEnSightReader> CreateEnSightGoldReader(
    const std::filesystem::path& casePath,
    bool readAllArrays) {
    vtkSmartPointer<vtkGenericEnSightReader> reader =
        vtkSmartPointer<vtkGenericEnSightReader>::New();
    reader->SetCaseFileName(casePath.string().c_str());
    if (readAllArrays) {
        reader->ReadAllVariablesOn();
    } else {
        reader->ReadAllVariablesOff();
    }
    reader->UpdateInformation();
    const int version = reader->GetEnSightVersion();
    if (version != vtkGenericEnSightReader::ENSIGHT_GOLD &&
        version != vtkGenericEnSightReader::ENSIGHT_GOLD_BINARY) {
        throw std::runtime_error("The .case source is not EnSight Gold: " +
                                 casePath.string());
    }
    if (readAllArrays) {
        reader->ReadAllVariablesOn();
        for (int i = 0; i < reader->GetNumberOfPointArrays(); ++i) {
            if (const char* name = reader->GetPointArrayName(i)) {
                reader->SetPointArrayStatus(name, 1);
            }
        }
        for (int i = 0; i < reader->GetNumberOfCellArrays(); ++i) {
            if (const char* name = reader->GetCellArrayName(i)) {
                reader->SetCellArrayStatus(name, 1);
            }
        }
    }
    return reader;
}

std::vector<double> EnSightGoldTimeSteps(vtkGenericEnSightReader* reader) {
    std::vector<double> values =
        reader == nullptr ? std::vector<double>()
                          : PipelineTimeSteps(reader->GetOutputInformation(0));
    if (values.empty() && reader != nullptr) {
        if (vtkDataArrayCollection* sets = reader->GetTimeSets()) {
            for (int setIndex = 0; setIndex < sets->GetNumberOfItems(); ++setIndex) {
                vtkDataArray* array = sets->GetItem(setIndex);
                if (array == nullptr) {
                    continue;
                }
                for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
                    const double value = array->GetComponent(i, 0);
                    if (IsFiniteTime(value)) {
                        values.push_back(value);
                    }
                }
            }
        }
    }
    if (values.empty()) {
        values.push_back(0.0);
    }
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end(), [](double lhs, double rhs) {
                     return std::abs(lhs - rhs) <= 1.0e-12;
                 }),
                 values.end());
    return values;
}

vtkSmartPointer<vtkDataSet> ReadEnSightGold(const DataSetFrame& frame,
                                           const ReadOptions& options) {
    vtkSmartPointer<vtkGenericEnSightReader> reader =
        CreateEnSightGoldReader(frame.path, options.readAllArrays);
    const double time = frame.readerTime.value_or(frame.time);
    reader->SetTimeValue(time);
    if (reader->UpdateTimeStep(time) == 0) {
        throw std::runtime_error("VTK failed to read EnSight Gold time " +
                                 std::to_string(time) + " from " + frame.path.string());
    }
    return CollapseDataObject(reader->GetOutputDataObject(0), options, frame.path);
}

bool HostIsLittleEndian() {
    const std::uint16_t value = 1;
    return *reinterpret_cast<const unsigned char*>(&value) == 1;
}

template <typename T>
T ByteSwapValue(T value) {
    std::array<unsigned char, sizeof(T)> bytes{};
    std::memcpy(bytes.data(), &value, sizeof(T));
    std::reverse(bytes.begin(), bytes.end());
    std::memcpy(&value, bytes.data(), sizeof(T));
    return value;
}

struct NekDescriptor {
    std::filesystem::path path;
    std::filesystem::path caseDirectory;
    std::string fileTemplate;
    int firstTimeStep = 0;
    int numberOfTimeSteps = 0;
    int integerConversions = 0;
};

struct NekFieldMeta {
    std::filesystem::path path;
    int wordSize = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int elementCount = 0;
    int globalElementCount = 0;
    int fileIndex = 0;
    int fileCount = 1;
    int nxyz = 0;
    int dimensions = 0;
    double time = 0.0;
    std::string code;
    bool littleEndian = true;
    std::uint64_t dataBaseOffset = 0;
    std::vector<std::string> labels;

    std::size_t pointCount() const {
        return static_cast<std::size_t>(elementCount) *
               static_cast<std::size_t>(nxyz);
    }
};

int ValidateNekPrintfTemplate(const std::string& format) {
    int conversions = 0;
    for (std::size_t i = 0; i < format.size(); ++i) {
        if (format[i] != '%') {
            continue;
        }
        if (i + 1 < format.size() && format[i + 1] == '%') {
            ++i;
            continue;
        }
        ++i;
        while (i < format.size() && std::string("-+ #0").find(format[i]) != std::string::npos) {
            ++i;
        }
        while (i < format.size() && std::isdigit(static_cast<unsigned char>(format[i]))) {
            ++i;
        }
        if (i < format.size() && format[i] == '.') {
            ++i;
            while (i < format.size() && std::isdigit(static_cast<unsigned char>(format[i]))) {
                ++i;
            }
        }
        if (i < format.size() &&
            std::string("hljztL").find(format[i]) != std::string::npos) {
            throw std::runtime_error(
                "Nek5000 filetemplate integer conversions must use the default int length.");
        }
        if (i >= format.size() || std::string("diuoxX").find(format[i]) == std::string::npos) {
            throw std::runtime_error(
                "Nek5000 filetemplate contains a non-integer or unsafe printf conversion.");
        }
        ++conversions;
    }
    if (conversions < 1 || conversions > 3) {
        throw std::runtime_error(
            "Nek5000 filetemplate must contain one timestep conversion, "
            "file/timestep conversions, or directory/file/timestep conversions.");
    }
    return conversions;
}

NekDescriptor ParseNekDescriptor(const std::filesystem::path& descriptorPath) {
    std::ifstream input(descriptorPath);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 descriptor: " +
                                 descriptorPath.string());
    }
    NekDescriptor descriptor;
    descriptor.path = std::filesystem::absolute(descriptorPath);
    descriptor.caseDirectory = descriptor.path.parent_path();
    bool hasTemplate = false;
    bool hasFirst = false;
    bool hasCount = false;
    std::string line;
    while (std::getline(input, line)) {
        line = StripComment(line);
        if (line.empty()) {
            continue;
        }
        const std::size_t colon = line.find(':');
        if (colon == std::string::npos) {
            continue;
        }
        const std::string key = ToLower(TrimAscii(line.substr(0, colon)));
        const std::string value = TrimAscii(line.substr(colon + 1));
        if (key == "filetemplate") {
            descriptor.fileTemplate = value;
            hasTemplate = !value.empty();
        } else if (key == "firsttimestep") {
            descriptor.firstTimeStep = std::stoi(value);
            hasFirst = true;
        } else if (key == "numtimesteps") {
            descriptor.numberOfTimeSteps = std::stoi(value);
            hasCount = true;
        } else if (key == "nek5000" || key == "endian" || key == "version") {
            continue;
        } else {
            throw std::runtime_error("Unsupported Nek5000 descriptor tag: " + key);
        }
    }
    if (!hasTemplate || !hasFirst || !hasCount || descriptor.numberOfTimeSteps <= 0) {
        throw std::runtime_error(
            "Nek5000 descriptor must define filetemplate, firsttimestep, and numtimesteps.");
    }
    descriptor.integerConversions = ValidateNekPrintfTemplate(descriptor.fileTemplate);
    return descriptor;
}

std::filesystem::path FormatNekFieldPath(const NekDescriptor& descriptor,
                                         int relativeFrameIndex,
                                         int fileIndex) {
    const int timeStep = descriptor.firstTimeStep + relativeFrameIndex;
    std::array<char, 4096> buffer{};
    int written = 0;
    if (descriptor.integerConversions == 3) {
        written = std::snprintf(buffer.data(), buffer.size(),
                                descriptor.fileTemplate.c_str(),
                                fileIndex, fileIndex, timeStep);
    } else if (descriptor.integerConversions == 2) {
        written = std::snprintf(buffer.data(), buffer.size(),
                                descriptor.fileTemplate.c_str(), fileIndex, timeStep);
    } else {
        written = std::snprintf(buffer.data(), buffer.size(),
                                descriptor.fileTemplate.c_str(), timeStep);
    }
    if (written < 0 || static_cast<std::size_t>(written) >= buffer.size()) {
        throw std::runtime_error("Nek5000 filetemplate produced an invalid path.");
    }
    std::filesystem::path path(buffer.data());
    if (path.is_relative()) {
        path = descriptor.caseDirectory / path;
    }
    return path;
}

float DecodeEndianTag(const std::array<char, 4>& bytes, bool littleEndian) {
    std::array<char, 4> ordered = bytes;
    if (littleEndian != HostIsLittleEndian()) {
        std::reverse(ordered.begin(), ordered.end());
    }
    float value = 0.0f;
    std::memcpy(&value, ordered.data(), sizeof(value));
    return value;
}

std::vector<std::string> SplitWhitespace(const std::string& value) {
    std::istringstream stream(value);
    std::vector<std::string> parts;
    std::string token;
    while (stream >> token) {
        parts.push_back(token);
    }
    return parts;
}

std::vector<std::string> NekLabelsForCode(std::string code, int dimensions) {
    std::transform(code.begin(), code.end(), code.begin(), [](unsigned char ch) {
        return static_cast<char>(std::toupper(ch));
    });
    std::vector<std::string> labels;
    if (code.find('X') != std::string::npos) {
        labels.push_back("x");
        labels.push_back("y");
        if (dimensions == 3) {
            labels.push_back("z");
        }
    }
    if (code.find('U') != std::string::npos) {
        labels.push_back("u");
        labels.push_back("v");
        if (dimensions == 3) {
            labels.push_back("w");
        }
    }
    if (code.find('P') != std::string::npos) {
        labels.push_back("p");
    }
    if (code.find('T') != std::string::npos) {
        labels.push_back("t");
    }
    const std::size_t scalar = code.find('S');
    if (scalar != std::string::npos) {
        std::size_t end = scalar + 1;
        while (end < code.size() && std::isdigit(static_cast<unsigned char>(code[end]))) {
            ++end;
        }
        const std::string countToken = code.substr(scalar + 1, end - scalar - 1);
        const int count = countToken.empty() ? 1 : std::stoi(countToken);
        for (int i = 0; i < count; ++i) {
            std::ostringstream name;
            name << 's' << std::setw(2) << std::setfill('0') << (i + 1);
            labels.push_back(name.str());
        }
    }
    return labels;
}

NekFieldMeta ReadNekHeader(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 field file: " + path.string());
    }
    std::array<char, 132> headerBytes{};
    std::array<char, 4> tagBytes{};
    input.read(headerBytes.data(), static_cast<std::streamsize>(headerBytes.size()));
    input.read(tagBytes.data(), static_cast<std::streamsize>(tagBytes.size()));
    if (!input) {
        throw std::runtime_error("Cannot read Nek5000 field header: " + path.string());
    }
    const std::vector<std::string> parts =
        SplitWhitespace(std::string(headerBytes.data(), headerBytes.size()));
    if (parts.size() < 12 || parts[0] != "#std") {
        throw std::runtime_error("Not a supported Nek5000 #std field file: " + path.string());
    }

    NekFieldMeta meta;
    meta.path = path;
    meta.wordSize = std::stoi(parts[1]);
    meta.nx = std::stoi(parts[2]);
    meta.ny = std::stoi(parts[3]);
    meta.nz = std::stoi(parts[4]);
    meta.elementCount = std::stoi(parts[5]);
    meta.globalElementCount = std::stoi(parts[6]);
    std::string timeToken = parts[7];
    std::replace(timeToken.begin(), timeToken.end(), 'D', 'E');
    std::replace(timeToken.begin(), timeToken.end(), 'd', 'E');
    meta.time = std::stod(timeToken);
    meta.fileIndex = std::stoi(parts[9]);
    meta.fileCount = std::max(1, std::stoi(parts[10]));
    meta.code = parts[11];
    meta.dimensions = meta.nz > 1 ? 3 : 2;
    meta.nxyz = meta.nx * meta.ny * meta.nz;
    const float littleTag = DecodeEndianTag(tagBytes, true);
    const float bigTag = DecodeEndianTag(tagBytes, false);
    if (std::abs(littleTag - 6.54321f) < 1.0e-4f) {
        meta.littleEndian = true;
    } else if (std::abs(bigTag - 6.54321f) < 1.0e-4f) {
        meta.littleEndian = false;
    } else {
        throw std::runtime_error("Cannot determine Nek5000 field endianness: " +
                                 path.string());
    }
    if (meta.wordSize != 4 && meta.wordSize != 8) {
        throw std::runtime_error("Nek5000 field word size must be 4 or 8 bytes: " +
                                 path.string());
    }
    if (meta.nx < 2 || meta.ny < 2 || meta.nz < 1 || meta.elementCount <= 0 ||
        meta.globalElementCount <= 0 || !IsFiniteTime(meta.time)) {
        throw std::runtime_error("Invalid Nek5000 field dimensions or metadata: " +
                                 path.string());
    }
    meta.dataBaseOffset = 136ull +
                          4ull * static_cast<std::uint64_t>(meta.elementCount);
    meta.labels = NekLabelsForCode(meta.code, meta.dimensions);
    return meta;
}

void ValidateNekPieceCompatible(const NekFieldMeta& reference,
                                const NekFieldMeta& current,
                                bool requireCode) {
    if ((requireCode && reference.wordSize != current.wordSize) ||
        reference.nx != current.nx ||
        reference.ny != current.ny || reference.nz != current.nz ||
        reference.globalElementCount != current.globalElementCount ||
        (requireCode && reference.fileCount != current.fileCount)) {
        throw std::runtime_error("Nek5000 field pieces have incompatible dimensions: " +
                                 current.path.string());
    }
    if (requireCode && reference.code != current.code) {
        throw std::runtime_error("Nek5000 field pieces have different rdcode values: " +
                                 current.path.string());
    }
}

std::vector<NekFieldMeta> ReadNekPieceSet(const NekDescriptor& descriptor,
                                          int relativeFrameIndex) {
    NekFieldMeta first =
        ReadNekHeader(FormatNekFieldPath(descriptor, relativeFrameIndex, 0));
    if (first.fileCount > 1 && descriptor.integerConversions < 2) {
        throw std::runtime_error(
            "Multi-file Nek5000 output requires filetemplate to contain file and timestep conversions.");
    }
    std::vector<NekFieldMeta> pieces;
    pieces.reserve(static_cast<std::size_t>(first.fileCount));
    pieces.push_back(first);
    for (int fileIndex = 1; fileIndex < first.fileCount; ++fileIndex) {
        NekFieldMeta piece =
            ReadNekHeader(FormatNekFieldPath(descriptor, relativeFrameIndex, fileIndex));
        ValidateNekPieceCompatible(first, piece, true);
        pieces.push_back(std::move(piece));
    }
    std::size_t totalElements = 0;
    for (const NekFieldMeta& piece : pieces) {
        totalElements += static_cast<std::size_t>(piece.elementCount);
    }
    if (totalElements != static_cast<std::size_t>(first.globalElementCount)) {
        throw std::runtime_error("Nek5000 field pieces do not cover the advertised global elements.");
    }
    return pieces;
}

std::vector<int> ReadNekElementMap(const NekFieldMeta& meta) {
    std::ifstream input(meta.path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 element map: " + meta.path.string());
    }
    input.seekg(136, std::ios::beg);
    std::vector<std::int32_t> raw(static_cast<std::size_t>(meta.elementCount));
    input.read(reinterpret_cast<char*>(raw.data()),
               static_cast<std::streamsize>(raw.size() * sizeof(std::int32_t)));
    if (!input) {
        throw std::runtime_error("Cannot read Nek5000 element map: " + meta.path.string());
    }
    std::vector<int> ids(raw.size(), 0);
    for (std::size_t i = 0; i < raw.size(); ++i) {
        std::int32_t value = raw[i];
        if (meta.littleEndian != HostIsLittleEndian()) {
            value = ByteSwapValue(value);
        }
        ids[i] = static_cast<int>(value);
    }
    return ids;
}

void ReadNekRawValues(std::ifstream& input,
                      const NekFieldMeta& meta,
                      std::uint64_t byteOffset,
                      std::size_t count,
                      double* destination) {
    input.seekg(static_cast<std::streamoff>(byteOffset), std::ios::beg);
    if (!input) {
        throw std::runtime_error("Cannot seek Nek5000 field data: " + meta.path.string());
    }
    constexpr std::size_t kChunk = 1u << 20;
    std::size_t written = 0;
    if (meta.wordSize == 4) {
        std::vector<float> buffer;
        while (written < count) {
            const std::size_t current = std::min(kChunk, count - written);
            buffer.resize(current);
            input.read(reinterpret_cast<char*>(buffer.data()),
                       static_cast<std::streamsize>(current * sizeof(float)));
            if (!input) {
                throw std::runtime_error("Cannot read Nek5000 float field data: " +
                                         meta.path.string());
            }
            for (std::size_t i = 0; i < current; ++i) {
                float value = buffer[i];
                if (meta.littleEndian != HostIsLittleEndian()) {
                    value = ByteSwapValue(value);
                }
                destination[written + i] = static_cast<double>(value);
            }
            written += current;
        }
        return;
    }
    std::vector<double> buffer;
    while (written < count) {
        const std::size_t current = std::min(kChunk, count - written);
        buffer.resize(current);
        input.read(reinterpret_cast<char*>(buffer.data()),
                   static_cast<std::streamsize>(current * sizeof(double)));
        if (!input) {
            throw std::runtime_error("Cannot read Nek5000 double field data: " +
                                     meta.path.string());
        }
        for (std::size_t i = 0; i < current; ++i) {
            double value = buffer[i];
            if (meta.littleEndian != HostIsLittleEndian()) {
                value = ByteSwapValue(value);
            }
            destination[written + i] = value;
        }
        written += current;
    }
}

struct NekVectorLayout {
    std::uint64_t blockOffset = 0;
    int component = 0;
};

std::optional<NekVectorLayout> ResolveNekVectorLayout(const NekFieldMeta& meta,
                                                       const std::string& label) {
    const std::string lower = ToLower(label);
    const int dimensions = meta.dimensions;
    if (lower == "x" || lower == "y" || lower == "z") {
        if (meta.code.find('X') == std::string::npos ||
            (lower == "z" && dimensions != 3)) {
            return std::nullopt;
        }
        return NekVectorLayout{meta.dataBaseOffset,
                               lower == "x" ? 0 : (lower == "y" ? 1 : 2)};
    }
    if (lower == "u" || lower == "v" || lower == "w") {
        if (meta.code.find('U') == std::string::npos ||
            (lower == "w" && dimensions != 3)) {
            return std::nullopt;
        }
        const std::uint64_t coordinateComponents =
            meta.code.find('X') == std::string::npos
                ? 0ull
                : static_cast<std::uint64_t>(dimensions);
        const std::uint64_t blockOffset =
            meta.dataBaseOffset + coordinateComponents *
                                      static_cast<std::uint64_t>(meta.elementCount) *
                                      static_cast<std::uint64_t>(meta.nxyz) *
                                      static_cast<std::uint64_t>(meta.wordSize);
        return NekVectorLayout{blockOffset,
                               lower == "u" ? 0 : (lower == "v" ? 1 : 2)};
    }
    return std::nullopt;
}

int NekLabelIndex(const NekFieldMeta& meta, const std::string& label) {
    const std::string lower = ToLower(label);
    for (std::size_t i = 0; i < meta.labels.size(); ++i) {
        if (meta.labels[i] == lower) {
            return static_cast<int>(i);
        }
    }
    throw std::runtime_error("Nek5000 field has no component '" + label +
                             "': " + meta.path.string());
}

void ReadNekComponent(const NekFieldMeta& meta,
                      const std::string& label,
                      std::vector<double>* values) {
    if (values == nullptr) {
        throw std::runtime_error("Internal error: null Nek5000 component destination.");
    }
    values->assign(meta.pointCount(), 0.0);
    std::ifstream input(meta.path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek5000 field data: " + meta.path.string());
    }
    if (const std::optional<NekVectorLayout> layout =
            ResolveNekVectorLayout(meta, label)) {
        const std::size_t nxyz = static_cast<std::size_t>(meta.nxyz);
        const std::size_t valuesPerElement =
            static_cast<std::size_t>(meta.dimensions) * nxyz;
        std::vector<double> raw(valuesPerElement);
        for (int element = 0; element < meta.elementCount; ++element) {
            const std::uint64_t offset =
                layout->blockOffset +
                static_cast<std::uint64_t>(element) *
                    static_cast<std::uint64_t>(valuesPerElement) *
                    static_cast<std::uint64_t>(meta.wordSize);
            ReadNekRawValues(input, meta, offset, valuesPerElement, raw.data());
            const std::size_t source =
                static_cast<std::size_t>(layout->component) * nxyz;
            const std::size_t target = static_cast<std::size_t>(element) * nxyz;
            std::copy_n(raw.data() + source, nxyz, values->data() + target);
        }
        return;
    }
    const int componentIndex = NekLabelIndex(meta, label);
    const std::uint64_t offset =
        meta.dataBaseOffset + static_cast<std::uint64_t>(componentIndex) *
                                  static_cast<std::uint64_t>(meta.elementCount) *
                                  static_cast<std::uint64_t>(meta.nxyz) *
                                  static_cast<std::uint64_t>(meta.wordSize);
    ReadNekRawValues(input, meta, offset, meta.pointCount(), values->data());
}

struct NekPieceLookup {
    std::vector<std::vector<int>> ids;
    std::vector<std::pair<std::size_t, std::size_t>> outputToPieceElement;
    std::vector<int> outputElementIds;
};

NekPieceLookup BuildNekPieceLookup(const std::vector<NekFieldMeta>& dataPieces) {
    NekPieceLookup lookup;
    struct PieceElement {
        int id = 0;
        std::size_t piece = 0;
        std::size_t local = 0;
    };
    std::vector<PieceElement> elements;
    for (std::size_t piece = 0; piece < dataPieces.size(); ++piece) {
        lookup.ids.push_back(ReadNekElementMap(dataPieces[piece]));
        for (std::size_t local = 0; local < lookup.ids.back().size(); ++local) {
            elements.push_back({lookup.ids.back()[local], piece, local});
        }
    }
    std::sort(elements.begin(), elements.end(), [](const PieceElement& lhs,
                                                    const PieceElement& rhs) {
        return lhs.id < rhs.id;
    });
    lookup.outputToPieceElement.reserve(elements.size());
    lookup.outputElementIds.reserve(elements.size());
    for (std::size_t output = 0; output < elements.size(); ++output) {
        if (output > 0 && elements[output - 1].id == elements[output].id) {
            throw std::runtime_error("Duplicate Nek5000 global element id " +
                                     std::to_string(elements[output].id) + ".");
        }
        lookup.outputToPieceElement.push_back(
            {elements[output].piece, elements[output].local});
        lookup.outputElementIds.push_back(elements[output].id);
    }
    return lookup;
}

std::unordered_map<int, std::pair<std::size_t, std::size_t>> BuildNekElementMap(
    const std::vector<NekFieldMeta>& pieces) {
    std::unordered_map<int, std::pair<std::size_t, std::size_t>> result;
    for (std::size_t piece = 0; piece < pieces.size(); ++piece) {
        const std::vector<int> ids = ReadNekElementMap(pieces[piece]);
        for (std::size_t local = 0; local < ids.size(); ++local) {
            if (!result.emplace(ids[local], std::make_pair(piece, local)).second) {
                throw std::runtime_error("Duplicate Nek5000 coordinate element id " +
                                         std::to_string(ids[local]) + ".");
            }
        }
    }
    return result;
}

void ScatterNekComponent(const std::vector<NekFieldMeta>& pieces,
                         const NekPieceLookup& outputLookup,
                         const std::string& label,
                         std::vector<double>* output) {
    if (output == nullptr || pieces.empty()) {
        throw std::runtime_error("Internal error while assembling Nek5000 component.");
    }
    const std::size_t nxyz = static_cast<std::size_t>(pieces.front().nxyz);
    output->assign(outputLookup.outputToPieceElement.size() * nxyz, 0.0);
    std::vector<std::vector<double>> pieceValues(pieces.size());
    for (std::size_t piece = 0; piece < pieces.size(); ++piece) {
        ReadNekComponent(pieces[piece], label, &pieceValues[piece]);
    }
    for (std::size_t outputElement = 0;
         outputElement < outputLookup.outputToPieceElement.size(); ++outputElement) {
        const auto [piece, local] = outputLookup.outputToPieceElement[outputElement];
        std::copy_n(pieceValues[piece].data() + local * nxyz, nxyz,
                    output->data() + outputElement * nxyz);
    }
}

void ScatterNekCoordinateComponent(const std::vector<NekFieldMeta>& coordinatePieces,
                                   const NekPieceLookup& outputLookup,
                                   const std::string& label,
                                   std::vector<double>* output) {
    if (output == nullptr || coordinatePieces.empty()) {
        throw std::runtime_error("Internal error while assembling Nek5000 coordinates.");
    }
    const std::size_t nxyz = static_cast<std::size_t>(coordinatePieces.front().nxyz);
    const auto coordinateMap = BuildNekElementMap(coordinatePieces);
    std::vector<std::vector<double>> pieceValues(coordinatePieces.size());
    for (std::size_t piece = 0; piece < coordinatePieces.size(); ++piece) {
        ReadNekComponent(coordinatePieces[piece], label, &pieceValues[piece]);
    }
    output->assign(outputLookup.outputElementIds.size() * nxyz, 0.0);
    for (std::size_t outputElement = 0;
         outputElement < outputLookup.outputElementIds.size(); ++outputElement) {
        const int id = outputLookup.outputElementIds[outputElement];
        const auto found = coordinateMap.find(id);
        if (found == coordinateMap.end()) {
            throw std::runtime_error("Nek5000 coordinate field is missing global element id " +
                                     std::to_string(id) + ".");
        }
        const auto [piece, local] = found->second;
        std::copy_n(pieceValues[piece].data() + local * nxyz, nxyz,
                    output->data() + outputElement * nxyz);
    }
}

double LegendreP(int order, double x) {
    if (order == 0) {
        return 1.0;
    }
    if (order == 1) {
        return x;
    }
    double previous = 1.0;
    double current = x;
    for (int k = 2; k <= order; ++k) {
        const double next =
            ((2.0 * static_cast<double>(k) - 1.0) * x * current -
             (static_cast<double>(k) - 1.0) * previous) /
            static_cast<double>(k);
        previous = current;
        current = next;
    }
    return current;
}

double LegendreDerivative(int order, double x) {
    if (order == 0) {
        return 0.0;
    }
    return static_cast<double>(order) *
           (x * LegendreP(order, x) - LegendreP(order - 1, x)) /
           (x * x - 1.0);
}

std::vector<double> GllNodes(int count) {
    if (count < 2) {
        throw std::runtime_error("GLL node count must be at least two.");
    }
    constexpr double kPi = 3.141592653589793238462643383279502884;
    const int order = count - 1;
    std::vector<double> nodes(static_cast<std::size_t>(count), 0.0);
    nodes.front() = -1.0;
    nodes.back() = 1.0;
    for (int i = 1; i < order; ++i) {
        double x = -std::cos(kPi * static_cast<double>(i) /
                             static_cast<double>(order));
        for (int iteration = 0; iteration < 60; ++iteration) {
            const double p = LegendreP(order, x);
            const double dp = LegendreDerivative(order, x);
            const double ddp =
                (2.0 * x * dp - static_cast<double>(order) *
                                      static_cast<double>(order + 1) * p) /
                (1.0 - x * x);
            const double delta = dp / ddp;
            x -= delta;
            if (std::abs(delta) < 2.0e-15) {
                break;
            }
        }
        nodes[static_cast<std::size_t>(i)] = x;
    }
    return nodes;
}

struct OutwardInterval {
    double lower = 0.0;
    double upper = 0.0;
};

double RoundDown(double value) {
    return std::nextafter(value, -std::numeric_limits<double>::infinity());
}

double RoundUp(double value) {
    return std::nextafter(value, std::numeric_limits<double>::infinity());
}

OutwardInterval IntervalPoint(double value) {
    return {value, value};
}

OutwardInterval IntervalAdd(const OutwardInterval& lhs,
                            const OutwardInterval& rhs) {
    return {RoundDown(lhs.lower + rhs.lower),
            RoundUp(lhs.upper + rhs.upper)};
}

OutwardInterval IntervalNegate(const OutwardInterval& value) {
    return {-value.upper, -value.lower};
}

OutwardInterval IntervalSubtract(const OutwardInterval& lhs,
                                 const OutwardInterval& rhs) {
    return IntervalAdd(lhs, IntervalNegate(rhs));
}

OutwardInterval IntervalMultiply(const OutwardInterval& lhs,
                                 const OutwardInterval& rhs) {
    const std::array<double, 4> products = {
        lhs.lower * rhs.lower, lhs.lower * rhs.upper,
        lhs.upper * rhs.lower, lhs.upper * rhs.upper};
    const auto minimum = std::min_element(products.begin(), products.end());
    const auto maximum = std::max_element(products.begin(), products.end());
    return {RoundDown(*minimum), RoundUp(*maximum)};
}

OutwardInterval IntervalDivide(const OutwardInterval& numerator,
                               const OutwardInterval& denominator) {
    if (denominator.lower <= 0.0 && denominator.upper >= 0.0) {
        throw std::runtime_error(
            "Cannot form a conservative Bernstein bound through a zero interval.");
    }
    const std::array<double, 4> quotients = {
        numerator.lower / denominator.lower,
        numerator.lower / denominator.upper,
        numerator.upper / denominator.lower,
        numerator.upper / denominator.upper};
    const auto minimum = std::min_element(quotients.begin(), quotients.end());
    const auto maximum = std::max_element(quotients.begin(), quotients.end());
    return {RoundDown(*minimum), RoundUp(*maximum)};
}

std::vector<OutwardInterval> MultiplyBernsteinByLinear(
    const std::vector<OutwardInterval>& polynomial,
    const std::array<OutwardInterval, 2>& linear) {
    if (polynomial.empty()) {
        throw std::runtime_error("Cannot multiply an empty Bernstein polynomial.");
    }
    const std::size_t degree = polynomial.size() - 1;
    std::vector<OutwardInterval> product(degree + 2);
    for (std::size_t k = 0; k <= degree + 1; ++k) {
        bool initialized = false;
        OutwardInterval coefficient{};
        if (k <= degree) {
            const OutwardInterval ratio = IntervalDivide(
                IntervalPoint(static_cast<double>(degree + 1 - k)),
                IntervalPoint(static_cast<double>(degree + 1)));
            coefficient = IntervalMultiply(
                IntervalMultiply(polynomial[k], linear[0]), ratio);
            initialized = true;
        }
        if (k > 0) {
            const OutwardInterval ratio = IntervalDivide(
                IntervalPoint(static_cast<double>(k)),
                IntervalPoint(static_cast<double>(degree + 1)));
            const OutwardInterval contribution = IntervalMultiply(
                IntervalMultiply(polynomial[k - 1], linear[1]), ratio);
            coefficient = initialized ? IntervalAdd(coefficient, contribution)
                                      : contribution;
        }
        product[k] = coefficient;
    }
    return product;
}

// On every open interval between adjacent interpolation nodes, the signs of
// all Lagrange basis functions are fixed. Therefore the Lebesgue function is
// one ordinary polynomial there. We construct each signed basis directly in
// the Bernstein basis on two half-intervals. Outward-rounded interval
// arithmetic encloses every control coefficient, and the Bernstein convex-hull
// property then supplies a rigorous upper bound without sampling the function.
double ConservativeLebesgueBound(const std::vector<double>& interpolationNodes) {
    const std::size_t count = interpolationNodes.size();
    if (count <= 2) {
        return 1.0;
    }
    double upperBound = 1.0;
    for (std::size_t gap = 0; gap + 1 < count; ++gap) {
        const double lowerNode = interpolationNodes[gap];
        const double upperNode = interpolationNodes[gap + 1];
        const double midpoint = lowerNode + 0.5 * (upperNode - lowerNode);
        const std::array<std::array<double, 2>, 2> subIntervals = {{
            {lowerNode, midpoint}, {midpoint, upperNode}}};
        for (const auto& subInterval : subIntervals) {
            std::vector<OutwardInterval> lebesgueCoefficients(
                count, IntervalPoint(0.0));
            for (std::size_t basisIndex = 0; basisIndex < count; ++basisIndex) {
                std::vector<OutwardInterval> basis(1, IntervalPoint(1.0));
                int sign = 1;
                for (std::size_t other = 0; other < count; ++other) {
                    if (other == basisIndex) {
                        continue;
                    }
                    const OutwardInterval denominator = IntervalSubtract(
                        IntervalPoint(interpolationNodes[basisIndex]),
                        IntervalPoint(interpolationNodes[other]));
                    const std::array<OutwardInterval, 2> factor = {
                        IntervalDivide(
                            IntervalSubtract(IntervalPoint(subInterval[0]),
                                             IntervalPoint(interpolationNodes[other])),
                            denominator),
                        IntervalDivide(
                            IntervalSubtract(IntervalPoint(subInterval[1]),
                                             IntervalPoint(interpolationNodes[other])),
                            denominator)};
                    basis = MultiplyBernsteinByLinear(basis, factor);

                    const bool numeratorNegative = other > gap;
                    const bool denominatorNegative = other > basisIndex;
                    if (numeratorNegative != denominatorNegative) {
                        sign = -sign;
                    }
                }
                for (std::size_t coefficient = 0; coefficient < count;
                     ++coefficient) {
                    const OutwardInterval signedCoefficient =
                        sign > 0 ? basis[coefficient]
                                 : IntervalNegate(basis[coefficient]);
                    lebesgueCoefficients[coefficient] = IntervalAdd(
                        lebesgueCoefficients[coefficient], signedCoefficient);
                }
            }
            for (const OutwardInterval& coefficient : lebesgueCoefficients) {
                upperBound = std::max(upperBound, coefficient.upper);
            }
        }
    }
    if (!std::isfinite(upperBound)) {
        throw std::runtime_error(
            "Cannot construct a finite conservative GLL Lebesgue bound.");
    }
    return RoundUp(upperBound);
}

std::vector<double> GllWeights(const std::vector<double>& nodes) {
    const int count = static_cast<int>(nodes.size());
    const int order = count - 1;
    std::vector<double> weights(nodes.size(), 0.0);
    for (int i = 0; i < count; ++i) {
        const double p = LegendreP(order, nodes[static_cast<std::size_t>(i)]);
        weights[static_cast<std::size_t>(i)] =
            2.0 / (static_cast<double>(order) *
                   static_cast<double>(order + 1) * p * p);
    }
    return weights;
}

std::vector<double> DerivativeMatrix(const std::vector<double>& nodes) {
    const int count = static_cast<int>(nodes.size());
    std::vector<double> bary(static_cast<std::size_t>(count), 1.0);
    for (int j = 0; j < count; ++j) {
        for (int k = 0; k < count; ++k) {
            if (j != k) {
                bary[static_cast<std::size_t>(j)] /=
                    nodes[static_cast<std::size_t>(j)] -
                    nodes[static_cast<std::size_t>(k)];
            }
        }
    }
    std::vector<double> matrix(static_cast<std::size_t>(count) *
                                   static_cast<std::size_t>(count),
                               0.0);
    for (int i = 0; i < count; ++i) {
        double sum = 0.0;
        for (int j = 0; j < count; ++j) {
            if (i == j) {
                continue;
            }
            const double value =
                bary[static_cast<std::size_t>(j)] /
                (bary[static_cast<std::size_t>(i)] *
                 (nodes[static_cast<std::size_t>(i)] -
                  nodes[static_cast<std::size_t>(j)]));
            matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(count) +
                   static_cast<std::size_t>(j)] = value;
            sum += value;
        }
        matrix[static_cast<std::size_t>(i) * static_cast<std::size_t>(count) +
               static_cast<std::size_t>(i)] = -sum;
    }
    return matrix;
}

std::vector<double> BarycentricWeights(const std::vector<double>& nodes) {
    std::vector<double> weights(nodes.size(), 1.0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            if (i != j) {
                weights[i] /= nodes[i] - nodes[j];
            }
        }
    }
    return weights;
}

void EvaluateNodalBasis(const std::vector<double>& nodes,
                        const std::vector<double>& barycentricWeights,
                        const std::vector<double>& nodalDerivativeMatrix,
                        double coordinate,
                        std::vector<double>* basis,
                        std::vector<double>* derivative) {
    if (basis == nullptr || derivative == nullptr ||
        nodes.size() != barycentricWeights.size() ||
        nodalDerivativeMatrix.size() != nodes.size() * nodes.size()) {
        throw std::runtime_error("Invalid nodal-basis evaluation workspace.");
    }
    basis->assign(nodes.size(), 0.0);
    derivative->assign(nodes.size(), 0.0);
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                                 std::max(1.0, std::abs(nodes[i]));
        if (std::abs(coordinate - nodes[i]) <= tolerance) {
            (*basis)[i] = 1.0;
            std::copy_n(nodalDerivativeMatrix.data() + i * nodes.size(),
                        nodes.size(), derivative->data());
            return;
        }
    }

    std::vector<double> scaled(nodes.size(), 0.0);
    double denominator = 0.0;
    double squaredDenominator = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double distance = coordinate - nodes[i];
        scaled[i] = barycentricWeights[i] / distance;
        denominator += scaled[i];
        squaredDenominator += scaled[i] / distance;
    }
    if (!std::isfinite(denominator) || denominator == 0.0) {
        throw std::runtime_error("Failed to evaluate a GLL barycentric basis.");
    }
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        const double value = scaled[i] / denominator;
        (*basis)[i] = value;
        (*derivative)[i] =
            value * (squaredDenominator / denominator -
                     1.0 / (coordinate - nodes[i]));
    }
}

double SquaredNorm3(const std::array<double, 3>& value) {
    return value[0] * value[0] + value[1] * value[1] + value[2] * value[2];
}

std::size_t NekLocalIndex(int nx, int ny, int i, int j, int k) {
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(ny) +
            static_cast<std::size_t>(j)) *
               static_cast<std::size_t>(nx) +
           static_cast<std::size_t>(i);
}

struct NekGeometryMetrics {
    std::vector<double> referenceWeights;
    std::vector<double> jacobians;
    std::vector<double> physicalWeights;
};

NekGeometryMetrics ComputeNekGeometryMetrics(vtkDataSet* dataSet,
                                             const Nek5000SpectralMetadata& meta,
                                             bool includeReferenceWeights,
                                             bool includeJacobians,
                                             bool includePhysicalWeights) {
    if (dataSet == nullptr ||
        static_cast<std::size_t>(dataSet->GetNumberOfPoints()) != meta.pointCount()) {
        throw std::runtime_error("Nek5000 geometry and GLL metadata point counts differ.");
    }
    const std::vector<double> nodesX = GllNodes(meta.nx);
    const std::vector<double> nodesY = GllNodes(meta.ny);
    const std::vector<double> nodesZ = meta.is3D() ? GllNodes(meta.nz)
                                                   : std::vector<double>{0.0};
    const std::vector<double> weightsX = GllWeights(nodesX);
    const std::vector<double> weightsY = GllWeights(nodesY);
    const std::vector<double> weightsZ = meta.is3D() ? GllWeights(nodesZ)
                                                     : std::vector<double>{1.0};
    const std::vector<double> dx = DerivativeMatrix(nodesX);
    const std::vector<double> dy = DerivativeMatrix(nodesY);
    const std::vector<double> dz = meta.is3D() ? DerivativeMatrix(nodesZ)
                                               : std::vector<double>();

    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(dataSet);
    vtkDataArray* coordinates =
        pointSet == nullptr || pointSet->GetPoints() == nullptr
            ? nullptr
            : pointSet->GetPoints()->GetData();
    if (coordinates == nullptr || coordinates->GetNumberOfComponents() != 3) {
        throw std::runtime_error(
            "Nek5000 geometry must expose a three-component point-coordinate array.");
    }
    vtkDoubleArray* doubleCoordinates = vtkDoubleArray::SafeDownCast(coordinates);
    const double* coordinateValues =
        doubleCoordinates == nullptr ? nullptr : doubleCoordinates->GetPointer(0);
    const auto coordinate = [coordinates, coordinateValues](std::size_t pointIndex,
                                                            int component) {
        return coordinateValues == nullptr
                   ? coordinates->GetComponent(static_cast<vtkIdType>(pointIndex), component)
                   : coordinateValues[3 * pointIndex + static_cast<std::size_t>(component)];
    };

    NekGeometryMetrics result;
    if (includeReferenceWeights) {
        result.referenceWeights.resize(meta.pointCount(), 0.0);
    }
    if (includeJacobians) {
        result.jacobians.resize(meta.pointCount(), 0.0);
    }
    if (includePhysicalWeights) {
        result.physicalWeights.resize(meta.pointCount(), 0.0);
    }
    const std::size_t nxyz = meta.pointsPerElement();
    for (std::size_t element = 0; element < meta.elementCount; ++element) {
        const std::size_t base = element * nxyz;
        for (int k = 0; k < meta.nz; ++k) {
            for (int j = 0; j < meta.ny; ++j) {
                for (int i = 0; i < meta.nx; ++i) {
                    std::array<double, 3> xr{};
                    std::array<double, 3> xs{};
                    std::array<double, 3> xt{};
                    for (int m = 0; m < meta.nx; ++m) {
                        const double d = dx[static_cast<std::size_t>(i) *
                                                static_cast<std::size_t>(meta.nx) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, m, j, k);
                        for (int component = 0; component < 3; ++component) {
                            xr[static_cast<std::size_t>(component)] +=
                                d * coordinate(index, component);
                        }
                    }
                    for (int m = 0; m < meta.ny; ++m) {
                        const double d = dy[static_cast<std::size_t>(j) *
                                                static_cast<std::size_t>(meta.ny) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, i, m, k);
                        for (int component = 0; component < 3; ++component) {
                            xs[static_cast<std::size_t>(component)] +=
                                d * coordinate(index, component);
                        }
                    }
                    if (meta.is3D()) {
                        for (int m = 0; m < meta.nz; ++m) {
                            const double d = dz[static_cast<std::size_t>(k) *
                                                    static_cast<std::size_t>(meta.nz) +
                                                static_cast<std::size_t>(m)];
                            const std::size_t index =
                                base + NekLocalIndex(meta.nx, meta.ny, i, j, m);
                            for (int component = 0; component < 3; ++component) {
                                xt[static_cast<std::size_t>(component)] +=
                                    d * coordinate(index, component);
                            }
                        }
                    }
                    double jacobian = 0.0;
                    if (meta.is3D()) {
                        jacobian =
                            xr[0] * (xs[1] * xt[2] - xt[1] * xs[2]) -
                            xs[0] * (xr[1] * xt[2] - xt[1] * xr[2]) +
                            xt[0] * (xr[1] * xs[2] - xs[1] * xr[2]);
                    } else {
                        const std::array<double, 3> cross = {
                            xr[1] * xs[2] - xr[2] * xs[1],
                            xr[2] * xs[0] - xr[0] * xs[2],
                            xr[0] * xs[1] - xr[1] * xs[0],
                        };
                        jacobian = std::sqrt(cross[0] * cross[0] + cross[1] * cross[1] +
                                             cross[2] * cross[2]);
                    }
                    if (!std::isfinite(jacobian) || std::abs(jacobian) <= 1.0e-30) {
                        throw std::runtime_error(
                            "Degenerate Nek5000 spectral-element Jacobian at element " +
                            std::to_string(element) + ".");
                    }
                    const double reference =
                        weightsX[static_cast<std::size_t>(i)] *
                        weightsY[static_cast<std::size_t>(j)] *
                        weightsZ[static_cast<std::size_t>(k)];
                    const std::size_t index =
                        base + NekLocalIndex(meta.nx, meta.ny, i, j, k);
                    if (includeReferenceWeights) {
                        result.referenceWeights[index] = reference;
                    }
                    if (includeJacobians) {
                        result.jacobians[index] = jacobian;
                    }
                    if (includePhysicalWeights) {
                        result.physicalWeights[index] =
                            reference * std::abs(jacobian);
                    }
                }
            }
        }
    }
    return result;
}

void AddDoublePointArray(vtkDataSet* dataSet,
                         const char* name,
                         const std::vector<double>& values) {
    vtkNew<vtkDoubleArray> array;
    array->SetName(name);
    array->SetNumberOfComponents(1);
    array->SetNumberOfTuples(static_cast<vtkIdType>(values.size()));
    std::copy(values.begin(), values.end(), array->GetPointer(0));
    dataSet->GetPointData()->AddArray(array);
}

void AttachNekQuadratureArrays(vtkDataSet* dataSet,
                               const NekGeometryMetrics& metrics) {
    if (!metrics.referenceWeights.empty()) {
        AddDoublePointArray(dataSet, kNek5000ReferenceWeightArray,
                            metrics.referenceWeights);
    }
    if (!metrics.jacobians.empty()) {
        AddDoublePointArray(dataSet, kNek5000JacobianDeterminantArray,
                            metrics.jacobians);
    }
    if (!metrics.physicalWeights.empty()) {
        AddDoublePointArray(dataSet, kNek5000QuadratureWeightArray,
                            metrics.physicalWeights);
    }
}

void AttachNekBasicMetadata(vtkDataSet* dataSet,
                            const Nek5000SpectralMetadata& meta) {
    vtkNew<vtkIntArray> dimensions;
    dimensions->SetName(kNek5000GllDimensionsField);
    dimensions->SetNumberOfComponents(3);
    dimensions->InsertNextTuple3(meta.nx, meta.ny, meta.nz);
    dataSet->GetFieldData()->AddArray(dimensions);

    vtkNew<vtkIntArray> count;
    count->SetName(kNek5000ElementCountField);
    count->InsertNextValue(static_cast<int>(meta.elementCount));
    dataSet->GetFieldData()->AddArray(count);

    vtkNew<vtkIntArray> ids;
    ids->SetName(kNek5000ElementIdsField);
    ids->SetNumberOfTuples(static_cast<vtkIdType>(meta.elementIds.size()));
    for (std::size_t i = 0; i < meta.elementIds.size(); ++i) {
        ids->SetValue(static_cast<vtkIdType>(i), meta.elementIds[i]);
    }
    dataSet->GetFieldData()->AddArray(ids);

    vtkNew<vtkStringArray> representation;
    representation->SetName(kNek5000RepresentationField);
    representation->InsertNextValue(
        "element-local GLL nodal values; interfaces intentionally unmerged");
    dataSet->GetFieldData()->AddArray(representation);
}

vtkSmartPointer<vtkDataSet> ReadNek5000(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    const NekDescriptor descriptor = ParseNekDescriptor(frame.path);
    if (frame.sourceIndex < 0 || frame.sourceIndex >= descriptor.numberOfTimeSteps) {
        throw std::runtime_error("Nek5000 frame index is out of range.");
    }
    const std::vector<NekFieldMeta> dataPieces =
        ReadNekPieceSet(descriptor, frame.sourceIndex);
    const auto piecesHaveCoordinates = [](const std::vector<NekFieldMeta>& pieces) {
        return !pieces.empty() &&
               std::all_of(pieces.begin(), pieces.end(),
                           [](const NekFieldMeta& piece) {
                               return piece.code.find('X') != std::string::npos;
                           });
    };
    const bool currentHasCoordinates = piecesHaveCoordinates(dataPieces);
    std::vector<NekFieldMeta> coordinatePieces;
    if (currentHasCoordinates) {
        coordinatePieces = dataPieces;
    } else if (frame.coordinateSourceIndex >= 0) {
        if (frame.coordinateSourceIndex >= descriptor.numberOfTimeSteps) {
            throw std::runtime_error(
                "Nek5000 coordinate frame index is out of range.");
        }
        coordinatePieces =
            ReadNekPieceSet(descriptor, frame.coordinateSourceIndex);
        if (!piecesHaveCoordinates(coordinatePieces)) {
            throw std::runtime_error(
                "Enumerated Nek5000 coordinate frame has no mesh block.");
        }
    } else {
        std::vector<int> candidates;
        candidates.reserve(static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        for (int candidate = frame.sourceIndex; candidate >= 0; --candidate) {
            candidates.push_back(candidate);
        }
        for (int candidate = frame.sourceIndex + 1;
             candidate < descriptor.numberOfTimeSteps; ++candidate) {
            candidates.push_back(candidate);
        }
        for (const int candidate : candidates) {
            std::vector<NekFieldMeta> pieces =
                ReadNekPieceSet(descriptor, candidate);
            if (piecesHaveCoordinates(pieces)) {
                coordinatePieces = std::move(pieces);
                break;
            }
        }
    }
    if (coordinatePieces.empty()) {
        throw std::runtime_error(
            "Nek5000 sequence does not contain a mesh-bearing field frame.");
    }
    ValidateNekPieceCompatible(dataPieces.front(), coordinatePieces.front(), false);

    const NekPieceLookup lookup = BuildNekPieceLookup(dataPieces);
    Nek5000SpectralMetadata metadata;
    metadata.nx = dataPieces.front().nx;
    metadata.ny = dataPieces.front().ny;
    metadata.nz = dataPieces.front().nz;
    metadata.elementCount = lookup.outputElementIds.size();
    metadata.elementIds = lookup.outputElementIds;
    if (metadata.pointCount() >
        static_cast<std::size_t>(std::numeric_limits<vtkIdType>::max())) {
        throw std::runtime_error("Nek5000 point count exceeds VTK limits.");
    }

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    ScatterNekCoordinateComponent(coordinatePieces, lookup, "x", &x);
    ScatterNekCoordinateComponent(coordinatePieces, lookup, "y", &y);
    if (metadata.is3D()) {
        ScatterNekCoordinateComponent(coordinatePieces, lookup, "z", &z);
    } else {
        z.assign(metadata.pointCount(), 0.0);
    }

    vtkNew<vtkDoubleArray> pointArray;
    pointArray->SetName("Points");
    pointArray->SetNumberOfComponents(3);
    pointArray->SetNumberOfTuples(static_cast<vtkIdType>(metadata.pointCount()));
    double* tuples = pointArray->GetPointer(0);
    for (std::size_t i = 0; i < metadata.pointCount(); ++i) {
        tuples[3 * i + 0] = x[i];
        tuples[3 * i + 1] = y[i];
        tuples[3 * i + 2] = z[i];
    }
    vtkNew<vtkPoints> points;
    points->SetData(pointArray);
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->SetPoints(points);

    if (options.readAllArrays) {
        const std::string code = dataPieces.front().code;
        if (code.find('U') != std::string::npos) {
            std::vector<double> u;
            std::vector<double> v;
            std::vector<double> w;
            ScatterNekComponent(dataPieces, lookup, "u", &u);
            ScatterNekComponent(dataPieces, lookup, "v", &v);
            if (metadata.is3D()) {
                ScatterNekComponent(dataPieces, lookup, "w", &w);
            } else {
                w.assign(metadata.pointCount(), 0.0);
            }
            vtkNew<vtkDoubleArray> velocity;
            velocity->SetName("Velocity");
            velocity->SetNumberOfComponents(3);
            velocity->SetNumberOfTuples(static_cast<vtkIdType>(metadata.pointCount()));
            for (std::size_t i = 0; i < metadata.pointCount(); ++i) {
                const double tuple[3] = {u[i], v[i], w[i]};
                velocity->SetTuple(static_cast<vtkIdType>(i), tuple);
            }
            grid->GetPointData()->AddArray(velocity);
        }

        const auto addScalar = [&](const std::string& label, const std::string& name) {
            std::vector<double> values;
            ScatterNekComponent(dataPieces, lookup, label, &values);
            vtkNew<vtkDoubleArray> array;
            array->SetName(name.c_str());
            array->SetNumberOfTuples(static_cast<vtkIdType>(values.size()));
            std::copy(values.begin(), values.end(), array->GetPointer(0));
            grid->GetPointData()->AddArray(array);
        };
        if (code.find('P') != std::string::npos) {
            addScalar("p", "Pressure");
        }
        if (code.find('T') != std::string::npos) {
            addScalar("t", "Temperature");
        }
        for (const std::string& label : dataPieces.front().labels) {
            if (!label.empty() && label.front() == 's') {
                std::string name = label;
                name.front() = 'S';
                addScalar(label, name);
            }
        }
    }

    vtkNew<vtkDoubleArray> time;
    time->SetName(kTimeValueArray);
    time->InsertNextValue(dataPieces.front().time);
    grid->GetFieldData()->AddArray(time);

    if (options.attachNek5000SpectralMetadata) {
        AttachNekBasicMetadata(grid, metadata);
        if (options.attachNek5000QuadratureArrays) {
            AttachNekQuadratureArrays(
                grid,
                ComputeNekGeometryMetrics(
                    grid, metadata,
                    options.attachNek5000ReferenceWeightArray,
                    options.attachNek5000JacobianDeterminantArray,
                    options.attachNek5000QuadratureWeightArray));
        }
    }
    return grid;
}

std::vector<DataSetFrame> ParsePvd(const std::filesystem::path& path,
                                   double fallbackPhysicalDt) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Cannot open PVD file: " + path.string());
    }
    std::ostringstream text;
    text << input.rdbuf();
    const std::regex tagPattern(R"(<\s*DataSet\b[^>]*>)", std::regex::icase);
    const std::regex filePattern(R"(file\s*=\s*(['"])(.*?)\1)", std::regex::icase);
    const std::regex timePattern(R"(timestep\s*=\s*(['"])(.*?)\1)", std::regex::icase);
    const std::string xml = text.str();
    std::vector<DataSetFrame> frames;
    for (std::sregex_iterator it(xml.begin(), xml.end(), tagPattern), end;
         it != end; ++it) {
        const std::string tag = it->str();
        std::smatch fileMatch;
        if (!std::regex_search(tag, fileMatch, filePattern)) {
            continue;
        }
        std::filesystem::path framePath(fileMatch[2].str());
        if (framePath.is_relative()) {
            framePath = path.parent_path() / framePath;
        }
        double time = static_cast<double>(frames.size()) * fallbackPhysicalDt;
        std::smatch timeMatch;
        if (std::regex_search(tag, timeMatch, timePattern)) {
            time = std::stod(timeMatch[2].str());
        }
        std::vector<DataSetFrame> expanded =
            EnumerateDataSetFrames(framePath, fallbackPhysicalDt);
        if (expanded.size() == 1) {
            expanded.front().time = time;
            if (expanded.front().kind != DataSetFrameKind::File) {
                expanded.front().readerTime = time;
            }
            frames.push_back(std::move(expanded.front()));
        } else {
            frames.insert(frames.end(), expanded.begin(), expanded.end());
        }
    }
    return frames;
}

long long LastNumericTag(const std::string& name) {
    const std::regex numbers(R"((\d+))");
    long long tag = -1;
    for (std::sregex_iterator it(name.begin(), name.end(), numbers), end;
         it != end; ++it) {
        tag = std::stoll((*it)[1].str());
    }
    return tag;
}

std::vector<DataSetFrame> EnumerateFolder(const std::filesystem::path& folder,
                                          double fallbackPhysicalDt) {
    struct TaggedFrame {
        long long tag = -1;
        std::string name;
        std::size_t localIndex = 0;
        DataSetFrame frame;
    };
    std::vector<TaggedFrame> tagged;
    std::unordered_set<std::string> nek5000Stems;
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
        if (entry.is_regular_file() &&
            LowerExtension(entry.path()) == ".nek5000") {
            nek5000Stems.insert(ToLower(entry.path().stem().string()));
        }
    }
    for (const auto& entry : std::filesystem::directory_iterator(folder)) {
        if (!entry.is_regular_file()) {
            continue;
        }
        const std::string extension = LowerExtension(entry.path());
        if (extension == ".nek" &&
            nek5000Stems.find(ToLower(entry.path().stem().string())) !=
                nek5000Stems.end()) {
            continue;
        }
        if (!IsDirectDataExtension(extension) && extension != ".nek5000" &&
            extension != ".nek" &&
            extension != ".case") {
            continue;
        }
        const std::vector<DataSetFrame> local =
            EnumerateDataSetFrames(entry.path(), fallbackPhysicalDt);
        for (std::size_t i = 0; i < local.size(); ++i) {
            tagged.push_back({LastNumericTag(entry.path().filename().string()),
                              entry.path().filename().string(), i, local[i]});
        }
    }
    std::sort(tagged.begin(), tagged.end(), [](const TaggedFrame& lhs,
                                               const TaggedFrame& rhs) {
        if (lhs.tag != rhs.tag) {
            return lhs.tag < rhs.tag;
        }
        if (lhs.name != rhs.name) {
            return lhs.name < rhs.name;
        }
        return lhs.localIndex < rhs.localIndex;
    });
    std::vector<DataSetFrame> frames;
    frames.reserve(tagged.size());
    for (std::size_t i = 0; i < tagged.size(); ++i) {
        DataSetFrame frame = std::move(tagged[i].frame);
        if (!frame.readerTime.has_value() && frame.kind == DataSetFrameKind::File) {
            frame.time = static_cast<double>(i) * fallbackPhysicalDt;
        }
        frames.push_back(std::move(frame));
    }
    return frames;
}

}  // namespace

struct Nek5000SpectralPointEvaluator::Impl {
    struct ElementBvhNode {
        std::array<double, 6> bounds{};
        double tolerance = 0.0;
        std::size_t begin = 0;
        std::size_t end = 0;
        std::size_t left = std::numeric_limits<std::size_t>::max();
        std::size_t right = std::numeric_limits<std::size_t>::max();

        bool IsLeaf() const noexcept {
            return left == std::numeric_limits<std::size_t>::max();
        }
    };

    vtkSmartPointer<vtkDataSet> dataSet;
    vtkSmartPointer<vtkDataArray> pointArray;
    Nek5000SpectralMetadata metadata;
    Nek5000PointEvaluationOptions options;
    int components = 0;
    std::array<std::vector<double>, 3> nodes;
    std::array<std::vector<double>, 3> barycentricWeights;
    std::array<std::vector<double>, 3> derivativeMatrices;
    vtkSmartPointer<vtkDataArray> coordinates;
    std::vector<std::array<double, 6>> elementBounds;
    std::vector<std::size_t> bvhElements;
    std::vector<ElementBvhNode> bvhNodes;
    std::size_t bvhRoot = std::numeric_limits<std::size_t>::max();
    vtkSmartPointer<vtkStaticCellLocator> cellLocator;
    std::size_t cellsPerElement = 0;
    std::shared_ptr<const Impl> geometryOwner;

    Impl(vtkDataSet* inputDataSet,
         vtkDataArray* inputArray,
         const Nek5000PointEvaluationOptions& inputOptions)
        : options(inputOptions) {
        if (inputDataSet == nullptr || inputArray == nullptr ||
            inputArray->GetNumberOfTuples() != inputDataSet->GetNumberOfPoints()) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation requires a complete point-data array.");
        }
        if (!std::isfinite(options.referenceTolerance) ||
            options.referenceTolerance < 0.0 ||
            !std::isfinite(options.relativePhysicalTolerance) ||
            options.relativePhysicalTolerance < 0.0 ||
            !std::isfinite(options.absolutePhysicalTolerance) ||
            options.absolutePhysicalTolerance < 0.0 ||
            options.maxNewtonIterations < 1) {
            throw std::runtime_error(
                "Invalid Nek5000 spectral point-evaluation options.");
        }
        dataSet = inputDataSet;
        pointArray = inputArray;
        metadata = GetNek5000SpectralMetadata(inputDataSet);
        components = inputArray->GetNumberOfComponents();
        if (components < 1) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation input has no components.");
        }

        nodes[0] = GllNodes(metadata.nx);
        nodes[1] = GllNodes(metadata.ny);
        nodes[2] = metadata.is3D() ? GllNodes(metadata.nz)
                                   : std::vector<double>{0.0};
        for (int axis = 0; axis < 3; ++axis) {
            barycentricWeights[static_cast<std::size_t>(axis)] =
                BarycentricWeights(nodes[static_cast<std::size_t>(axis)]);
            derivativeMatrices[static_cast<std::size_t>(axis)] =
                nodes[static_cast<std::size_t>(axis)].size() == 1
                    ? std::vector<double>{0.0}
                    : DerivativeMatrix(nodes[static_cast<std::size_t>(axis)]);
        }

        vtkPointSet* pointSet = vtkPointSet::SafeDownCast(inputDataSet);
        if (pointSet == nullptr || pointSet->GetPoints() == nullptr ||
            pointSet->GetPoints()->GetData() == nullptr) {
            throw std::runtime_error(
                "Nek5000 spectral interpolation requires explicit point coordinates.");
        }
        coordinates = pointSet->GetPoints()->GetData();

        const std::size_t nxyz = metadata.pointsPerElement();
        const std::array<std::vector<double>, 3> unitNodes = {
            [&]() {
                std::vector<double> values(nodes[0].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[0][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
            [&]() {
                std::vector<double> values(nodes[1].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[1][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
            [&]() {
                if (!metadata.is3D()) {
                    return std::vector<double>{0.0};
                }
                std::vector<double> values(nodes[2].size(), 0.0);
                for (std::size_t i = 0; i < values.size(); ++i) {
                    values[i] = 0.5 * (nodes[2][i] + 1.0);
                }
                values.front() = 0.0;
                values.back() = 1.0;
                return values;
            }(),
        };
        const double lebesgueProduct =
            ConservativeLebesgueBound(nodes[0]) *
            ConservativeLebesgueBound(nodes[1]) *
            (metadata.is3D() ? ConservativeLebesgueBound(nodes[2]) : 1.0);
        if (!std::isfinite(lebesgueProduct)) {
            throw std::runtime_error(
                "Cannot construct finite conservative Nek5000 element bounds.");
        }

        elementBounds.resize(metadata.elementCount);
        std::atomic<bool> invalidCoordinates{false};
        vtkSMPTools::For(
            static_cast<vtkIdType>(0),
            static_cast<vtkIdType>(metadata.elementCount),
            static_cast<vtkIdType>(32),
            [&](vtkIdType begin, vtkIdType end) {
                const std::size_t nxy =
                    static_cast<std::size_t>(metadata.nx) *
                    static_cast<std::size_t>(metadata.ny);
                std::vector<std::array<double, 3>> lowerFaces(nxy);
                std::vector<std::array<double, 3>> upperFaces(nxy);
                const auto interpolate = [](double lower, double upper,
                                            double parameter) {
                    if (parameter == 0.0) {
                        return lower;
                    }
                    if (parameter == 1.0) {
                        return upper;
                    }
                    return std::fma(parameter, upper - lower, lower);
                };

                for (vtkIdType elementId = begin; elementId < end; ++elementId) {
                    const std::size_t element =
                        static_cast<std::size_t>(elementId);
                    const std::size_t base = element * nxyz;
                    std::array<std::array<double, 3>, 8> corners{};
                    bool elementValid = true;
                    for (int cornerK = 0; cornerK < 2; ++cornerK) {
                        const int k = metadata.is3D() && cornerK != 0
                                          ? metadata.nz - 1
                                          : 0;
                        for (int cornerJ = 0; cornerJ < 2; ++cornerJ) {
                            const int j = cornerJ == 0 ? 0 : metadata.ny - 1;
                            for (int cornerI = 0; cornerI < 2; ++cornerI) {
                                const int i = cornerI == 0 ? 0 : metadata.nx - 1;
                                const std::size_t corner =
                                    static_cast<std::size_t>(cornerI +
                                                             2 * cornerJ +
                                                             4 * cornerK);
                                const vtkIdType pointIndex = static_cast<vtkIdType>(
                                    base + NekLocalIndex(metadata.nx, metadata.ny,
                                                         i, j, k));
                                coordinates->GetTuple(pointIndex,
                                                      corners[corner].data());
                                for (double coordinate : corners[corner]) {
                                    elementValid = elementValid &&
                                                   std::isfinite(coordinate);
                                }
                            }
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }

                    std::array<std::array<double, 3>, 8> localCorners{};
                    std::array<double, 3> cornerMinimum = {
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
                    std::array<double, 3> cornerMaximum = {
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
                    for (std::size_t corner = 0; corner < corners.size(); ++corner) {
                        for (int axis = 0; axis < 3; ++axis) {
                            const std::size_t component =
                                static_cast<std::size_t>(axis);
                            localCorners[corner][component] =
                                corners[corner][component] - corners[0][component];
                            cornerMinimum[component] =
                                std::min(cornerMinimum[component],
                                         corners[corner][component]);
                            cornerMaximum[component] =
                                std::max(cornerMaximum[component],
                                         corners[corner][component]);
                        }
                    }

                    for (int j = 0; j < metadata.ny; ++j) {
                        const double v = unitNodes[1][static_cast<std::size_t>(j)];
                        for (int i = 0; i < metadata.nx; ++i) {
                            const double u =
                                unitNodes[0][static_cast<std::size_t>(i)];
                            const std::size_t faceIndex =
                                static_cast<std::size_t>(j) *
                                    static_cast<std::size_t>(metadata.nx) +
                                static_cast<std::size_t>(i);
                            for (int axis = 0; axis < 3; ++axis) {
                                const std::size_t component =
                                    static_cast<std::size_t>(axis);
                                const double lowerLow = interpolate(
                                    localCorners[0][component],
                                    localCorners[1][component], u);
                                const double lowerHigh = interpolate(
                                    localCorners[2][component],
                                    localCorners[3][component], u);
                                const double upperLow = interpolate(
                                    localCorners[4][component],
                                    localCorners[5][component], u);
                                const double upperHigh = interpolate(
                                    localCorners[6][component],
                                    localCorners[7][component], u);
                                lowerFaces[faceIndex][component] =
                                    interpolate(lowerLow, lowerHigh, v);
                                upperFaces[faceIndex][component] =
                                    interpolate(upperLow, upperHigh, v);
                            }
                        }
                    }

                    std::array<double, 3> maximumResidual{};
                    std::array<double, 3> maximumLocalMagnitude{};
                    std::array<double, 3> nodeMinimum = {
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity(),
                        std::numeric_limits<double>::infinity()};
                    std::array<double, 3> nodeMaximum = {
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity(),
                        -std::numeric_limits<double>::infinity()};
                    for (int k = 0; k < metadata.nz && elementValid; ++k) {
                        const double w = metadata.is3D()
                                             ? unitNodes[2][static_cast<std::size_t>(k)]
                                             : 0.0;
                        for (int j = 0; j < metadata.ny && elementValid; ++j) {
                            for (int i = 0; i < metadata.nx; ++i) {
                                const std::size_t faceIndex =
                                    static_cast<std::size_t>(j) *
                                        static_cast<std::size_t>(metadata.nx) +
                                    static_cast<std::size_t>(i);
                                const std::size_t local = NekLocalIndex(
                                    metadata.nx, metadata.ny, i, j, k);
                                std::array<double, 3> point{};
                                coordinates->GetTuple(
                                    static_cast<vtkIdType>(base + local),
                                    point.data());
                                for (int axis = 0; axis < 3; ++axis) {
                                    const std::size_t component =
                                        static_cast<std::size_t>(axis);
                                    const double localPoint =
                                        point[component] - corners[0][component];
                                    const double trilinearPoint = interpolate(
                                        lowerFaces[faceIndex][component],
                                        upperFaces[faceIndex][component], w);
                                    if (!std::isfinite(point[component]) ||
                                        !std::isfinite(localPoint) ||
                                        !std::isfinite(trilinearPoint)) {
                                        elementValid = false;
                                        break;
                                    }
                                    nodeMinimum[component] =
                                        std::min(nodeMinimum[component],
                                                 point[component]);
                                    nodeMaximum[component] =
                                        std::max(nodeMaximum[component],
                                                 point[component]);
                                    maximumLocalMagnitude[component] = std::max(
                                        maximumLocalMagnitude[component],
                                        std::abs(localPoint));
                                    maximumResidual[component] = std::max(
                                        maximumResidual[component],
                                        std::abs(localPoint - trilinearPoint));
                                }
                            }
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }

                    std::array<double, 6> bounds{};
                    // The residual path contains three nested affine
                    // interpolations plus origin shifts. This factor is far
                    // above the corresponding floating-point gamma bound, so
                    // the value added below also encloses their roundoff.
                    constexpr double kResidualRoundoffFactor = 4096.0;
                    for (int axis = 0; axis < 3; ++axis) {
                        const std::size_t component =
                            static_cast<std::size_t>(axis);
                        const double residualPadding =
                            kResidualRoundoffFactor *
                            std::numeric_limits<double>::epsilon() *
                            maximumLocalMagnitude[component];
                        const double safeResidual = RoundUp(
                            maximumResidual[component] + residualPadding);
                        const double margin = RoundUp(
                            lebesgueProduct * safeResidual);
                        const double lower = RoundDown(
                            cornerMinimum[component] - margin);
                        const double upper = RoundUp(
                            cornerMaximum[component] + margin);
                        bounds[2 * component] =
                            std::min(lower, nodeMinimum[component]);
                        bounds[2 * component + 1] =
                            std::max(upper, nodeMaximum[component]);
                        if (!std::isfinite(bounds[2 * component]) ||
                            !std::isfinite(bounds[2 * component + 1])) {
                            elementValid = false;
                        }
                    }
                    if (!elementValid) {
                        invalidCoordinates.store(true, std::memory_order_relaxed);
                        continue;
                    }
                    elementBounds[element] = bounds;
                }
            });
        if (invalidCoordinates.load(std::memory_order_relaxed)) {
            throw std::runtime_error(
                "Nek5000 spectral coordinates must be finite and bounded.");
        }
        BuildElementBvh();

        cellsPerElement =
            static_cast<std::size_t>(metadata.nx - 1) *
            static_cast<std::size_t>(metadata.ny - 1) *
            static_cast<std::size_t>(metadata.is3D() ? metadata.nz - 1 : 1);
        if (cellsPerElement > 0 &&
            static_cast<std::size_t>(inputDataSet->GetNumberOfCells()) ==
                metadata.elementCount * cellsPerElement) {
            cellLocator = vtkSmartPointer<vtkStaticCellLocator>::New();
            cellLocator->SetDataSet(inputDataSet);
            cellLocator->BuildLocator();
        }
    }

    Impl(std::shared_ptr<const Impl> inputGeometry,
         vtkDataArray* inputArray)
        : geometryOwner(std::move(inputGeometry)) {
        if (geometryOwner == nullptr || inputArray == nullptr ||
            inputArray->GetNumberOfTuples() !=
                geometryOwner->dataSet->GetNumberOfPoints()) {
            throw std::runtime_error(
                "Rebound Nek5000 spectral field has an incompatible point count.");
        }
        pointArray = inputArray;
        components = inputArray->GetNumberOfComponents();
        options = geometryOwner->options;
        if (components < 1) {
            throw std::runtime_error(
                "Rebound Nek5000 spectral field has no components.");
        }
    }

    const Impl& Geometry() const {
        return geometryOwner == nullptr ? *this : *geometryOwner;
    }

    double PhysicalTolerance(std::size_t element) const {
        const auto& bounds = Geometry().elementBounds[element];
        const double dx = bounds[1] - bounds[0];
        const double dy = bounds[3] - bounds[2];
        const double dz = bounds[5] - bounds[4];
        return std::max(options.absolutePhysicalTolerance,
                        options.relativePhysicalTolerance *
                            std::sqrt(dx * dx + dy * dy + dz * dz));
    }

    bool BoundsContain(std::size_t element, const double point[3]) const {
        const auto& bounds = Geometry().elementBounds[element];
        const double tolerance = PhysicalTolerance(element);
        for (int axis = 0; axis < 3; ++axis) {
            if (point[axis] < bounds[2 * static_cast<std::size_t>(axis)] - tolerance ||
                point[axis] > bounds[2 * static_cast<std::size_t>(axis) + 1] + tolerance) {
                return false;
            }
        }
        return true;
    }

    std::size_t BuildElementBvhNode(std::size_t begin, std::size_t end) {
        ElementBvhNode node;
        node.bounds = {
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        std::array<double, 3> minimumCenter = {
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity(),
        };
        std::array<double, 3> maximumCenter = {
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
            -std::numeric_limits<double>::infinity(),
        };
        for (std::size_t offset = begin; offset < end; ++offset) {
            const std::size_t element = bvhElements[offset];
            const auto& bounds = elementBounds[element];
            node.tolerance = std::max(node.tolerance,
                                      PhysicalTolerance(element));
            for (int axis = 0; axis < 3; ++axis) {
                const std::size_t lo = 2 * static_cast<std::size_t>(axis);
                const std::size_t hi = lo + 1;
                node.bounds[lo] = std::min(node.bounds[lo], bounds[lo]);
                node.bounds[hi] = std::max(node.bounds[hi], bounds[hi]);
                const double center = 0.5 * (bounds[lo] + bounds[hi]);
                minimumCenter[static_cast<std::size_t>(axis)] =
                    std::min(minimumCenter[static_cast<std::size_t>(axis)], center);
                maximumCenter[static_cast<std::size_t>(axis)] =
                    std::max(maximumCenter[static_cast<std::size_t>(axis)], center);
            }
        }

        const std::size_t nodeIndex = bvhNodes.size();
        bvhNodes.push_back(node);
        constexpr std::size_t kLeafElementCount = 8;
        if (end - begin <= kLeafElementCount) {
            bvhNodes[nodeIndex].begin = begin;
            bvhNodes[nodeIndex].end = end;
            return nodeIndex;
        }

        int splitAxis = 0;
        double splitExtent = maximumCenter[0] - minimumCenter[0];
        for (int axis = 1; axis < 3; ++axis) {
            const double extent = maximumCenter[static_cast<std::size_t>(axis)] -
                                  minimumCenter[static_cast<std::size_t>(axis)];
            if (extent > splitExtent) {
                splitExtent = extent;
                splitAxis = axis;
            }
        }
        const std::size_t middle = begin + (end - begin) / 2;
        std::nth_element(
            bvhElements.begin() + static_cast<std::ptrdiff_t>(begin),
            bvhElements.begin() + static_cast<std::ptrdiff_t>(middle),
            bvhElements.begin() + static_cast<std::ptrdiff_t>(end),
            [&](std::size_t lhs, std::size_t rhs) {
                const std::size_t lo =
                    2 * static_cast<std::size_t>(splitAxis);
                const double lhsCenter =
                    0.5 * (elementBounds[lhs][lo] + elementBounds[lhs][lo + 1]);
                const double rhsCenter =
                    0.5 * (elementBounds[rhs][lo] + elementBounds[rhs][lo + 1]);
                return lhsCenter < rhsCenter ||
                       (lhsCenter == rhsCenter && lhs < rhs);
            });
        const std::size_t left = BuildElementBvhNode(begin, middle);
        const std::size_t right = BuildElementBvhNode(middle, end);
        bvhNodes[nodeIndex].left = left;
        bvhNodes[nodeIndex].right = right;
        return nodeIndex;
    }

    void BuildElementBvh() {
        bvhElements.resize(metadata.elementCount);
        std::iota(bvhElements.begin(), bvhElements.end(), 0);
        if (bvhElements.empty()) {
            return;
        }
        if (bvhElements.size() <=
            std::numeric_limits<std::size_t>::max() / 2) {
            bvhNodes.reserve(2 * bvhElements.size());
        }
        bvhRoot = BuildElementBvhNode(0, bvhElements.size());
    }

    bool BvhNodeContains(const ElementBvhNode& node,
                         const double point[3]) const {
        for (int axis = 0; axis < 3; ++axis) {
            const std::size_t lo = 2 * static_cast<std::size_t>(axis);
            if (point[axis] < node.bounds[lo] - node.tolerance ||
                point[axis] > node.bounds[lo + 1] + node.tolerance) {
                return false;
            }
        }
        return true;
    }

    void CollectBvhCandidates(const double point[3],
                              std::vector<std::size_t>* candidates) const {
        if (candidates == nullptr) {
            return;
        }
        candidates->clear();
        const Impl& geometry = Geometry();
        if (geometry.bvhRoot == std::numeric_limits<std::size_t>::max()) {
            return;
        }
        thread_local std::vector<std::size_t> stack;
        stack.clear();
        stack.push_back(geometry.bvhRoot);
        while (!stack.empty()) {
            const std::size_t nodeIndex = stack.back();
            stack.pop_back();
            const ElementBvhNode& node = geometry.bvhNodes[nodeIndex];
            if (!geometry.BvhNodeContains(node, point)) {
                continue;
            }
            if (node.IsLeaf()) {
                candidates->insert(
                    candidates->end(),
                    geometry.bvhElements.begin() +
                        static_cast<std::ptrdiff_t>(node.begin),
                    geometry.bvhElements.begin() +
                        static_cast<std::ptrdiff_t>(node.end));
                continue;
            }
            stack.push_back(node.right);
            stack.push_back(node.left);
        }
    }

    void EvaluateTensor(std::size_t element,
                        const std::array<double, 3>& reference,
                        std::array<double, 3>* mapped,
                        std::array<double, 3>* xr,
                        std::array<double, 3>* xs,
                        std::array<double, 3>* xt,
                        double* tuple) const {
        const Impl& geometry = Geometry();
        std::array<std::vector<double>, 3> basis;
        std::array<std::vector<double>, 3> derivative;
        for (int axis = 0; axis < 3; ++axis) {
            EvaluateNodalBasis(geometry.nodes[static_cast<std::size_t>(axis)],
                               geometry.barycentricWeights[static_cast<std::size_t>(axis)],
                               geometry.derivativeMatrices[static_cast<std::size_t>(axis)],
                               reference[static_cast<std::size_t>(axis)],
                               &basis[static_cast<std::size_t>(axis)],
                               &derivative[static_cast<std::size_t>(axis)]);
        }
        if (mapped != nullptr) {
            mapped->fill(0.0);
        }
        if (xr != nullptr) {
            xr->fill(0.0);
        }
        if (xs != nullptr) {
            xs->fill(0.0);
        }
        if (xt != nullptr) {
            xt->fill(0.0);
        }
        if (tuple != nullptr) {
            std::fill_n(tuple, components, 0.0);
        }

        const std::size_t base = element * geometry.metadata.pointsPerElement();
        for (int k = 0; k < geometry.metadata.nz; ++k) {
            for (int j = 0; j < geometry.metadata.ny; ++j) {
                for (int i = 0; i < geometry.metadata.nx; ++i) {
                    const double lx = basis[0][static_cast<std::size_t>(i)];
                    const double ly = basis[1][static_cast<std::size_t>(j)];
                    const double lz = basis[2][static_cast<std::size_t>(k)];
                    const double weight = lx * ly * lz;
                    const double weightR =
                        derivative[0][static_cast<std::size_t>(i)] * ly * lz;
                    const double weightS =
                        lx * derivative[1][static_cast<std::size_t>(j)] * lz;
                    const double weightT =
                        lx * ly * derivative[2][static_cast<std::size_t>(k)];
                    const std::size_t local =
                        NekLocalIndex(geometry.metadata.nx, geometry.metadata.ny,
                                      i, j, k);
                    const std::size_t pointIndex = base + local;
                    for (int axis = 0; axis < 3; ++axis) {
                        const double coordinate =
                            geometry.coordinates->GetComponent(
                                static_cast<vtkIdType>(pointIndex), axis);
                        if (mapped != nullptr) {
                            (*mapped)[static_cast<std::size_t>(axis)] +=
                                weight * coordinate;
                        }
                        if (xr != nullptr) {
                            (*xr)[static_cast<std::size_t>(axis)] +=
                                weightR * coordinate;
                        }
                        if (xs != nullptr) {
                            (*xs)[static_cast<std::size_t>(axis)] +=
                                weightS * coordinate;
                        }
                        if (xt != nullptr) {
                            (*xt)[static_cast<std::size_t>(axis)] +=
                                weightT * coordinate;
                        }
                    }
                    if (tuple != nullptr) {
                        for (int component = 0; component < components; ++component) {
                            tuple[component] +=
                                weight * pointArray->GetComponent(
                                             static_cast<vtkIdType>(pointIndex), component);
                        }
                    }
                }
            }
        }
    }

    static double Dot(const std::array<double, 3>& lhs,
                      const std::array<double, 3>& rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    }

    static std::array<double, 3> Cross(const std::array<double, 3>& lhs,
                                       const std::array<double, 3>& rhs) {
        return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
                lhs[2] * rhs[0] - lhs[0] * rhs[2],
                lhs[0] * rhs[1] - lhs[1] * rhs[0]};
    }

    bool NewtonFrom(std::size_t element,
                    const double point[3],
                    std::array<double, 3> reference,
                    std::array<double, 3>* result,
                    double* finalResidual) const {
        const Nek5000SpectralMetadata& geometryMetadata = Geometry().metadata;
        const double tolerance = PhysicalTolerance(element);
        const double toleranceSquared = tolerance * tolerance;
        double residualSquared = std::numeric_limits<double>::infinity();
        for (int iteration = 0; iteration < options.maxNewtonIterations;
             ++iteration) {
            std::array<double, 3> mapped{};
            std::array<double, 3> xr{};
            std::array<double, 3> xs{};
            std::array<double, 3> xt{};
            EvaluateTensor(element, reference, &mapped, &xr, &xs, &xt, nullptr);
            const std::array<double, 3> residual = {
                mapped[0] - point[0], mapped[1] - point[1], mapped[2] - point[2]};
            residualSquared = SquaredNorm3(residual);
            if (residualSquared <= toleranceSquared) {
                break;
            }

            std::array<double, 3> delta{};
            if (geometryMetadata.is3D()) {
                const double determinant = Dot(xr, Cross(xs, xt));
                if (!std::isfinite(determinant) ||
                    std::abs(determinant) <= 1.0e-30) {
                    return false;
                }
                delta[0] = Dot(residual, Cross(xs, xt)) / determinant;
                delta[1] = Dot(xr, Cross(residual, xt)) / determinant;
                delta[2] = Dot(xr, Cross(xs, residual)) / determinant;
            } else {
                const double g00 = Dot(xr, xr);
                const double g01 = Dot(xr, xs);
                const double g11 = Dot(xs, xs);
                const double determinant = g00 * g11 - g01 * g01;
                if (!std::isfinite(determinant) || determinant <= 1.0e-30) {
                    return false;
                }
                const double rhs0 = Dot(xr, residual);
                const double rhs1 = Dot(xs, residual);
                delta[0] = (g11 * rhs0 - g01 * rhs1) / determinant;
                delta[1] = (-g01 * rhs0 + g00 * rhs1) / determinant;
            }

            double scale = 1.0;
            bool accepted = false;
            for (int lineSearch = 0; lineSearch < 10; ++lineSearch) {
                std::array<double, 3> candidate = {
                    reference[0] - scale * delta[0],
                    reference[1] - scale * delta[1],
                    geometryMetadata.is3D() ? reference[2] - scale * delta[2] : 0.0};
                if (std::abs(candidate[0]) > 4.0 ||
                    std::abs(candidate[1]) > 4.0 ||
                    (geometryMetadata.is3D() && std::abs(candidate[2]) > 4.0)) {
                    scale *= 0.5;
                    continue;
                }
                std::array<double, 3> candidateMapped{};
                EvaluateTensor(element, candidate, &candidateMapped, nullptr,
                               nullptr, nullptr, nullptr);
                const std::array<double, 3> candidateResidual = {
                    candidateMapped[0] - point[0],
                    candidateMapped[1] - point[1],
                    candidateMapped[2] - point[2]};
                if (SquaredNorm3(candidateResidual) < residualSquared) {
                    reference = candidate;
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!accepted) {
                return false;
            }
        }

        std::array<double, 3> mapped{};
        EvaluateTensor(element, reference, &mapped, nullptr, nullptr, nullptr,
                       nullptr);
        const std::array<double, 3> residual = {
            mapped[0] - point[0], mapped[1] - point[1], mapped[2] - point[2]};
        residualSquared = SquaredNorm3(residual);
        if (finalResidual != nullptr) {
            *finalResidual = std::sqrt(residualSquared);
        }
        const double referenceLimit = 1.0 + options.referenceTolerance;
        const bool insideReference =
            std::abs(reference[0]) <= referenceLimit &&
            std::abs(reference[1]) <= referenceLimit &&
            (!geometryMetadata.is3D() || std::abs(reference[2]) <= referenceLimit);
        if (!insideReference ||
            residualSquared > PhysicalTolerance(element) * PhysicalTolerance(element)) {
            return false;
        }
        if (result != nullptr) {
            *result = reference;
        }
        return true;
    }

    bool Locate(std::size_t element,
                const double point[3],
                std::array<double, 3>* reference) const {
        const Impl& geometry = Geometry();
        if (!BoundsContain(element, point)) {
            return false;
        }

        const std::size_t nxyz = geometry.metadata.pointsPerElement();
        const std::size_t base = element * nxyz;
        std::size_t nearest = 0;
        double nearestDistance = std::numeric_limits<double>::infinity();
        for (std::size_t local = 0; local < nxyz; ++local) {
            const vtkIdType pointIndex = static_cast<vtkIdType>(base + local);
            const double dx = geometry.coordinates->GetComponent(pointIndex, 0) - point[0];
            const double dy = geometry.coordinates->GetComponent(pointIndex, 1) - point[1];
            const double dz = geometry.coordinates->GetComponent(pointIndex, 2) - point[2];
            const double distance = dx * dx + dy * dy + dz * dz;
            if (distance < nearestDistance) {
                nearestDistance = distance;
                nearest = local;
            }
        }
        const int nearestI = static_cast<int>(nearest %
                                               static_cast<std::size_t>(geometry.metadata.nx));
        const int nearestJ = static_cast<int>(
            (nearest / static_cast<std::size_t>(geometry.metadata.nx)) %
            static_cast<std::size_t>(geometry.metadata.ny));
        const int nearestK = static_cast<int>(
            nearest / (static_cast<std::size_t>(geometry.metadata.nx) *
                       static_cast<std::size_t>(geometry.metadata.ny)));

        std::array<std::array<double, 3>, 3> initial = {
            std::array<double, 3>{0.0, 0.0, 0.0},
            std::array<double, 3>{geometry.nodes[0][static_cast<std::size_t>(nearestI)],
                                  geometry.nodes[1][static_cast<std::size_t>(nearestJ)],
                                  geometry.metadata.is3D()
                                      ? geometry.nodes[2][static_cast<std::size_t>(nearestK)]
                                      : 0.0},
            std::array<double, 3>{0.0, 0.0, 0.0},
        };
        const auto& bounds = geometry.elementBounds[element];
        for (int axis = 0; axis < (geometry.metadata.is3D() ? 3 : 2); ++axis) {
            const double extent =
                bounds[2 * static_cast<std::size_t>(axis) + 1] -
                bounds[2 * static_cast<std::size_t>(axis)];
            if (extent > 0.0) {
                initial[2][static_cast<std::size_t>(axis)] =
                    std::max(-1.0, std::min(
                                       1.0,
                                       2.0 * (point[axis] -
                                              bounds[2 * static_cast<std::size_t>(axis)]) /
                                               extent -
                                           1.0));
            }
        }

        bool found = false;
        double bestResidual = std::numeric_limits<double>::infinity();
        std::array<double, 3> best{};
        for (const auto& guess : initial) {
            std::array<double, 3> candidate{};
            double residual = 0.0;
            if (NewtonFrom(element, point, guess, &candidate, &residual) &&
                residual < bestResidual) {
                bestResidual = residual;
                best = candidate;
                found = true;
            }
        }
        if (found && reference != nullptr) {
            *reference = best;
        }
        return found;
    }

    bool Evaluate(const double point[3],
                  double* tuple,
                  int* elementHint,
                  std::array<double, 3>* referencePoint) const {
        if (point == nullptr || tuple == nullptr || !std::isfinite(point[0]) ||
            !std::isfinite(point[1]) || !std::isfinite(point[2])) {
            return false;
        }
        const Impl& geometry = Geometry();
        const int hint = elementHint == nullptr ? -1 : *elementHint;
        const auto tryElement = [&](std::size_t element) {
            std::array<double, 3> reference{};
            if (!Locate(element, point, &reference)) {
                return false;
            }
            EvaluateTensor(element, reference, nullptr, nullptr, nullptr, nullptr,
                           tuple);
            if (elementHint != nullptr) {
                *elementHint = static_cast<int>(element);
            }
            if (referencePoint != nullptr) {
                *referencePoint = reference;
            }
            return true;
        };

        if (hint >= 0 &&
            static_cast<std::size_t>(hint) < geometry.metadata.elementCount &&
            tryElement(static_cast<std::size_t>(hint))) {
            return true;
        }
        int locatorElement = -1;
        if (geometry.cellLocator != nullptr) {
            thread_local vtkSmartPointer<vtkGenericCell> locatorCell =
                vtkSmartPointer<vtkGenericCell>::New();
            thread_local std::vector<double> locatorWeights;
            locatorWeights.resize(static_cast<std::size_t>(
                std::max<vtkIdType>(1, geometry.dataSet->GetMaxCellSize())));
            double query[3] = {point[0], point[1], point[2]};
            double parametricCoordinates[3] = {};
            int subId = 0;
            const vtkIdType cellId = geometry.cellLocator->FindCell(
                query, 0.0, locatorCell, subId, parametricCoordinates,
                locatorWeights.data());
            if (cellId >= 0) {
                locatorElement = static_cast<int>(
                    static_cast<std::size_t>(cellId) / geometry.cellsPerElement);
                if (locatorElement != hint && locatorElement >= 0 &&
                    static_cast<std::size_t>(locatorElement) <
                        geometry.metadata.elementCount &&
                    tryElement(static_cast<std::size_t>(locatorElement))) {
                    return true;
                }
            }
        }
        thread_local std::vector<std::size_t> bvhCandidates;
        geometry.CollectBvhCandidates(point, &bvhCandidates);
        for (const std::size_t element : bvhCandidates) {
            if (static_cast<int>(element) == hint ||
                static_cast<int>(element) == locatorElement) {
                continue;
            }
            if (tryElement(element)) {
                return true;
            }
        }
        if (elementHint != nullptr) {
            *elementHint = -1;
        }
        if (referencePoint != nullptr) {
            referencePoint->fill(0.0);
        }
        return false;
    }
};

Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    vtkDataSet* dataSet,
    vtkDataArray* pointArray,
    const Nek5000PointEvaluationOptions& options)
    : impl_(std::make_shared<Impl>(dataSet, pointArray, options)) {}

Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    std::shared_ptr<const Impl> impl)
    : impl_(std::move(impl)) {}

Nek5000SpectralPointEvaluator::~Nek5000SpectralPointEvaluator() = default;
Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    const Nek5000SpectralPointEvaluator&) noexcept = default;
Nek5000SpectralPointEvaluator& Nek5000SpectralPointEvaluator::operator=(
    const Nek5000SpectralPointEvaluator&) noexcept = default;
Nek5000SpectralPointEvaluator::Nek5000SpectralPointEvaluator(
    Nek5000SpectralPointEvaluator&&) noexcept = default;
Nek5000SpectralPointEvaluator& Nek5000SpectralPointEvaluator::operator=(
    Nek5000SpectralPointEvaluator&&) noexcept = default;

int Nek5000SpectralPointEvaluator::numberOfComponents() const noexcept {
    return impl_ == nullptr ? 0 : impl_->components;
}

std::array<double, 6> Nek5000SpectralPointEvaluator::bounds() const noexcept {
    if (impl_ == nullptr) {
        return {};
    }
    const Impl& geometry = impl_->Geometry();
    if (geometry.bvhRoot == std::numeric_limits<std::size_t>::max() ||
        geometry.bvhRoot >= geometry.bvhNodes.size()) {
        return {};
    }
    return geometry.bvhNodes[geometry.bvhRoot].bounds;
}

Nek5000SpectralPointEvaluator
Nek5000SpectralPointEvaluator::WithPointArray(vtkDataArray* pointArray) const {
    if (impl_ == nullptr) {
        throw std::runtime_error(
            "Cannot rebind an empty Nek5000 spectral point evaluator.");
    }
    std::shared_ptr<const Impl> geometry =
        impl_->geometryOwner == nullptr ? impl_ : impl_->geometryOwner;
    return Nek5000SpectralPointEvaluator(
        std::make_shared<Impl>(std::move(geometry), pointArray));
}

bool Nek5000SpectralPointEvaluator::EvaluateReference(
    std::size_t elementIndex,
    const double referencePoint[3],
    double* tuple) const {
    if (impl_ == nullptr || referencePoint == nullptr || tuple == nullptr) {
        return false;
    }
    const Impl& geometry = impl_->Geometry();
    if (elementIndex >= geometry.metadata.elementCount) {
        return false;
    }
    const double tolerance = geometry.options.referenceTolerance;
    std::array<double, 3> reference = {
        referencePoint[0], referencePoint[1],
        geometry.metadata.is3D() ? referencePoint[2] : 0.0};
    for (int axis = 0; axis < (geometry.metadata.is3D() ? 3 : 2); ++axis) {
        if (!std::isfinite(reference[static_cast<std::size_t>(axis)]) ||
            reference[static_cast<std::size_t>(axis)] < -1.0 - tolerance ||
            reference[static_cast<std::size_t>(axis)] > 1.0 + tolerance) {
            return false;
        }
        reference[static_cast<std::size_t>(axis)] = std::clamp(
            reference[static_cast<std::size_t>(axis)], -1.0, 1.0);
    }
    impl_->EvaluateTensor(elementIndex, reference, nullptr, nullptr, nullptr,
                          nullptr, tuple);
    return true;
}

bool Nek5000SpectralPointEvaluator::Evaluate(const double physicalPoint[3],
                                             double* tuple,
                                             int* elementHint) const {
    return impl_ != nullptr &&
           impl_->Evaluate(physicalPoint, tuple, elementHint, nullptr);
}

bool Nek5000SpectralPointEvaluator::EvaluateWithReference(
    const double physicalPoint[3],
    double* tuple,
    int* elementHint,
    std::array<double, 3>* referencePoint) const {
    return impl_ != nullptr && impl_->Evaluate(physicalPoint, tuple, elementHint,
                                                referencePoint);
}

bool Nek5000SpectralPointEvaluator::EvaluateGeometryReference(
    std::size_t elementIndex,
    const double referencePoint[3],
    std::array<double, 3>* physicalPoint,
    std::array<double, 3>* derivativeR,
    std::array<double, 3>* derivativeS,
    std::array<double, 3>* derivativeT) const {
    if (impl_ == nullptr || referencePoint == nullptr || physicalPoint == nullptr) {
        return false;
    }
    const Impl& geometry = impl_->Geometry();
    if (elementIndex >= geometry.metadata.elementCount) {
        return false;
    }
    std::array<double, 3> reference = {
        referencePoint[0], referencePoint[1],
        geometry.metadata.is3D() ? referencePoint[2] : 0.0};
    const double tolerance = geometry.options.referenceTolerance;
    for (int axis = 0; axis < (geometry.metadata.is3D() ? 3 : 2); ++axis) {
        if (!std::isfinite(reference[static_cast<std::size_t>(axis)]) ||
            reference[static_cast<std::size_t>(axis)] < -1.0 - tolerance ||
            reference[static_cast<std::size_t>(axis)] > 1.0 + tolerance) {
            return false;
        }
        reference[static_cast<std::size_t>(axis)] = std::clamp(
            reference[static_cast<std::size_t>(axis)], -1.0, 1.0);
    }
    impl_->EvaluateTensor(elementIndex, reference, physicalPoint, derivativeR,
                          derivativeS, derivativeT, nullptr);
    return true;
}

std::string ToString(DataSetFrameKind kind) {
    switch (kind) {
    case DataSetFrameKind::File:
        return "file";
    case DataSetFrameKind::Nek5000TimeStep:
        return "nek5000_time_step";
    case DataSetFrameKind::EnSightGoldTimeStep:
        return "ensight_gold_time_step";
    case DataSetFrameKind::OpenFoamTimeStep:
        return "openfoam_time_step";
    case DataSetFrameKind::MultiZoneManifestTimeStep:
        return "two_zone_manifest_time_step";
    }
    throw std::runtime_error("Unknown dataset frame kind.");
}

bool IsMultiZoneMeshManifest(const std::filesystem::path& sourcePath) {
    if (LowerExtension(sourcePath) != ".scpm") {
        return false;
    }
    try {
        ParseMultiZoneManifest(sourcePath);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

bool IsSupportedDataSource(const std::filesystem::path& sourcePath) {
    if (std::filesystem::is_directory(sourcePath)) {
        if (ResolveOpenFoamCase(sourcePath).has_value()) {
            return true;
        }
        for (const auto& entry :
             std::filesystem::directory_iterator(sourcePath)) {
            if (!entry.is_regular_file()) {
                continue;
            }
            const std::string extension = LowerExtension(entry.path());
            if (IsDirectDataExtension(extension) || extension == ".pvd" ||
                extension == ".nek5000" || extension == ".nek" ||
                extension == ".case") {
                return true;
            }
        }
        return false;
    }
    if (!std::filesystem::is_regular_file(
            streamcenterplus::vtkhdf::RealPath(sourcePath))) {
        return false;
    }
    const std::string extension = LowerExtension(sourcePath);
    if (extension == ".scpm") {
        return IsMultiZoneMeshManifest(sourcePath);
    }
    return IsDirectDataExtension(extension) || extension == ".pvd" ||
           (extension == ".nek5000" || extension == ".nek") ||
           extension == ".case" ||
           extension == ".foam";
}

std::vector<DataSetFrame> EnumerateDataSetFrames(
    const std::filesystem::path& sourcePath,
    double fallbackPhysicalDt) {
    if (!std::isfinite(fallbackPhysicalDt) || fallbackPhysicalDt <= 0.0) {
        throw std::runtime_error("Fallback physical timestep must be positive.");
    }
    if (const auto foam = ResolveOpenFoamCase(sourcePath)) {
        vtkSmartPointer<vtkPOpenFOAMReader> reader = CreateOpenFoamReader(*foam, true);
        const std::vector<double> times = OpenFoamTimeSteps(reader);
        std::vector<DataSetFrame> frames;
        frames.reserve(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            DataSetFrame frame;
            frame.path = *foam;
            frame.time = times[i];
            frame.kind = DataSetFrameKind::OpenFoamTimeStep;
            frame.readerTime = times[i];
            frame.sourceIndex = static_cast<int>(i);
            frames.push_back(std::move(frame));
        }
        return frames;
    }
    if (std::filesystem::is_directory(sourcePath)) {
        std::vector<DataSetFrame> frames =
            EnumerateFolder(sourcePath, fallbackPhysicalDt);
        if (frames.empty()) {
            throw std::runtime_error("No supported dataset frames were found in: " +
                                     sourcePath.string());
        }
        return frames;
    }
    const std::filesystem::path realPath =
        streamcenterplus::vtkhdf::RealPath(sourcePath);
    if (!std::filesystem::is_regular_file(realPath)) {
        throw std::runtime_error("Dataset source does not exist: " +
                                 sourcePath.string());
    }
    const std::string extension = LowerExtension(sourcePath);
    if (extension == ".scpm") {
        return EnumerateMultiZoneManifestFrames(sourcePath, fallbackPhysicalDt);
    }
    if (extension == ".pvd") {
        return ParsePvd(sourcePath, fallbackPhysicalDt);
    }
    if (extension == ".nek5000" || extension == ".nek") {
        const NekDescriptor descriptor = ParseNekDescriptor(sourcePath);
        std::vector<DataSetFrame> frames;
        frames.reserve(static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        std::vector<bool> meshBearingFrames;
        meshBearingFrames.reserve(
            static_cast<std::size_t>(descriptor.numberOfTimeSteps));
        for (int i = 0; i < descriptor.numberOfTimeSteps; ++i) {
            const NekFieldMeta meta =
                ReadNekHeader(FormatNekFieldPath(descriptor, i, 0));
            DataSetFrame frame;
            frame.path = descriptor.path;
            frame.time = meta.time;
            frame.kind = DataSetFrameKind::Nek5000TimeStep;
            frame.readerTime = meta.time;
            frame.sourceIndex = i;
            frames.push_back(std::move(frame));
            meshBearingFrames.push_back(meta.code.find('X') != std::string::npos);
        }

        int previousMeshFrame = -1;
        for (int i = 0; i < descriptor.numberOfTimeSteps; ++i) {
            if (meshBearingFrames[static_cast<std::size_t>(i)]) {
                previousMeshFrame = i;
            }
            frames[static_cast<std::size_t>(i)].coordinateSourceIndex =
                previousMeshFrame;
        }
        int nextMeshFrame = -1;
        for (int i = descriptor.numberOfTimeSteps - 1; i >= 0; --i) {
            if (meshBearingFrames[static_cast<std::size_t>(i)]) {
                nextMeshFrame = i;
            }
            DataSetFrame& frame = frames[static_cast<std::size_t>(i)];
            if (frame.coordinateSourceIndex < 0) {
                frame.coordinateSourceIndex = nextMeshFrame;
            }
            if (frame.coordinateSourceIndex >= 0) {
                frame.coordinatePath = FormatNekFieldPath(
                    descriptor, frame.coordinateSourceIndex, 0);
            }
        }
        if (nextMeshFrame < 0) {
            throw std::runtime_error(
                "Nek5000 sequence does not contain a mesh-bearing field frame.");
        }
        return frames;
    }
    if (extension == ".case") {
        vtkSmartPointer<vtkGenericEnSightReader> reader =
            CreateEnSightGoldReader(sourcePath, true);
        const std::vector<double> times = EnSightGoldTimeSteps(reader);
        std::vector<DataSetFrame> frames;
        frames.reserve(times.size());
        for (std::size_t i = 0; i < times.size(); ++i) {
            DataSetFrame frame;
            frame.path = std::filesystem::absolute(sourcePath);
            frame.time = times[i];
            frame.kind = DataSetFrameKind::EnSightGoldTimeStep;
            frame.readerTime = times[i];
            frame.sourceIndex = static_cast<int>(i);
            frames.push_back(std::move(frame));
        }
        return frames;
    }
    if (!IsDirectDataExtension(extension)) {
        throw std::runtime_error("Unsupported dataset source: " + sourcePath.string());
    }
    const std::vector<std::filesystem::path> expanded =
        IsTecplotExtension(extension)
            ? std::vector<std::filesystem::path>{sourcePath}
            : streamcenterplus::vtkhdf::ExpandTimeSteps(sourcePath);
    const bool isVtkHdf = streamcenterplus::vtkhdf::IsVtkHdfPath(sourcePath);
    std::vector<double> nativeTimes;
    if (isVtkHdf) {
        nativeTimes = streamcenterplus::vtkhdf::TimeValues(sourcePath);
        if (nativeTimes.size() != expanded.size() ||
            !std::all_of(nativeTimes.cbegin(), nativeTimes.cend(), IsFiniteTime)) {
            throw std::runtime_error(
                "VTKHDF transient time metadata is missing, non-finite, or inconsistent "
                "with its frame count: " + sourcePath.string());
        }
    }
    std::vector<DataSetFrame> frames;
    frames.reserve(expanded.size());
    for (std::size_t i = 0; i < expanded.size(); ++i) {
        DataSetFrame frame;
        frame.path = expanded[i];
        frame.time = isVtkHdf ? nativeTimes[i] : static_cast<double>(i) * fallbackPhysicalDt;
        frame.kind = DataSetFrameKind::File;
        // File frames normally have no reader-time value, and folder
        // enumeration uses that fact to assign its fallback cadence. Preserve
        // a VTKHDF file's native time explicitly so a directory source cannot
        // overwrite non-uniform transient times.
        if (isVtkHdf) {
            frame.readerTime = nativeTimes[i];
        }
        frame.sourceIndex = static_cast<int>(i);
        frames.push_back(std::move(frame));
    }
    return frames;
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const std::filesystem::path& sourcePath,
                                        const ReadOptions& options) {
    std::vector<DataSetFrame> frames = EnumerateDataSetFrames(sourcePath);
    if (frames.empty()) {
        throw std::runtime_error("Dataset source contains no frames: " +
                                 sourcePath.string());
    }
    return ReadDataSet(frames.front(), options);
}

vtkSmartPointer<vtkDataSet> ReadDataSet(const DataSetFrame& frame,
                                        const ReadOptions& options) {
    switch (frame.kind) {
    case DataSetFrameKind::File:
        return ReadOrdinaryVtk(frame.path, options);
    case DataSetFrameKind::Nek5000TimeStep:
        return ReadNek5000(frame, options);
    case DataSetFrameKind::EnSightGoldTimeStep:
        return ReadEnSightGold(frame, options);
    case DataSetFrameKind::OpenFoamTimeStep:
        return ReadOpenFoam(frame, options);
    case DataSetFrameKind::MultiZoneManifestTimeStep:
        return ReadMultiZoneManifestFrame(frame, options);
    }
    throw std::runtime_error("Unsupported dataset frame kind.");
}

bool HasNek5000SpectralMetadata(vtkDataSet* dataSet) {
    if (dataSet == nullptr || dataSet->GetFieldData() == nullptr) {
        return false;
    }
    vtkDataArray* dimensions =
        dataSet->GetFieldData()->GetArray(kNek5000GllDimensionsField);
    vtkDataArray* count =
        dataSet->GetFieldData()->GetArray(kNek5000ElementCountField);
    return dimensions != nullptr && dimensions->GetNumberOfComponents() == 3 &&
           dimensions->GetNumberOfTuples() == 1 && count != nullptr &&
           count->GetNumberOfTuples() == 1;
}

Nek5000SpectralMetadata GetNek5000SpectralMetadata(vtkDataSet* dataSet) {
    if (!HasNek5000SpectralMetadata(dataSet)) {
        throw std::runtime_error("Dataset has no Streamcenter+ Nek5000 spectral metadata.");
    }
    Nek5000SpectralMetadata metadata;
    vtkDataArray* dimensions =
        dataSet->GetFieldData()->GetArray(kNek5000GllDimensionsField);
    metadata.nx = static_cast<int>(dimensions->GetComponent(0, 0));
    metadata.ny = static_cast<int>(dimensions->GetComponent(0, 1));
    metadata.nz = static_cast<int>(dimensions->GetComponent(0, 2));
    metadata.elementCount = static_cast<std::size_t>(
        dataSet->GetFieldData()->GetArray(kNek5000ElementCountField)
            ->GetComponent(0, 0));
    if (vtkDataArray* ids =
            dataSet->GetFieldData()->GetArray(kNek5000ElementIdsField)) {
        metadata.elementIds.resize(static_cast<std::size_t>(ids->GetNumberOfTuples()));
        for (vtkIdType i = 0; i < ids->GetNumberOfTuples(); ++i) {
            metadata.elementIds[static_cast<std::size_t>(i)] =
                static_cast<int>(ids->GetComponent(i, 0));
        }
    }
    if (metadata.nx < 2 || metadata.ny < 2 || metadata.nz < 1 ||
        metadata.elementCount == 0 || metadata.pointCount() !=
                                            static_cast<std::size_t>(dataSet->GetNumberOfPoints())) {
        throw std::runtime_error("Dataset has inconsistent Nek5000 spectral metadata.");
    }
    return metadata;
}

namespace {

constexpr const char* kSemAwareDisplayRepresentation =
    "SEM-aware uniform-reference sampled display proxy; approximate VTK "
    "visualization/interoperability only";
constexpr const char* kSemAwareOperationDisplayRepresentation =
    "SEM-aware evaluator-sampled operation display approximation; canonical "
    "Nek5000 GLL source remains immutable and unmerged";

vtkSmartPointer<vtkPoints> NewDoublePoints() {
    vtkSmartPointer<vtkDoubleArray> coordinates =
        vtkSmartPointer<vtkDoubleArray>::New();
    coordinates->SetNumberOfComponents(3);
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetData(coordinates);
    return points;
}

void MarkNek5000DisplayProxy(vtkDataSet* dataSet,
                             const char* representation,
                             const char* operation = nullptr,
                             double maximumReferenceSpacing = 0.0,
                             const Nek5000SpectralMetadata* metadata = nullptr) {
    if (dataSet == nullptr) {
        throw std::runtime_error("Cannot mark a null Nek5000 display proxy.");
    }
    if (dataSet->GetFieldData() == nullptr) {
        vtkNew<vtkFieldData> createdFieldData;
        dataSet->SetFieldData(createdFieldData);
    }
    vtkFieldData* fieldData = dataSet->GetFieldData();
    fieldData->RemoveArray(kNek5000DisplayProxyField);
    fieldData->RemoveArray(kNek5000DisplayRepresentationField);
    fieldData->RemoveArray(kNek5000DisplayOperationField);
    fieldData->RemoveArray(kNek5000DisplayReferenceSpacingField);
    fieldData->RemoveArray(kNek5000DisplaySourceElementIdsField);
    vtkNew<vtkStringArray> marker;
    marker->SetName(kNek5000DisplayProxyField);
    marker->InsertNextValue("1");
    fieldData->AddArray(marker);

    vtkNew<vtkStringArray> description;
    description->SetName(kNek5000DisplayRepresentationField);
    description->InsertNextValue(representation);
    fieldData->AddArray(description);

    if (operation != nullptr) {
        vtkNew<vtkStringArray> operationName;
        operationName->SetName(kNek5000DisplayOperationField);
        operationName->InsertNextValue(operation);
        fieldData->AddArray(operationName);
    }
    if (maximumReferenceSpacing > 0.0 &&
        std::isfinite(maximumReferenceSpacing)) {
        vtkNew<vtkDoubleArray> spacing;
        spacing->SetName(kNek5000DisplayReferenceSpacingField);
        spacing->SetNumberOfComponents(1);
        spacing->InsertNextValue(maximumReferenceSpacing);
        fieldData->AddArray(spacing);
    }
    if (metadata != nullptr) {
        vtkNew<vtkIntArray> sourceElements;
        sourceElements->SetName(kNek5000DisplaySourceElementIdsField);
        sourceElements->SetNumberOfComponents(1);
        sourceElements->SetNumberOfTuples(
            static_cast<vtkIdType>(metadata->elementIds.size()));
        for (std::size_t index = 0; index < metadata->elementIds.size(); ++index) {
            sourceElements->SetValue(static_cast<vtkIdType>(index),
                                     metadata->elementIds[index]);
        }
        fieldData->AddArray(sourceElements);
    }
}

bool IsNek5000DisplayDiagnosticArray(const char* name) {
    if (name == nullptr) {
        return true;
    }
    return std::strcmp(name, kNek5000ReferenceWeightArray) == 0 ||
           std::strcmp(name, kNek5000JacobianDeterminantArray) == 0 ||
           std::strcmp(name, kNek5000QuadratureWeightArray) == 0 ||
           std::strcmp(name, kNek5000SpectralElementIdArray) == 0 ||
           std::strcmp(name, kNek5000DisplayElementIdArray) == 0 ||
           std::strcmp(name, kNek5000DisplayElementIndexArray) == 0 ||
           std::strcmp(name, kNek5000DisplayLocalGllIndexArray) == 0 ||
           std::strcmp(name, kNek5000DisplayReferenceCoordinatesArray) == 0 ||
           std::strcmp(name, kNek5000DisplayFaceIdArray) == 0;
}

int SampleCountForGllAxis(int gllCount,
                           const Nek5000SpectralDisplayOptions& options) {
    if (gllCount < 1 || options.samplesPerGllInterval < 1 ||
        options.maxSamplesPerAxis < 2) {
        throw std::runtime_error(
            "Invalid Nek5000 spectral display sampling options.");
    }
    if (gllCount == 1) {
        return 1;
    }
    const long long desired =
        static_cast<long long>(gllCount - 1) *
            static_cast<long long>(options.samplesPerGllInterval) +
        1LL;
    return static_cast<int>(std::min<long long>(
        desired, static_cast<long long>(options.maxSamplesPerAxis)));
}

double UniformReferenceCoordinate(int index, int count) {
    if (count <= 1) {
        return 0.0;
    }
    if (index == 0) {
        return -1.0;
    }
    if (index == count - 1) {
        return 1.0;
    }
    return -1.0 + 2.0 * static_cast<double>(index) /
                      static_cast<double>(count - 1);
}

std::size_t SampleLocalIndex(int nx, int ny, int i, int j, int k) {
    return NekLocalIndex(nx, ny, i, j, k);
}

struct NekDisplayCellProvenance {
    vtkSmartPointer<vtkIntArray> elementIds = vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> elementIndices =
        vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> faceIds = vtkSmartPointer<vtkIntArray>::New();

    NekDisplayCellProvenance() {
        elementIds->SetName(kNek5000DisplayElementIdArray);
        elementIndices->SetName(kNek5000DisplayElementIndexArray);
        faceIds->SetName(kNek5000DisplayFaceIdArray);
    }

    void Add(const Nek5000SpectralMetadata& metadata,
             std::size_t elementIndex,
             int faceId) {
        elementIds->InsertNextValue(metadata.elementIds.at(elementIndex));
        elementIndices->InsertNextValue(static_cast<int>(elementIndex));
        faceIds->InsertNextValue(faceId);
    }

    void Attach(vtkDataSet* dataSet) const {
        dataSet->GetCellData()->AddArray(elementIds);
        dataSet->GetCellData()->AddArray(elementIndices);
        dataSet->GetCellData()->AddArray(faceIds);
    }
};

class NekDisplayPointWriter {
public:
    NekDisplayPointWriter(vtkPointSet* output,
                          vtkDataSet* canonicalDataSet,
                          const Nek5000SpectralMetadata& metadata,
                          const Nek5000SpectralPointEvaluator& geometryEvaluator,
                          bool interpolateFields)
        : output_(output),
          metadata_(metadata),
          geometryEvaluator_(geometryEvaluator),
          gllNodes_{GllNodes(metadata.nx), GllNodes(metadata.ny),
                    metadata.is3D() ? GllNodes(metadata.nz)
                                    : std::vector<double>{0.0}} {
        if (output_ == nullptr || output_->GetPoints() == nullptr ||
            canonicalDataSet == nullptr) {
            throw std::runtime_error(
                "Cannot create a Nek5000 display proxy without points.");
        }
        elementIds_ = vtkSmartPointer<vtkIntArray>::New();
        elementIds_->SetName(kNek5000DisplayElementIdArray);
        elementIndices_ = vtkSmartPointer<vtkIntArray>::New();
        elementIndices_->SetName(kNek5000DisplayElementIndexArray);
        localGllIndices_ = vtkSmartPointer<vtkIntArray>::New();
        localGllIndices_->SetName(kNek5000DisplayLocalGllIndexArray);
        referenceCoordinates_ = vtkSmartPointer<vtkDoubleArray>::New();
        referenceCoordinates_->SetName(kNek5000DisplayReferenceCoordinatesArray);
        referenceCoordinates_->SetNumberOfComponents(3);
        faceIds_ = vtkSmartPointer<vtkIntArray>::New();
        faceIds_->SetName(kNek5000DisplayFaceIdArray);

        if (interpolateFields) {
            vtkPointData* pointData = canonicalDataSet->GetPointData();
            for (int index = 0; pointData != nullptr &&
                                index < pointData->GetNumberOfArrays();
                 ++index) {
                vtkDataArray* input = pointData->GetArray(index);
                const char* name = input == nullptr ? nullptr : input->GetName();
                if (input == nullptr || name == nullptr ||
                    input->GetNumberOfComponents() < 1 ||
                    input->GetNumberOfTuples() !=
                        canonicalDataSet->GetNumberOfPoints() ||
                    IsNek5000DisplayDiagnosticArray(name)) {
                    continue;
                }
                vtkSmartPointer<vtkDoubleArray> values =
                    vtkSmartPointer<vtkDoubleArray>::New();
                values->SetName(name);
                values->SetNumberOfComponents(input->GetNumberOfComponents());
                output_->GetPointData()->AddArray(values);
                fields_.emplace_back(geometryEvaluator_.WithPointArray(input),
                                     std::move(values));
            }
        }
        output_->GetPointData()->AddArray(elementIds_);
        output_->GetPointData()->AddArray(elementIndices_);
        output_->GetPointData()->AddArray(localGllIndices_);
        output_->GetPointData()->AddArray(referenceCoordinates_);
        output_->GetPointData()->AddArray(faceIds_);
    }

    vtkIdType Add(std::size_t elementIndex,
                  const std::array<double, 3>& referencePoint,
                  int knownLocalGllIndex = -1,
                  int faceId = -1) {
        std::array<double, 3> physicalPoint{};
        if (!geometryEvaluator_.EvaluateReference(elementIndex,
                                                  referencePoint.data(),
                                                  physicalPoint.data())) {
            throw std::runtime_error(
                "Failed to evaluate a Nek5000 display-proxy coordinate.");
        }
        const vtkIdType id = output_->GetPoints()->InsertNextPoint(
            physicalPoint.data());
        for (Field& field : fields_) {
            std::vector<double> tuple(
                static_cast<std::size_t>(field.evaluator.numberOfComponents()));
            if (!field.evaluator.EvaluateReference(elementIndex,
                                                    referencePoint.data(),
                                                    tuple.data())) {
                throw std::runtime_error(
                    "Failed to evaluate a Nek5000 display-proxy field.");
            }
            field.values->InsertNextTuple(tuple.data());
        }
        elementIds_->InsertNextValue(metadata_.elementIds.at(elementIndex));
        elementIndices_->InsertNextValue(static_cast<int>(elementIndex));
        localGllIndices_->InsertNextValue(
            knownLocalGllIndex >= 0 ? knownLocalGllIndex
                                    : LocalGllIndex(referencePoint));
        referenceCoordinates_->InsertNextTuple(referencePoint.data());
        faceIds_->InsertNextValue(faceId);
        return id;
    }

private:
    struct Field {
        Field(Nek5000SpectralPointEvaluator inputEvaluator,
              vtkSmartPointer<vtkDoubleArray> inputValues)
            : evaluator(std::move(inputEvaluator)),
              values(std::move(inputValues)) {}

        Nek5000SpectralPointEvaluator evaluator;
        vtkSmartPointer<vtkDoubleArray> values;
    };

    int NodeIndex(int axis, double coordinate) const {
        const std::vector<double>& nodes = gllNodes_.at(
            static_cast<std::size_t>(axis));
        for (std::size_t index = 0; index < nodes.size(); ++index) {
            if (std::abs(nodes[index] - coordinate) <=
                64.0 * std::numeric_limits<double>::epsilon()) {
                return static_cast<int>(index);
            }
        }
        return -1;
    }

    int LocalGllIndex(const std::array<double, 3>& referencePoint) const {
        const int i = NodeIndex(0, referencePoint[0]);
        const int j = NodeIndex(1, referencePoint[1]);
        const int k = metadata_.is3D() ? NodeIndex(2, referencePoint[2]) : 0;
        if (i < 0 || j < 0 || k < 0) {
            return -1;
        }
        return static_cast<int>(NekLocalIndex(metadata_.nx, metadata_.ny, i, j, k));
    }

    vtkPointSet* output_ = nullptr;
    const Nek5000SpectralMetadata& metadata_;
    Nek5000SpectralPointEvaluator geometryEvaluator_;
    std::array<std::vector<double>, 3> gllNodes_;
    std::vector<Field> fields_;
    vtkSmartPointer<vtkIntArray> elementIds_;
    vtkSmartPointer<vtkIntArray> elementIndices_;
    vtkSmartPointer<vtkIntArray> localGllIndices_;
    vtkSmartPointer<vtkDoubleArray> referenceCoordinates_;
    vtkSmartPointer<vtkIntArray> faceIds_;
};

int NekFaceCount(const Nek5000SpectralMetadata& metadata) {
    return metadata.is3D() ? 6 : 4;
}

std::array<double, 3> FaceReference(const Nek5000SpectralMetadata& metadata,
                                    int faceId,
                                    double first,
                                    double second = 0.0) {
    if (!metadata.is3D()) {
        switch (faceId) {
        case 0:
            return {-1.0, first, 0.0};
        case 1:
            return {1.0, first, 0.0};
        case 2:
            return {first, -1.0, 0.0};
        case 3:
            return {first, 1.0, 0.0};
        default:
            throw std::runtime_error("Invalid 2D Nek5000 face id.");
        }
    }
    switch (faceId) {
    case 0:
        return {-1.0, first, second};
    case 1:
        return {1.0, first, second};
    case 2:
        return {first, -1.0, second};
    case 3:
        return {first, 1.0, second};
    case 4:
        return {first, second, -1.0};
    case 5:
        return {first, second, 1.0};
    default:
        throw std::runtime_error("Invalid 3D Nek5000 face id.");
    }
}

std::array<int, 2> FaceSampleDimensions(const Nek5000SpectralMetadata& metadata,
                                        int faceId,
                                        int samplesX,
                                        int samplesY,
                                        int samplesZ) {
    if (!metadata.is3D()) {
        return (faceId == 0 || faceId == 1)
                   ? std::array<int, 2>{samplesY, 1}
                   : std::array<int, 2>{samplesX, 1};
    }
    switch (faceId) {
    case 0:
    case 1:
        return {samplesY, samplesZ};
    case 2:
    case 3:
        return {samplesX, samplesZ};
    case 4:
    case 5:
        return {samplesX, samplesY};
    default:
        throw std::runtime_error("Invalid Nek5000 face id.");
    }
}

std::vector<std::size_t> FaceLocalIndices(
    const Nek5000SpectralMetadata& metadata,
    int faceId) {
    std::vector<std::size_t> indices;
    const auto add = [&](int i, int j, int k) {
        indices.push_back(NekLocalIndex(metadata.nx, metadata.ny, i, j, k));
    };
    if (!metadata.is3D()) {
        if (faceId == 0 || faceId == 1) {
            const int i = faceId == 0 ? 0 : metadata.nx - 1;
            for (int j = 0; j < metadata.ny; ++j) {
                add(i, j, 0);
            }
        } else {
            const int j = faceId == 2 ? 0 : metadata.ny - 1;
            for (int i = 0; i < metadata.nx; ++i) {
                add(i, j, 0);
            }
        }
        return indices;
    }
    if (faceId == 0 || faceId == 1) {
        const int i = faceId == 0 ? 0 : metadata.nx - 1;
        for (int k = 0; k < metadata.nz; ++k) {
            for (int j = 0; j < metadata.ny; ++j) {
                add(i, j, k);
            }
        }
    } else if (faceId == 2 || faceId == 3) {
        const int j = faceId == 2 ? 0 : metadata.ny - 1;
        for (int k = 0; k < metadata.nz; ++k) {
            for (int i = 0; i < metadata.nx; ++i) {
                add(i, j, k);
            }
        }
    } else {
        const int k = faceId == 4 ? 0 : metadata.nz - 1;
        for (int j = 0; j < metadata.ny; ++j) {
            for (int i = 0; i < metadata.nx; ++i) {
                add(i, j, k);
            }
        }
    }
    return indices;
}

struct NekDisplayFace {
    std::size_t elementIndex = 0;
    int faceId = -1;
    std::vector<std::array<double, 3>> sortedCoordinates;
    // Coordinates alone cannot distinguish an ordinary shared face from an
    // overlapping or periodic copy.  This outward direction makes the
    // display-only exterior decision conservative.
    std::array<double, 3> outwardNormal{};
    bool hasOutwardNormal = false;
};

std::array<double, 3> SubtractCoordinates(
    const std::array<double, 3>& lhs,
    const std::array<double, 3>& rhs) {
    return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

std::array<double, 3> CrossCoordinates(const std::array<double, 3>& lhs,
                                        const std::array<double, 3>& rhs) {
    return {lhs[1] * rhs[2] - lhs[2] * rhs[1],
            lhs[2] * rhs[0] - lhs[0] * rhs[2],
            lhs[0] * rhs[1] - lhs[1] * rhs[0]};
}

double DotCoordinates(const std::array<double, 3>& lhs,
                      const std::array<double, 3>& rhs) {
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

double CoordinateNorm(const std::array<double, 3>& value) {
    return std::sqrt(DotCoordinates(value, value));
}

bool CoordinateLess(const std::array<double, 3>& lhs,
                    const std::array<double, 3>& rhs) {
    for (int component = 0; component < 3; ++component) {
        if (lhs[static_cast<std::size_t>(component)] !=
            rhs[static_cast<std::size_t>(component)]) {
            return lhs[static_cast<std::size_t>(component)] <
                   rhs[static_cast<std::size_t>(component)];
        }
    }
    return false;
}

std::string FaceSignature(const NekDisplayFace& face, double tolerance) {
    std::ostringstream stream;
    stream << face.sortedCoordinates.size() << ':';
    for (const std::array<double, 3>& point : face.sortedCoordinates) {
        for (const double coordinate : point) {
            if (!std::isfinite(coordinate)) {
                throw std::runtime_error(
                    "Nek5000 face matching requires finite coordinates.");
            }
            const long double scaled =
                static_cast<long double>(coordinate) /
                static_cast<long double>(tolerance);
            if (std::abs(scaled) >
                static_cast<long double>(std::numeric_limits<long long>::max() / 4)) {
                stream << std::setprecision(17) << coordinate << ',';
            } else {
                stream << std::llround(scaled) << ',';
            }
        }
    }
    return stream.str();
}

bool FacesCoincide(const NekDisplayFace& lhs,
                   const NekDisplayFace& rhs,
                   double tolerance) {
    if (lhs.sortedCoordinates.size() != rhs.sortedCoordinates.size()) {
        return false;
    }
    for (std::size_t point = 0; point < lhs.sortedCoordinates.size(); ++point) {
        for (int component = 0; component < 3; ++component) {
            if (std::abs(lhs.sortedCoordinates[point][static_cast<std::size_t>(component)] -
                         rhs.sortedCoordinates[point][static_cast<std::size_t>(component)]) >
                tolerance) {
                return false;
            }
        }
    }
    return true;
}

bool FacesHaveOpposedOutwardNormals(const NekDisplayFace& lhs,
                                    const NekDisplayFace& rhs) {
    // A true element interface has opposing outward normals.  If the normals
    // are unavailable or ambiguous, keep the faces visible rather than
    // silently erasing a possible periodic or overlapping boundary.
    return lhs.hasOutwardNormal && rhs.hasOutwardNormal &&
           DotCoordinates(lhs.outwardNormal, rhs.outwardNormal) < -0.95;
}

std::vector<NekDisplayFace> BuildNekDisplayFaces(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata) {
    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(canonicalDataSet);
    if (pointSet == nullptr || pointSet->GetPoints() == nullptr) {
        throw std::runtime_error(
            "Nek5000 display proxy requires explicit canonical point coordinates.");
    }
    const std::size_t pointsPerElement = metadata.pointsPerElement();
    std::vector<NekDisplayFace> faces;
    faces.reserve(metadata.elementCount *
                  static_cast<std::size_t>(NekFaceCount(metadata)));
    for (std::size_t element = 0; element < metadata.elementCount; ++element) {
        std::array<double, 3> elementCenter{};
        for (std::size_t localIndex = 0; localIndex < pointsPerElement;
             ++localIndex) {
            std::array<double, 3> point{};
            pointSet->GetPoint(static_cast<vtkIdType>(
                                   element * pointsPerElement + localIndex),
                               point.data());
            for (int axis = 0; axis < 3; ++axis) {
                elementCenter[static_cast<std::size_t>(axis)] +=
                    point[static_cast<std::size_t>(axis)];
            }
        }
        for (double& coordinate : elementCenter) {
            coordinate /= static_cast<double>(pointsPerElement);
        }
        for (int faceId = 0; faceId < NekFaceCount(metadata); ++faceId) {
            NekDisplayFace face;
            face.elementIndex = element;
            face.faceId = faceId;
            const std::vector<std::size_t> localIndices =
                FaceLocalIndices(metadata, faceId);
            std::vector<std::array<double, 3>> faceCoordinates;
            faceCoordinates.reserve(localIndices.size());
            for (const std::size_t localIndex : localIndices) {
                std::array<double, 3> point{};
                pointSet->GetPoint(static_cast<vtkIdType>(
                                       element * pointsPerElement + localIndex),
                                   point.data());
                faceCoordinates.push_back(point);
            }
            if (metadata.is3D()) {
                const std::array<int, 2> dimensions = FaceSampleDimensions(
                    metadata, faceId, metadata.nx, metadata.ny, metadata.nz);
                if (dimensions[0] >= 2 && dimensions[1] >= 2) {
                    const std::array<double, 3>& origin = faceCoordinates.front();
                    const std::array<double, 3>& firstAxis = faceCoordinates.at(
                        static_cast<std::size_t>(dimensions[0] - 1));
                    const std::array<double, 3>& secondAxis = faceCoordinates.at(
                        static_cast<std::size_t>((dimensions[1] - 1) * dimensions[0]));
                    std::array<double, 3> normal = CrossCoordinates(
                        SubtractCoordinates(firstAxis, origin),
                        SubtractCoordinates(secondAxis, origin));
                    std::array<double, 3> faceCenter{};
                    for (const std::array<double, 3>& point : faceCoordinates) {
                        for (int axis = 0; axis < 3; ++axis) {
                            faceCenter[static_cast<std::size_t>(axis)] +=
                                point[static_cast<std::size_t>(axis)];
                        }
                    }
                    for (double& coordinate : faceCenter) {
                        coordinate /= static_cast<double>(faceCoordinates.size());
                    }
                    const std::array<double, 3> towardFace =
                        SubtractCoordinates(faceCenter, elementCenter);
                    const double normalLength = CoordinateNorm(normal);
                    const double towardFaceLength = CoordinateNorm(towardFace);
                    const double orientation = DotCoordinates(normal, towardFace);
                    if (std::isfinite(normalLength) &&
                        std::isfinite(towardFaceLength) &&
                        std::isfinite(orientation) &&
                        normalLength > 0.0 && towardFaceLength > 0.0 &&
                        std::abs(orientation) >
                            1.0e-10 * normalLength * towardFaceLength) {
                        if (orientation < 0.0) {
                            for (double& component : normal) {
                                component = -component;
                            }
                        }
                        for (int axis = 0; axis < 3; ++axis) {
                            face.outwardNormal[static_cast<std::size_t>(axis)] =
                                normal[static_cast<std::size_t>(axis)] / normalLength;
                        }
                        face.hasOutwardNormal = true;
                    }
                }
            }
            face.sortedCoordinates = std::move(faceCoordinates);
            std::sort(face.sortedCoordinates.begin(), face.sortedCoordinates.end(),
                      CoordinateLess);
            faces.push_back(std::move(face));
        }
    }
    return faces;
}

std::vector<NekDisplayFace> ExteriorNekDisplayFaces(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata,
    const std::vector<NekDisplayFace>& faces) {
    double bounds[6] = {};
    canonicalDataSet->GetBounds(bounds);
    const double scale = std::max(
        {1.0, std::abs(bounds[1] - bounds[0]), std::abs(bounds[3] - bounds[2]),
         std::abs(bounds[5] - bounds[4])});
    // Matching is only an exterior-face visibility decision.  It never changes
    // canonical point IDs, and an uncertain match remains visible rather than
    // being silently collapsed.
    const double tolerance = 1.0e-11 * scale;
    std::unordered_map<std::string, std::vector<std::size_t>> buckets;
    buckets.reserve(faces.size());
    for (std::size_t index = 0; index < faces.size(); ++index) {
        buckets[FaceSignature(faces[index], tolerance)].push_back(index);
    }
    std::vector<bool> interior(faces.size(), false);
    for (const auto& [unused, candidates] : buckets) {
        (void)unused;
        // Non-manifold or duplicate-coordinate buckets are ambiguous.  Hide a
        // face only for exactly one mutually opposed pair; retain periodic or
        // overlapping candidates rather than treating them as interior by
        // coordinate coincidence alone.
        if (candidates.size() != 2) {
            continue;
        }
        const NekDisplayFace& first = faces[candidates[0]];
        const NekDisplayFace& second = faces[candidates[1]];
        if (first.elementIndex != second.elementIndex &&
            FacesCoincide(first, second, tolerance) &&
            FacesHaveOpposedOutwardNormals(first, second)) {
            interior[candidates[0]] = true;
            interior[candidates[1]] = true;
        }
    }
    std::vector<NekDisplayFace> exterior;
    exterior.reserve(faces.size());
    for (std::size_t index = 0; index < faces.size(); ++index) {
        if (!interior[index]) {
            exterior.push_back(faces[index]);
        }
    }
    return exterior;
}

void AttachNekSpectralSurfaceNormals(
    vtkPolyData* surface,
    const Nek5000SpectralMetadata& metadata,
    const Nek5000SpectralPointEvaluator& geometryEvaluator) {
    if (surface == nullptr || surface->GetPointData() == nullptr) {
        throw std::runtime_error(
            "Cannot prepare Nek5000 PBR normals without a surface point-data table.");
    }
    vtkDataArray* elementIndices = surface->GetPointData()->GetArray(
        kNek5000DisplayElementIndexArray);
    vtkDataArray* references = surface->GetPointData()->GetArray(
        kNek5000DisplayReferenceCoordinatesArray);
    vtkDataArray* faceIds = surface->GetPointData()->GetArray(
        kNek5000DisplayFaceIdArray);
    if (elementIndices == nullptr || references == nullptr ||
        references->GetNumberOfComponents() != 3 ||
        elementIndices->GetNumberOfTuples() != surface->GetNumberOfPoints() ||
        references->GetNumberOfTuples() != surface->GetNumberOfPoints()) {
        throw std::runtime_error(
            "Nek5000 PBR surface lost element/reference provenance.");
    }

    vtkNew<vtkDoubleArray> normals;
    normals->SetName(kNek5000SpectralSurfaceNormalArray);
    normals->SetNumberOfComponents(3);
    normals->SetNumberOfTuples(surface->GetNumberOfPoints());
    std::vector<std::array<double, 3>> elementCenters(metadata.elementCount);
    std::vector<bool> haveCenter(metadata.elementCount, false);

    for (vtkIdType pointId = 0; pointId < surface->GetNumberOfPoints(); ++pointId) {
        const int elementIndex = static_cast<int>(std::lround(
            elementIndices->GetComponent(pointId, 0)));
        if (elementIndex < 0 ||
            static_cast<std::size_t>(elementIndex) >= metadata.elementCount) {
            throw std::runtime_error(
                "Nek5000 PBR surface contains an invalid source element index.");
        }
        std::array<double, 3> reference{};
        references->GetTuple(pointId, reference.data());
        std::array<double, 3> physical{};
        std::array<double, 3> derivativeR{};
        std::array<double, 3> derivativeS{};
        std::array<double, 3> derivativeT{};
        if (!geometryEvaluator.EvaluateGeometryReference(
                static_cast<std::size_t>(elementIndex), reference.data(),
                &physical, &derivativeR, &derivativeS, &derivativeT)) {
            throw std::runtime_error(
                "Failed to evaluate a Nek5000 spectral PBR surface normal.");
        }

        std::array<double, 3> normal{};
        const int faceId = faceIds == nullptr ? -1 : static_cast<int>(std::lround(
            faceIds->GetComponent(pointId, 0)));
        if (!metadata.is3D()) {
            normal = CrossCoordinates(derivativeR, derivativeS);
        } else {
            switch (faceId) {
            case 0:
                normal = CrossCoordinates(derivativeS, derivativeT);
                for (double& component : normal) {
                    component = -component;
                }
                break;
            case 1:
                normal = CrossCoordinates(derivativeS, derivativeT);
                break;
            case 2:
                normal = CrossCoordinates(derivativeT, derivativeR);
                for (double& component : normal) {
                    component = -component;
                }
                break;
            case 3:
                normal = CrossCoordinates(derivativeT, derivativeR);
                break;
            case 4:
                normal = CrossCoordinates(derivativeR, derivativeS);
                for (double& component : normal) {
                    component = -component;
                }
                break;
            case 5:
                normal = CrossCoordinates(derivativeR, derivativeS);
                break;
            default:
                throw std::runtime_error(
                    "Nek5000 PBR surface contains no exterior-face provenance.");
            }
            const std::size_t element = static_cast<std::size_t>(elementIndex);
            if (!haveCenter[element]) {
                constexpr double centerReference[3] = {0.0, 0.0, 0.0};
                std::array<double, 3> ignoredR{};
                std::array<double, 3> ignoredS{};
                std::array<double, 3> ignoredT{};
                if (!geometryEvaluator.EvaluateGeometryReference(
                        element, centerReference, &elementCenters[element],
                        &ignoredR, &ignoredS, &ignoredT)) {
                    throw std::runtime_error(
                        "Failed to evaluate a Nek5000 spectral element centre.");
                }
                haveCenter[element] = true;
            }
            const std::array<double, 3> outward = SubtractCoordinates(
                physical, elementCenters[element]);
            if (DotCoordinates(normal, outward) < 0.0) {
                for (double& component : normal) {
                    component = -component;
                }
            }
        }
        const double magnitude = CoordinateNorm(normal);
        if (!std::isfinite(magnitude) || magnitude <= 1.0e-30) {
            throw std::runtime_error(
                "Nek5000 PBR surface has a degenerate spectral mapping normal.");
        }
        for (double& component : normal) {
            component /= magnitude;
        }
        normals->SetTuple(pointId, normal.data());
    }
    surface->GetPointData()->SetNormals(normals);
}

vtkSmartPointer<vtkDataSet> BuildNekDisplayVolume(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata,
    const Nek5000SpectralPointEvaluator& geometryEvaluator,
    int samplesX,
    int samplesY,
    int samplesZ) {
    vtkSmartPointer<vtkUnstructuredGrid> volume =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
    volume->SetPoints(NewDoublePoints());
    NekDisplayPointWriter points(volume, canonicalDataSet, metadata,
                                 geometryEvaluator, true);
    NekDisplayCellProvenance cells;
    const std::size_t pointsPerElement =
        static_cast<std::size_t>(samplesX) * static_cast<std::size_t>(samplesY) *
        static_cast<std::size_t>(samplesZ);
    for (std::size_t element = 0; element < metadata.elementCount; ++element) {
        const vtkIdType base = volume->GetNumberOfPoints();
        for (int k = 0; k < samplesZ; ++k) {
            for (int j = 0; j < samplesY; ++j) {
                for (int i = 0; i < samplesX; ++i) {
                    points.Add(element,
                               {UniformReferenceCoordinate(i, samplesX),
                                UniformReferenceCoordinate(j, samplesY),
                                metadata.is3D()
                                    ? UniformReferenceCoordinate(k, samplesZ)
                                    : 0.0});
                }
            }
        }
        (void)pointsPerElement;
        if (metadata.is3D()) {
            for (int k = 0; k < samplesZ - 1; ++k) {
                for (int j = 0; j < samplesY - 1; ++j) {
                    for (int i = 0; i < samplesX - 1; ++i) {
                        vtkIdType ids[8] = {
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j, k)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j, k)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j + 1, k)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j + 1, k)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j, k + 1)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j, k + 1)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j + 1, k + 1)),
                            base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j + 1, k + 1))};
                        volume->InsertNextCell(VTK_HEXAHEDRON, 8, ids);
                        cells.Add(metadata, element, -1);
                    }
                }
            }
        } else {
            for (int j = 0; j < samplesY - 1; ++j) {
                for (int i = 0; i < samplesX - 1; ++i) {
                    vtkIdType ids[4] = {
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j + 1, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j + 1, 0))};
                    volume->InsertNextCell(VTK_QUAD, 4, ids);
                    cells.Add(metadata, element, -1);
                }
            }
        }
    }
    cells.Attach(volume);
    MarkNek5000DisplayProxy(volume, kSemAwareDisplayRepresentation, nullptr,
                            0.0, &metadata);
    return volume;
}

vtkSmartPointer<vtkPolyData> BuildNekDisplaySurface(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata,
    const Nek5000SpectralPointEvaluator& geometryEvaluator,
    const std::vector<NekDisplayFace>& exteriorFaces,
    int samplesX,
    int samplesY,
    int samplesZ) {
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->SetPoints(NewDoublePoints());
    NekDisplayPointWriter points(surface, canonicalDataSet, metadata,
                                 geometryEvaluator, true);
    NekDisplayCellProvenance cells;

    if (!metadata.is3D()) {
        // A two-dimensional SEM domain is itself the visible surface. Its
        // element boundaries have no area, so retain sampled element quads.
        vtkNew<vtkCellArray> polygons;
        for (std::size_t element = 0; element < metadata.elementCount; ++element) {
            const vtkIdType base = surface->GetNumberOfPoints();
            for (int j = 0; j < samplesY; ++j) {
                for (int i = 0; i < samplesX; ++i) {
                    points.Add(element,
                               {UniformReferenceCoordinate(i, samplesX),
                                UniformReferenceCoordinate(j, samplesY), 0.0});
                }
            }
            for (int j = 0; j < samplesY - 1; ++j) {
                for (int i = 0; i < samplesX - 1; ++i) {
                    vtkIdType ids[4] = {
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i + 1, j + 1, 0)),
                        base + static_cast<vtkIdType>(SampleLocalIndex(samplesX, samplesY, i, j + 1, 0))};
                    polygons->InsertNextCell(4, ids);
                    cells.Add(metadata, element, -1);
                }
            }
        }
        surface->SetPolys(polygons);
    } else {
        vtkNew<vtkCellArray> polygons;
        for (const NekDisplayFace& face : exteriorFaces) {
            const std::array<int, 2> dimensions = FaceSampleDimensions(
                metadata, face.faceId, samplesX, samplesY, samplesZ);
            const vtkIdType base = surface->GetNumberOfPoints();
            for (int second = 0; second < dimensions[1]; ++second) {
                for (int first = 0; first < dimensions[0]; ++first) {
                    points.Add(face.elementIndex,
                               FaceReference(
                                   metadata, face.faceId,
                                   UniformReferenceCoordinate(first, dimensions[0]),
                                   UniformReferenceCoordinate(second, dimensions[1])),
                               -1,
                               face.faceId);
                }
            }
            for (int second = 0; second < dimensions[1] - 1; ++second) {
                for (int first = 0; first < dimensions[0] - 1; ++first) {
                    const vtkIdType lower =
                        base + static_cast<vtkIdType>(second * dimensions[0] + first);
                    vtkIdType ids[4] = {lower, lower + 1,
                                         lower + dimensions[0] + 1,
                                         lower + dimensions[0]};
                    polygons->InsertNextCell(4, ids);
                    cells.Add(metadata, face.elementIndex, face.faceId);
                }
            }
        }
        surface->SetPolys(polygons);
    }
    cells.Attach(surface);
    // Exterior PBR normals use the exact spectral mapping derivatives at the
    // evaluator-sampled surface points, not normals inferred from GLL
    // subcells or from the raw nodal cloud.
    AttachNekSpectralSurfaceNormals(surface, metadata, geometryEvaluator);
    MarkNek5000DisplayProxy(surface, kSemAwareDisplayRepresentation, nullptr,
                            0.0, &metadata);
    return surface;
}

vtkSmartPointer<vtkPolyData> BuildNekElementGrid(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata,
    const Nek5000SpectralPointEvaluator& geometryEvaluator,
    const std::vector<NekDisplayFace>& allFaces,
    int samplesX,
    int samplesY,
    int samplesZ) {
    vtkSmartPointer<vtkPolyData> grid = vtkSmartPointer<vtkPolyData>::New();
    grid->SetPoints(NewDoublePoints());
    NekDisplayPointWriter points(grid, canonicalDataSet, metadata,
                                 geometryEvaluator, false);
    vtkNew<vtkCellArray> lines;
    NekDisplayCellProvenance cells;
    const std::vector<double> nodesX = GllNodes(metadata.nx);
    const std::vector<double> nodesY = GllNodes(metadata.ny);
    const std::vector<double> nodesZ = metadata.is3D() ? GllNodes(metadata.nz)
                                                        : std::vector<double>{0.0};

    const auto addLine = [&](std::size_t elementIndex,
                             int faceId,
                             const std::vector<std::array<double, 3>>& references) {
        std::vector<vtkIdType> ids;
        ids.reserve(references.size());
        for (const std::array<double, 3>& reference : references) {
            ids.push_back(points.Add(elementIndex, reference, -1, faceId));
        }
        lines->InsertNextCell(static_cast<vtkIdType>(ids.size()), ids.data());
        cells.Add(metadata, elementIndex, faceId);
    };

    for (const NekDisplayFace& face : allFaces) {
        const std::array<int, 2> dimensions = FaceSampleDimensions(
            metadata, face.faceId, samplesX, samplesY, samplesZ);
        const std::vector<double>* firstNodes = nullptr;
        const std::vector<double>* secondNodes = nullptr;
        if (!metadata.is3D()) {
            firstNodes = (face.faceId == 0 || face.faceId == 1) ? &nodesY : &nodesX;
        } else if (face.faceId == 0 || face.faceId == 1) {
            firstNodes = &nodesY;
            secondNodes = &nodesZ;
        } else if (face.faceId == 2 || face.faceId == 3) {
            firstNodes = &nodesX;
            secondNodes = &nodesZ;
        } else {
            firstNodes = &nodesX;
            secondNodes = &nodesY;
        }
        if (!metadata.is3D()) {
            std::vector<std::array<double, 3>> reference;
            reference.reserve(static_cast<std::size_t>(dimensions[0]));
            for (int first = 0; first < dimensions[0]; ++first) {
                reference.push_back(FaceReference(
                    metadata, face.faceId,
                    UniformReferenceCoordinate(first, dimensions[0])));
            }
            addLine(face.elementIndex, face.faceId, reference);
            continue;
        }
        for (const double fixedFirst : *firstNodes) {
            std::vector<std::array<double, 3>> reference;
            reference.reserve(static_cast<std::size_t>(dimensions[1]));
            for (int second = 0; second < dimensions[1]; ++second) {
                reference.push_back(FaceReference(
                    metadata, face.faceId, fixedFirst,
                    UniformReferenceCoordinate(second, dimensions[1])));
            }
            addLine(face.elementIndex, face.faceId, reference);
        }
        for (const double fixedSecond : *secondNodes) {
            std::vector<std::array<double, 3>> reference;
            reference.reserve(static_cast<std::size_t>(dimensions[0]));
            for (int first = 0; first < dimensions[0]; ++first) {
                reference.push_back(FaceReference(
                    metadata, face.faceId,
                    UniformReferenceCoordinate(first, dimensions[0]), fixedSecond));
            }
            addLine(face.elementIndex, face.faceId, reference);
        }
    }
    grid->SetLines(lines);
    cells.Attach(grid);
    MarkNek5000DisplayProxy(grid, kSemAwareDisplayRepresentation, nullptr,
                            0.0, &metadata);
    return grid;
}

vtkSmartPointer<vtkPolyData> BuildNekGllNodes(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralMetadata& metadata) {
    vtkPointSet* canonicalPoints = vtkPointSet::SafeDownCast(canonicalDataSet);
    if (canonicalPoints == nullptr || canonicalPoints->GetPoints() == nullptr) {
        throw std::runtime_error("Nek5000 GLL overlay requires explicit points.");
    }
    vtkSmartPointer<vtkPolyData> nodes = vtkSmartPointer<vtkPolyData>::New();
    nodes->SetPoints(NewDoublePoints());
    vtkNew<vtkCellArray> vertices;
    vtkNew<vtkIntArray> elementIds;
    vtkNew<vtkIntArray> elementIndices;
    vtkNew<vtkIntArray> localIndices;
    vtkNew<vtkDoubleArray> referenceCoordinates;
    elementIds->SetName(kNek5000DisplayElementIdArray);
    elementIndices->SetName(kNek5000DisplayElementIndexArray);
    localIndices->SetName(kNek5000DisplayLocalGllIndexArray);
    referenceCoordinates->SetName(kNek5000DisplayReferenceCoordinatesArray);
    referenceCoordinates->SetNumberOfComponents(3);
    const std::vector<double> nodesX = GllNodes(metadata.nx);
    const std::vector<double> nodesY = GllNodes(metadata.ny);
    const std::vector<double> nodesZ = metadata.is3D() ? GllNodes(metadata.nz)
                                                        : std::vector<double>{0.0};
    const std::size_t pointsPerElement = metadata.pointsPerElement();
    for (std::size_t element = 0; element < metadata.elementCount; ++element) {
        for (int k = 0; k < metadata.nz; ++k) {
            for (int j = 0; j < metadata.ny; ++j) {
                for (int i = 0; i < metadata.nx; ++i) {
                    const std::size_t local = NekLocalIndex(metadata.nx, metadata.ny,
                                                            i, j, k);
                    std::array<double, 3> coordinate{};
                    canonicalPoints->GetPoint(static_cast<vtkIdType>(
                                                  element * pointsPerElement + local),
                                              coordinate.data());
                    const vtkIdType id = nodes->GetPoints()->InsertNextPoint(
                        coordinate.data());
                    vertices->InsertNextCell(1, &id);
                    elementIds->InsertNextValue(metadata.elementIds[element]);
                    elementIndices->InsertNextValue(static_cast<int>(element));
                    localIndices->InsertNextValue(static_cast<int>(local));
                    const double reference[3] = {
                        nodesX[static_cast<std::size_t>(i)],
                        nodesY[static_cast<std::size_t>(j)],
                        nodesZ[static_cast<std::size_t>(k)]};
                    referenceCoordinates->InsertNextTuple(reference);
                }
            }
        }
    }
    nodes->SetVerts(vertices);
    nodes->GetPointData()->AddArray(elementIds);
    nodes->GetPointData()->AddArray(elementIndices);
    nodes->GetPointData()->AddArray(localIndices);
    nodes->GetPointData()->AddArray(referenceCoordinates);
    MarkNek5000DisplayProxy(nodes, kSemAwareDisplayRepresentation, nullptr,
                            0.0, &metadata);
    return nodes;
}

}  // namespace

bool IsNek5000SpectralDisplayProxy(vtkDataSet* dataSet) {
    if (dataSet == nullptr || dataSet->GetFieldData() == nullptr) {
        return false;
    }
    vtkStringArray* marker = vtkStringArray::SafeDownCast(
        dataSet->GetFieldData()->GetAbstractArray(kNek5000DisplayProxyField));
    return marker != nullptr && marker->GetNumberOfTuples() == 1 &&
           marker->GetValue(0) == "1";
}

Nek5000SpectralDisplayProxy CreateNek5000SpectralDisplayProxy(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralDisplayOptions& options) {
    if (canonicalDataSet == nullptr ||
        !HasNek5000SpectralMetadata(canonicalDataSet) ||
        IsNek5000SpectralDisplayProxy(canonicalDataSet)) {
        throw std::runtime_error(
            "Nek5000 display proxies require canonical element-local GLL data.");
    }
    const Nek5000SpectralMetadata metadata =
        GetNek5000SpectralMetadata(canonicalDataSet);
    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(canonicalDataSet);
    vtkDataArray* coordinates =
        pointSet == nullptr || pointSet->GetPoints() == nullptr
            ? nullptr
            : pointSet->GetPoints()->GetData();
    if (coordinates == nullptr) {
        throw std::runtime_error(
            "Nek5000 display proxies require a canonical coordinate array.");
    }
    const int samplesX = SampleCountForGllAxis(metadata.nx, options);
    const int samplesY = SampleCountForGllAxis(metadata.ny, options);
    const int samplesZ = metadata.is3D()
                             ? SampleCountForGllAxis(metadata.nz, options)
                             : 1;
    Nek5000SpectralPointEvaluator geometryEvaluator(canonicalDataSet, coordinates);
    Nek5000SpectralDisplayProxy proxy;
    if (options.buildVolumeProxy) {
        proxy.volume = BuildNekDisplayVolume(canonicalDataSet, metadata,
                                             geometryEvaluator, samplesX, samplesY,
                                             samplesZ);
    }
    std::vector<NekDisplayFace> allFaces;
    if (options.buildExteriorSurface || options.buildElementGrid) {
        allFaces = BuildNekDisplayFaces(canonicalDataSet, metadata);
    }
    if (options.buildExteriorSurface) {
        proxy.exteriorSurface = BuildNekDisplaySurface(
            canonicalDataSet, metadata, geometryEvaluator,
            metadata.is3D()
                ? ExteriorNekDisplayFaces(canonicalDataSet, metadata, allFaces)
                : allFaces,
            samplesX, samplesY, samplesZ);
    }
    if (options.buildElementGrid) {
        proxy.elementGrid = BuildNekElementGrid(canonicalDataSet, metadata,
                                                geometryEvaluator, allFaces,
                                                samplesX, samplesY, samplesZ);
    }
    if (options.buildGllNodes) {
        proxy.gllNodes = BuildNekGllNodes(canonicalDataSet, metadata);
    }
    return proxy;
}

namespace {

const char* NekSpectralOperationName(Nek5000SpectralOperationKind kind) {
    switch (kind) {
    case Nek5000SpectralOperationKind::Crop:
        return "crop";
    case Nek5000SpectralOperationKind::Clip:
        return "clip";
    case Nek5000SpectralOperationKind::Slice:
        return "slice";
    case Nek5000SpectralOperationKind::Contour:
        return "contour";
    }
    throw std::runtime_error("Unknown Native Nek5000 spectral operation.");
}

vtkSmartPointer<vtkDataSet> DeepCopyNekOperationOutput(vtkDataSet* source) {
    if (source == nullptr) {
        return nullptr;
    }
    vtkDataSet* rawCopy = vtkDataSet::SafeDownCast(source->NewInstance());
    if (rawCopy == nullptr) {
        return nullptr;
    }
    vtkSmartPointer<vtkDataSet> copy;
    copy.TakeReference(rawCopy);
    copy->DeepCopy(source);
    return copy;
}

vtkSmartPointer<vtkPolyData> DeepCopyNekOperationSurface(vtkPolyData* source) {
    if (source == nullptr) {
        return nullptr;
    }
    vtkSmartPointer<vtkPolyData> copy = vtkSmartPointer<vtkPolyData>::New();
    copy->DeepCopy(source);
    return copy;
}

vtkSmartPointer<vtkDataSet> RemoveUnusedNekOperationPoints(
    vtkDataSet* source) {
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(source);
    if (grid == nullptr) {
        return DeepCopyNekOperationOutput(source);
    }
    // vtkClipDataSet may leave rejected input points in its output point
    // table.  Removing only unused points preserves the non-merging contract
    // while ensuring bounds and downstream display filters describe the
    // actual evaluator-sampled operation geometry.
    vtkNew<vtkRemoveUnusedPoints> removeUnused;
    removeUnused->SetInputData(grid);
    removeUnused->Update();
    return DeepCopyNekOperationOutput(removeUnused->GetOutput());
}

std::array<double, 6> NormalizedNekCropBounds(
    const std::array<double, 6>& input) {
    std::array<double, 6> bounds = input;
    for (int axis = 0; axis < 3; ++axis) {
        const std::size_t lower = static_cast<std::size_t>(axis * 2);
        const std::size_t upper = lower + 1;
        if (!std::isfinite(bounds[lower]) || !std::isfinite(bounds[upper])) {
            throw std::runtime_error(
                "Native Nek5000 crop bounds must be finite.");
        }
        if (bounds[upper] < bounds[lower]) {
            std::swap(bounds[lower], bounds[upper]);
        }
        if (!(bounds[upper] > bounds[lower])) {
            const double scale = std::max(
                {1.0, std::abs(bounds[lower]), std::abs(bounds[upper])});
            const double padding = scale * 1.0e-9;
            bounds[lower] -= padding;
            bounds[upper] += padding;
        }
    }
    return bounds;
}

void ValidateNekPlane(const std::array<double, 3>& origin,
                      const std::array<double, 3>& normal) {
    double lengthSquared = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        if (!std::isfinite(origin[static_cast<std::size_t>(axis)]) ||
            !std::isfinite(normal[static_cast<std::size_t>(axis)])) {
            throw std::runtime_error(
                "Native Nek5000 operation plane must be finite.");
        }
        lengthSquared += normal[static_cast<std::size_t>(axis)] *
                         normal[static_cast<std::size_t>(axis)];
    }
    if (!std::isfinite(lengthSquared) ||
        lengthSquared <= std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error(
            "Native Nek5000 operation plane normal must be non-zero.");
    }
}

void RequireNekOperationPointProvenance(vtkDataSet* dataSet) {
    if (dataSet == nullptr || dataSet->GetPointData() == nullptr) {
        throw std::runtime_error(
            "Native Nek5000 operation produced no point-data provenance.");
    }
    vtkDataArray* elementIndices = dataSet->GetPointData()->GetArray(
        kNek5000DisplayElementIndexArray);
    vtkDataArray* references = dataSet->GetPointData()->GetArray(
        kNek5000DisplayReferenceCoordinatesArray);
    if (elementIndices == nullptr || references == nullptr ||
        elementIndices->GetNumberOfTuples() != dataSet->GetNumberOfPoints() ||
        references->GetNumberOfComponents() != 3 ||
        references->GetNumberOfTuples() != dataSet->GetNumberOfPoints()) {
        throw std::runtime_error(
            "Native Nek5000 operation lost unmerged element/reference provenance.");
    }
}

vtkDataArray* PrepareNekContourScalar(
    vtkDataSet* sampledData,
    const Nek5000SpectralOperationRequest& request) {
    if (sampledData == nullptr || sampledData->GetPointData() == nullptr ||
        request.contourArrayName.empty()) {
        throw std::runtime_error(
            "Native Nek5000 contour requires a canonical point-data field.");
    }
    vtkPointData* pointData = sampledData->GetPointData();
    vtkDataArray* input = pointData->GetArray(request.contourArrayName.c_str());
    if (input == nullptr || input->GetNumberOfComponents() < 1 ||
        input->GetNumberOfTuples() != sampledData->GetNumberOfPoints()) {
        throw std::runtime_error(
            "Native Nek5000 contour field is not available on evaluator samples: " +
            request.contourArrayName);
    }
    if (request.contourComponent == 0 && input->GetNumberOfComponents() == 1) {
        return input;
    }
    if (request.contourComponent < -1 ||
        request.contourComponent >= input->GetNumberOfComponents()) {
        throw std::runtime_error(
            "Native Nek5000 contour component is outside the source field range.");
    }

    constexpr const char* kContourScalarName =
        "SCP_Nek5000SpectralOperationContourScalar";
    pointData->RemoveArray(kContourScalarName);
    vtkNew<vtkDoubleArray> scalar;
    scalar->SetName(kContourScalarName);
    scalar->SetNumberOfComponents(1);
    scalar->SetNumberOfTuples(input->GetNumberOfTuples());
    for (vtkIdType tuple = 0; tuple < input->GetNumberOfTuples(); ++tuple) {
        double value = 0.0;
        if (request.contourComponent >= 0) {
            value = input->GetComponent(tuple, request.contourComponent);
        } else {
            for (int component = 0;
                 component < input->GetNumberOfComponents(); ++component) {
                const double componentValue = input->GetComponent(tuple, component);
                value += componentValue * componentValue;
            }
            value = std::sqrt(value);
        }
        if (!std::isfinite(value)) {
            throw std::runtime_error(
                "Native Nek5000 contour field contains a non-finite evaluator sample.");
        }
        scalar->SetValue(tuple, value);
    }
    pointData->AddArray(scalar);
    return pointData->GetArray(kContourScalarName);
}

vtkSmartPointer<vtkDataSet> ApplyNekSpectralOperation(
    vtkDataSet* sampledInput,
    const Nek5000SpectralOperationRequest& request) {
    if (sampledInput == nullptr || sampledInput->GetNumberOfPoints() <= 0) {
        throw std::runtime_error(
            "Native Nek5000 operation has no evaluator-sampled input.");
    }
    vtkSmartPointer<vtkDataSet> output;
    if (request.kind == Nek5000SpectralOperationKind::Crop) {
        const std::array<double, 6> bounds =
            NormalizedNekCropBounds(request.cropBounds);
        vtkNew<vtkBox> box;
        box->SetBounds(bounds.data());
        vtkNew<vtkClipDataSet> crop;
        vtkNew<vtkNonMergingPointLocator> locator;
        crop->SetLocator(locator);
        crop->SetInputData(sampledInput);
        crop->SetClipFunction(box);
        crop->SetValue(0.0);
        crop->InsideOutOn();
        crop->Update();
        output = DeepCopyNekOperationOutput(
            vtkDataSet::SafeDownCast(crop->GetOutputDataObject(0)));
    } else if (request.kind == Nek5000SpectralOperationKind::Clip) {
        ValidateNekPlane(request.planeOrigin, request.planeNormal);
        vtkNew<vtkPlane> plane;
        plane->SetOrigin(request.planeOrigin.data());
        plane->SetNormal(request.planeNormal.data());
        vtkNew<vtkClipDataSet> clip;
        vtkNew<vtkNonMergingPointLocator> locator;
        clip->SetLocator(locator);
        clip->SetInputData(sampledInput);
        clip->SetClipFunction(plane);
        clip->SetValue(0.0);
        clip->InsideOutOn();
        clip->Update();
        output = DeepCopyNekOperationOutput(
            vtkDataSet::SafeDownCast(clip->GetOutputDataObject(0)));
    } else if (request.kind == Nek5000SpectralOperationKind::Slice) {
        ValidateNekPlane(request.planeOrigin, request.planeNormal);
        vtkNew<vtkPlane> plane;
        plane->SetOrigin(request.planeOrigin.data());
        plane->SetNormal(request.planeNormal.data());
        vtkNew<vtkCutter> cutter;
        vtkNew<vtkNonMergingPointLocator> locator;
        cutter->SetLocator(locator);
        cutter->SetInputData(sampledInput);
        cutter->SetCutFunction(plane);
        cutter->Update();
        output = DeepCopyNekOperationOutput(
            vtkDataSet::SafeDownCast(cutter->GetOutputDataObject(0)));
    } else if (request.kind == Nek5000SpectralOperationKind::Contour) {
        if (request.contourValues.empty()) {
            throw std::runtime_error(
                "Native Nek5000 contour requires at least one finite isovalue.");
        }
        vtkDataArray* scalar = PrepareNekContourScalar(sampledInput, request);
        if (scalar == nullptr || scalar->GetName() == nullptr) {
            throw std::runtime_error(
                "Native Nek5000 contour could not prepare an evaluator-sampled scalar.");
        }
        vtkNew<vtkContourFilter> contour;
        vtkNew<vtkNonMergingPointLocator> locator;
        contour->SetLocator(locator);
        contour->SetInputData(sampledInput);
        contour->SetInputArrayToProcess(0, 0, 0,
                                        vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                        scalar->GetName());
        contour->SetNumberOfContours(
            static_cast<int>(request.contourValues.size()));
        for (std::size_t index = 0; index < request.contourValues.size(); ++index) {
            const double value = request.contourValues[index];
            if (!std::isfinite(value)) {
                throw std::runtime_error(
                    "Native Nek5000 contour isovalue must be finite.");
            }
            contour->SetValue(static_cast<int>(index), value);
        }
        contour->Update();
        output = DeepCopyNekOperationOutput(
            vtkDataSet::SafeDownCast(contour->GetOutputDataObject(0)));
    }
    if (output == nullptr || output->GetNumberOfPoints() <= 0) {
        throw std::runtime_error(
            std::string("Native Nek5000 ") + NekSpectralOperationName(request.kind) +
            " produced no visible evaluator-sampled data.");
    }
    output = RemoveUnusedNekOperationPoints(output);
    if (output == nullptr || output->GetNumberOfPoints() <= 0) {
        throw std::runtime_error(
            std::string("Native Nek5000 ") + NekSpectralOperationName(request.kind) +
            " lost all evaluator-sampled points while removing unused VTK points.");
    }
    RequireNekOperationPointProvenance(output);
    return output;
}

vtkSmartPointer<vtkPolyData> PrepareNekSampledPbrSurface(
    vtkDataSet* sampledData) {
    vtkSmartPointer<vtkPolyData> surface;
    if (vtkPolyData* direct = vtkPolyData::SafeDownCast(sampledData)) {
        surface = DeepCopyNekOperationSurface(direct);
    } else {
        vtkNew<vtkDataSetSurfaceFilter> extract;
        extract->SetInputData(sampledData);
        extract->Update();
        surface = DeepCopyNekOperationSurface(extract->GetOutput());
    }
    if (surface == nullptr || surface->GetNumberOfPoints() <= 0) {
        throw std::runtime_error(
            "Native Nek5000 PBR preparation produced no sampled surface.");
    }
    vtkNew<vtkPolyDataNormals> normals;
    normals->SetInputData(surface);
    normals->ComputePointNormalsOn();
    normals->ComputeCellNormalsOff();
    normals->SplittingOff();
    normals->ConsistencyOff();
    normals->AutoOrientNormalsOff();
    normals->FlipNormalsOff();
    normals->NonManifoldTraversalOff();
    normals->Update();
    vtkSmartPointer<vtkPolyData> pbrSurface =
        DeepCopyNekOperationSurface(normals->GetOutput());
    if (pbrSurface == nullptr || pbrSurface->GetPointData() == nullptr ||
        pbrSurface->GetPointData()->GetNormals() == nullptr) {
        throw std::runtime_error(
            "Native Nek5000 PBR preparation did not produce sampled normals.");
    }
    RequireNekOperationPointProvenance(pbrSurface);
    return pbrSurface;
}

}  // namespace

Nek5000SpectralOperationResult CreateNek5000SpectralOperation(
    vtkDataSet* canonicalDataSet,
    const Nek5000SpectralOperationRequest& request,
    const Nek5000SpectralOperationOptions& options) {
    return CreateNek5000SpectralOperationPipeline(canonicalDataSet,
                                                    {request}, options);
}

Nek5000SpectralOperationResult CreateNek5000SpectralOperationPipeline(
    vtkDataSet* canonicalDataSet,
    const std::vector<Nek5000SpectralOperationRequest>& requests,
    const Nek5000SpectralOperationOptions& options) {
    if (canonicalDataSet == nullptr ||
        !HasNek5000SpectralMetadata(canonicalDataSet) ||
        IsNek5000SpectralDisplayProxy(canonicalDataSet)) {
        throw std::runtime_error(
            "Native Nek5000 operations require canonical element-local GLL data.");
    }
    if (options.samplesPerGllInterval < 1 || options.maxSamplesPerAxis < 2) {
        throw std::runtime_error(
            "Native Nek5000 operation sampling options are invalid.");
    }
    const Nek5000SpectralMetadata metadata =
        GetNek5000SpectralMetadata(canonicalDataSet);
    Nek5000SpectralDisplayOptions sampling;
    sampling.samplesPerGllInterval = options.samplesPerGllInterval;
    sampling.maxSamplesPerAxis = options.maxSamplesPerAxis;
    sampling.buildVolumeProxy = true;
    sampling.buildExteriorSurface = options.preparePbrSurface && requests.empty();
    sampling.buildElementGrid = false;
    sampling.buildGllNodes = false;
    const int samplesX = SampleCountForGllAxis(metadata.nx, sampling);
    const int samplesY = SampleCountForGllAxis(metadata.ny, sampling);
    const int samplesZ = metadata.is3D()
        ? SampleCountForGllAxis(metadata.nz, sampling)
        : 1;
    const double maximumReferenceEdgeLength = std::max(
        2.0 / static_cast<double>(std::max(1, samplesX - 1)),
        std::max(2.0 / static_cast<double>(std::max(1, samplesY - 1)),
                 metadata.is3D()
                     ? 2.0 / static_cast<double>(std::max(1, samplesZ - 1))
                     : 0.0));
    Nek5000SpectralDisplayProxy proxy =
        CreateNek5000SpectralDisplayProxy(canonicalDataSet, sampling);
    vtkSmartPointer<vtkDataSet> current = proxy.volume;
    if (current == nullptr) {
        throw std::runtime_error(
            "Native Nek5000 operation sampling produced no volume proxy.");
    }

    std::ostringstream operationLabel;
    if (requests.empty()) {
        operationLabel << "surface";
    }
    for (std::size_t index = 0; index < requests.size(); ++index) {
        if (index != 0) {
            operationLabel << '>';
        }
        operationLabel << NekSpectralOperationName(requests[index].kind);
        current = ApplyNekSpectralOperation(current, requests[index]);
    }
    MarkNek5000DisplayProxy(current, kSemAwareOperationDisplayRepresentation,
                            operationLabel.str().c_str(),
                            maximumReferenceEdgeLength, &metadata);

    Nek5000SpectralOperationResult result;
    result.data = current;
    result.maximumReferenceEdgeLength = maximumReferenceEdgeLength;
    if (options.preparePbrSurface) {
        if (requests.empty()) {
            result.pbrSurface = proxy.exteriorSurface;
            if (result.pbrSurface == nullptr) {
                throw std::runtime_error(
                    "Native Nek5000 PBR preparation produced no exterior surface.");
            }
            MarkNek5000DisplayProxy(result.pbrSurface,
                                    kSemAwareOperationDisplayRepresentation,
                                    "pbr_surface",
                                    maximumReferenceEdgeLength, &metadata);
        } else {
            result.pbrSurface = PrepareNekSampledPbrSurface(current);
            MarkNek5000DisplayProxy(result.pbrSurface,
                                    kSemAwareOperationDisplayRepresentation,
                                    (operationLabel.str() + ">pbr_surface").c_str(),
                                    maximumReferenceEdgeLength, &metadata);
        }
    }
    return result;
}

bool ProbeNek5000SpectralField(
    vtkDataSet* canonicalDataSet,
    vtkDataArray* pointArray,
    const double physicalPoint[3],
    Nek5000SpectralProbeResult* result,
    const Nek5000PointEvaluationOptions& options) {
    if (result == nullptr || canonicalDataSet == nullptr || pointArray == nullptr ||
        !HasNek5000SpectralMetadata(canonicalDataSet) ||
        IsNek5000SpectralDisplayProxy(canonicalDataSet) || physicalPoint == nullptr) {
        return false;
    }
    Nek5000SpectralPointEvaluator evaluator(canonicalDataSet, pointArray, options);
    std::vector<double> values(static_cast<std::size_t>(
        evaluator.numberOfComponents()));
    int elementIndex = -1;
    std::array<double, 3> reference{};
    if (!evaluator.EvaluateWithReference(physicalPoint, values.data(), &elementIndex,
                                         &reference)) {
        return false;
    }
    const Nek5000SpectralMetadata metadata =
        GetNek5000SpectralMetadata(canonicalDataSet);
    if (elementIndex < 0 ||
        static_cast<std::size_t>(elementIndex) >= metadata.elementIds.size()) {
        return false;
    }
    result->elementIndex = elementIndex;
    result->elementId = metadata.elementIds[static_cast<std::size_t>(elementIndex)];
    result->referenceCoordinates = reference;
    result->values = std::move(values);
    return true;
}

bool HasSameNek5000SpectralGeometry(vtkDataSet* lhs,
                                    vtkDataSet* rhs,
                                    double tolerance) {
    if (!std::isfinite(tolerance) || tolerance < 0.0) {
        throw std::runtime_error(
            "Nek5000 spectral geometry tolerance must be finite and non-negative.");
    }
    if (lhs == nullptr || rhs == nullptr ||
        !HasNek5000SpectralMetadata(lhs) ||
        !HasNek5000SpectralMetadata(rhs)) {
        return false;
    }
    const Nek5000SpectralMetadata lhsMetadata =
        GetNek5000SpectralMetadata(lhs);
    const Nek5000SpectralMetadata rhsMetadata =
        GetNek5000SpectralMetadata(rhs);
    if (lhsMetadata.nx != rhsMetadata.nx ||
        lhsMetadata.ny != rhsMetadata.ny ||
        lhsMetadata.nz != rhsMetadata.nz ||
        lhsMetadata.elementCount != rhsMetadata.elementCount ||
        lhsMetadata.elementIds != rhsMetadata.elementIds ||
        lhs->GetNumberOfPoints() != rhs->GetNumberOfPoints()) {
        return false;
    }

    vtkPointSet* lhsPointSet = vtkPointSet::SafeDownCast(lhs);
    vtkPointSet* rhsPointSet = vtkPointSet::SafeDownCast(rhs);
    vtkDoubleArray* lhsCoordinates =
        lhsPointSet == nullptr || lhsPointSet->GetPoints() == nullptr
            ? nullptr
            : vtkDoubleArray::SafeDownCast(lhsPointSet->GetPoints()->GetData());
    vtkDoubleArray* rhsCoordinates =
        rhsPointSet == nullptr || rhsPointSet->GetPoints() == nullptr
            ? nullptr
            : vtkDoubleArray::SafeDownCast(rhsPointSet->GetPoints()->GetData());
    if (lhsCoordinates != nullptr && rhsCoordinates != nullptr &&
        lhsCoordinates->GetNumberOfComponents() == 3 &&
        rhsCoordinates->GetNumberOfComponents() == 3) {
        const std::size_t valueCount =
            3 * static_cast<std::size_t>(lhs->GetNumberOfPoints());
        const double* lhsValues = lhsCoordinates->GetPointer(0);
        const double* rhsValues = rhsCoordinates->GetPointer(0);
        for (std::size_t value = 0; value < valueCount; ++value) {
            if (!std::isfinite(lhsValues[value]) ||
                !std::isfinite(rhsValues[value]) ||
                (tolerance == 0.0
                     ? lhsValues[value] != rhsValues[value]
                     : std::abs(lhsValues[value] - rhsValues[value]) >
                           tolerance)) {
                return false;
            }
        }
        return true;
    }

    double lhsPoint[3] = {};
    double rhsPoint[3] = {};
    for (vtkIdType pointId = 0; pointId < lhs->GetNumberOfPoints(); ++pointId) {
        lhs->GetPoint(pointId, lhsPoint);
        rhs->GetPoint(pointId, rhsPoint);
        for (int axis = 0; axis < 3; ++axis) {
            if (!std::isfinite(lhsPoint[axis]) ||
                !std::isfinite(rhsPoint[axis])) {
                return false;
            }
            if (tolerance == 0.0) {
                if (lhsPoint[axis] != rhsPoint[axis]) {
                    return false;
                }
            } else if (std::abs(lhsPoint[axis] - rhsPoint[axis]) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

vtkSmartPointer<vtkDoubleArray> ComputeNek5000SpectralGradient(
    vtkDataSet* dataSet,
    vtkDataArray* pointArray,
    const std::string& outputName) {
    if (dataSet == nullptr || pointArray == nullptr ||
        pointArray->GetNumberOfTuples() != dataSet->GetNumberOfPoints()) {
        throw std::runtime_error(
            "Nek5000 spectral gradient requires a complete point-data array.");
    }
    const Nek5000SpectralMetadata meta = GetNek5000SpectralMetadata(dataSet);
    const int components = pointArray->GetNumberOfComponents();
    if (components <= 0) {
        throw std::runtime_error("Nek5000 spectral gradient input has no components.");
    }
    const std::vector<double> nodesX = GllNodes(meta.nx);
    const std::vector<double> nodesY = GllNodes(meta.ny);
    const std::vector<double> nodesZ = meta.is3D() ? GllNodes(meta.nz)
                                                   : std::vector<double>{0.0};
    const std::vector<double> dx = DerivativeMatrix(nodesX);
    const std::vector<double> dy = DerivativeMatrix(nodesY);
    const std::vector<double> dz = meta.is3D() ? DerivativeMatrix(nodesZ)
                                               : std::vector<double>();

    vtkPointSet* pointSet = vtkPointSet::SafeDownCast(dataSet);
    vtkDataArray* coordinates =
        pointSet == nullptr || pointSet->GetPoints() == nullptr
            ? nullptr
            : pointSet->GetPoints()->GetData();
    if (coordinates == nullptr || coordinates->GetNumberOfComponents() != 3) {
        throw std::runtime_error(
            "Nek5000 geometry must expose a three-component point-coordinate array.");
    }
    vtkDoubleArray* doubleCoordinates = vtkDoubleArray::SafeDownCast(coordinates);
    const double* coordinateValues =
        doubleCoordinates == nullptr ? nullptr : doubleCoordinates->GetPointer(0);
    const auto coordinate = [coordinates, coordinateValues](std::size_t pointIndex,
                                                            int component) {
        return coordinateValues == nullptr
                   ? coordinates->GetComponent(static_cast<vtkIdType>(pointIndex), component)
                   : coordinateValues[3 * pointIndex + static_cast<std::size_t>(component)];
    };

    vtkSmartPointer<vtkDoubleArray> output =
        vtkSmartPointer<vtkDoubleArray>::New();
    const std::string name = outputName.empty()
                                 ? std::string(pointArray->GetName() == nullptr
                                                   ? "SpectralGradient"
                                                   : pointArray->GetName()) +
                                       "_SpectralGradient"
                                 : outputName;
    output->SetName(name.c_str());
    output->SetNumberOfComponents(3 * components);
    output->SetNumberOfTuples(static_cast<vtkIdType>(meta.pointCount()));
    output->FillValue(0.0);

    const std::size_t nxyz = meta.pointsPerElement();
    std::vector<double> qr(static_cast<std::size_t>(components), 0.0);
    std::vector<double> qs(static_cast<std::size_t>(components), 0.0);
    std::vector<double> qt(static_cast<std::size_t>(components), 0.0);
    for (std::size_t element = 0; element < meta.elementCount; ++element) {
        const std::size_t base = element * nxyz;
        for (int k = 0; k < meta.nz; ++k) {
            for (int j = 0; j < meta.ny; ++j) {
                for (int i = 0; i < meta.nx; ++i) {
                    std::array<double, 3> xr{};
                    std::array<double, 3> xs{};
                    std::array<double, 3> xt{};
                    std::fill(qr.begin(), qr.end(), 0.0);
                    std::fill(qs.begin(), qs.end(), 0.0);
                    std::fill(qt.begin(), qt.end(), 0.0);
                    for (int m = 0; m < meta.nx; ++m) {
                        const double d = dx[static_cast<std::size_t>(i) *
                                                static_cast<std::size_t>(meta.nx) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, m, j, k);
                        for (int axis = 0; axis < 3; ++axis) {
                            xr[static_cast<std::size_t>(axis)] +=
                                d * coordinate(index, axis);
                        }
                        for (int component = 0; component < components; ++component) {
                            qr[static_cast<std::size_t>(component)] +=
                                d * pointArray->GetComponent(static_cast<vtkIdType>(index),
                                                             component);
                        }
                    }
                    for (int m = 0; m < meta.ny; ++m) {
                        const double d = dy[static_cast<std::size_t>(j) *
                                                static_cast<std::size_t>(meta.ny) +
                                            static_cast<std::size_t>(m)];
                        const std::size_t index =
                            base + NekLocalIndex(meta.nx, meta.ny, i, m, k);
                        for (int axis = 0; axis < 3; ++axis) {
                            xs[static_cast<std::size_t>(axis)] +=
                                d * coordinate(index, axis);
                        }
                        for (int component = 0; component < components; ++component) {
                            qs[static_cast<std::size_t>(component)] +=
                                d * pointArray->GetComponent(static_cast<vtkIdType>(index),
                                                             component);
                        }
                    }
                    if (meta.is3D()) {
                        for (int m = 0; m < meta.nz; ++m) {
                            const double d = dz[static_cast<std::size_t>(k) *
                                                    static_cast<std::size_t>(meta.nz) +
                                                static_cast<std::size_t>(m)];
                            const std::size_t index =
                                base + NekLocalIndex(meta.nx, meta.ny, i, j, m);
                            for (int axis = 0; axis < 3; ++axis) {
                                xt[static_cast<std::size_t>(axis)] +=
                                    d * coordinate(index, axis);
                            }
                            for (int component = 0; component < components; ++component) {
                                qt[static_cast<std::size_t>(component)] +=
                                    d * pointArray->GetComponent(
                                            static_cast<vtkIdType>(index), component);
                            }
                        }
                    }

                    const std::size_t outputIndex =
                        base + NekLocalIndex(meta.nx, meta.ny, i, j, k);
                    if (meta.is3D()) {
                        const double determinant =
                            xr[0] * (xs[1] * xt[2] - xt[1] * xs[2]) -
                            xs[0] * (xr[1] * xt[2] - xt[1] * xr[2]) +
                            xt[0] * (xr[1] * xs[2] - xs[1] * xr[2]);
                        if (!std::isfinite(determinant) ||
                            std::abs(determinant) <= 1.0e-30) {
                            throw std::runtime_error(
                                "Degenerate Nek5000 Jacobian while computing a spectral gradient.");
                        }
                        const double rX = (xs[1] * xt[2] - xt[1] * xs[2]) / determinant;
                        const double sX = (xt[1] * xr[2] - xr[1] * xt[2]) / determinant;
                        const double tX = (xr[1] * xs[2] - xs[1] * xr[2]) / determinant;
                        const double rY = (xt[0] * xs[2] - xs[0] * xt[2]) / determinant;
                        const double sY = (xr[0] * xt[2] - xt[0] * xr[2]) / determinant;
                        const double tY = (xs[0] * xr[2] - xr[0] * xs[2]) / determinant;
                        const double rZ = (xs[0] * xt[1] - xt[0] * xs[1]) / determinant;
                        const double sZ = (xt[0] * xr[1] - xr[0] * xt[1]) / determinant;
                        const double tZ = (xr[0] * xs[1] - xs[0] * xr[1]) / determinant;
                        for (int component = 0; component < components; ++component) {
                            const std::size_t c = static_cast<std::size_t>(component);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 0,
                                                 qr[c] * rX + qs[c] * sX + qt[c] * tX);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 1,
                                                 qr[c] * rY + qs[c] * sY + qt[c] * tY);
                            output->SetComponent(static_cast<vtkIdType>(outputIndex),
                                                 3 * component + 2,
                                                 qr[c] * rZ + qs[c] * sZ + qt[c] * tZ);
                        }
                    } else {
                        const double g00 = xr[0] * xr[0] + xr[1] * xr[1] + xr[2] * xr[2];
                        const double g01 = xr[0] * xs[0] + xr[1] * xs[1] + xr[2] * xs[2];
                        const double g11 = xs[0] * xs[0] + xs[1] * xs[1] + xs[2] * xs[2];
                        const double determinant = g00 * g11 - g01 * g01;
                        if (!std::isfinite(determinant) || determinant <= 1.0e-30) {
                            throw std::runtime_error(
                                "Degenerate 2-D Nek5000 metric while computing a spectral gradient.");
                        }
                        for (int component = 0; component < components; ++component) {
                            const std::size_t c = static_cast<std::size_t>(component);
                            const double coefficientR = (g11 * qr[c] - g01 * qs[c]) / determinant;
                            const double coefficientS = (-g01 * qr[c] + g00 * qs[c]) / determinant;
                            for (int axis = 0; axis < 3; ++axis) {
                                output->SetComponent(
                                    static_cast<vtkIdType>(outputIndex),
                                    3 * component + axis,
                                    coefficientR * xr[static_cast<std::size_t>(axis)] +
                                        coefficientS * xs[static_cast<std::size_t>(axis)]);
                            }
                        }
                    }
                }
            }
        }
    }
    return output;
}

}  // namespace streamcenterplus::io
