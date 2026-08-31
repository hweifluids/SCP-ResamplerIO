#include "streamcenterplus_resample/Core.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace streamcenterplus::resample {

namespace {

std::string Trim(const std::string& value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

std::string StripQuotes(std::string value) {
    value = Trim(value);
    if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                              (value.front() == '\'' && value.back() == '\''))) {
        return value.substr(1, value.size() - 2);
    }
    return value;
}

std::string StripInlineComment(const std::string& line) {
    bool inSingle = false;
    bool inDouble = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char ch = line[i];
        if (ch == '\'' && !inDouble) {
            inSingle = !inSingle;
        } else if (ch == '"' && !inSingle) {
            inDouble = !inDouble;
        } else if (ch == '#' && !inSingle && !inDouble) {
            return Trim(line.substr(0, i));
        }
    }
    return Trim(line);
}

std::string NormalizeKey(std::string key) {
    key = ToLower(Trim(key));
    std::replace(key.begin(), key.end(), '-', '_');
    return key;
}

std::vector<std::string> SplitList(std::string text) {
    text = Trim(text);
    if (text.size() >= 2 && ((text.front() == '[' && text.back() == ']') ||
                             (text.front() == '(' && text.back() == ')'))) {
        text = text.substr(1, text.size() - 2);
    }
    std::replace(text.begin(), text.end(), ';', ',');
    std::vector<std::string> values;
    if (text.find(',') == std::string::npos) {
        std::stringstream stream(text);
        std::string token;
        while (stream >> token) {
            token = StripQuotes(token);
            if (!token.empty()) {
                values.push_back(token);
            }
        }
    } else {
        std::stringstream stream(text);
        std::string token;
        while (std::getline(stream, token, ',')) {
            token = StripQuotes(token);
            if (!token.empty()) {
                values.push_back(token);
            }
        }
    }
    return values;
}

std::map<std::string, std::string> ParseDeck(const std::filesystem::path& deckFile) {
    std::ifstream input(deckFile);
    if (!input) {
        throw std::runtime_error("Cannot open resample .resample.scpo file: " + deckFile.string());
    }

    std::map<std::string, std::string> raw;
    std::vector<std::string> blockStack;
    std::string line;
    while (std::getline(input, line)) {
        line = StripInlineComment(line);
        if (line.empty()) {
            continue;
        }

        while (!line.empty() && line.front() == '}') {
            if (!blockStack.empty()) {
                blockStack.pop_back();
            }
            line = Trim(line.substr(1));
        }
        if (line.empty()) {
            continue;
        }

        if (line.back() == '{') {
            const std::string blockName = NormalizeKey(line.substr(0, line.size() - 1));
            if (!blockName.empty()) {
                blockStack.push_back(blockName);
            }
            continue;
        }

        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            continue;
        }

        std::string key = NormalizeKey(line.substr(0, eq));
        std::string value = Trim(line.substr(eq + 1));
        int bracketBalance = static_cast<int>(std::count(value.begin(), value.end(), '[')) -
                             static_cast<int>(std::count(value.begin(), value.end(), ']'));
        while (bracketBalance > 0 && std::getline(input, line)) {
            line = StripInlineComment(line);
            value += " " + line;
            bracketBalance += static_cast<int>(std::count(line.begin(), line.end(), '['));
            bracketBalance -= static_cast<int>(std::count(line.begin(), line.end(), ']'));
        }

        std::string fullKey;
        for (const std::string& block : blockStack) {
            fullKey += block + ".";
        }
        fullKey += key;
        raw[NormalizeKey(fullKey)] = value;
    }
    return raw;
}

bool TryRaw(const std::map<std::string, std::string>& raw, const std::string& key, std::string* out) {
    const std::string normalized = NormalizeKey(key);
    if (const auto it = raw.find(normalized); it != raw.end()) {
        if (out != nullptr) {
            *out = it->second;
        }
        return true;
    }
    if (const auto it = raw.find("resample." + normalized); it != raw.end()) {
        if (out != nullptr) {
            *out = it->second;
        }
        return true;
    }
    return false;
}

bool TryAliases(const std::map<std::string, std::string>& raw,
                const std::initializer_list<const char*> keys,
                std::string* out) {
    for (const char* key : keys) {
        if (TryRaw(raw, key, out)) {
            return true;
        }
    }
    return false;
}

bool TryString(const std::map<std::string, std::string>& raw,
               const std::initializer_list<const char*> keys,
               std::string* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    *out = StripQuotes(value);
    return true;
}

bool TryBoolText(const std::string& rawValue, bool* out) {
    if (out == nullptr) {
        return false;
    }
    const std::string value = ToLower(StripQuotes(rawValue));
    if (value == "1" || value == "true" || value == "yes" || value == "on" ||
        value == "enabled" || value == "enable") {
        *out = true;
        return true;
    }
    if (value == "0" || value == "false" || value == "no" || value == "off" ||
        value == "disabled" || value == "disable") {
        *out = false;
        return true;
    }
    return false;
}

bool TryBool(const std::map<std::string, std::string>& raw,
             const std::initializer_list<const char*> keys,
             bool* out) {
    std::string value;
    return TryAliases(raw, keys, &value) && TryBoolText(value, out);
}

bool TryInt(const std::map<std::string, std::string>& raw,
            const std::initializer_list<const char*> keys,
            int* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    try {
        *out = std::stoi(StripQuotes(value));
        return true;
    } catch (...) {
        return false;
    }
}

bool TryDouble(const std::map<std::string, std::string>& raw,
               const std::initializer_list<const char*> keys,
               double* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    value = ToLower(StripQuotes(value));
    if (value == "auto" || value == "default" || value == "none" || value.empty()) {
        *out = 0.0;
        return true;
    }
    try {
        *out = std::stod(value);
        return true;
    } catch (...) {
        return false;
    }
}

bool TryVector3i(const std::map<std::string, std::string>& raw,
                 const std::initializer_list<const char*> keys,
                 std::array<int, 3>* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    const std::vector<std::string> tokens = SplitList(value);
    if (tokens.size() != 3) {
        return false;
    }
    try {
        *out = {std::stoi(tokens[0]), std::stoi(tokens[1]), std::stoi(tokens[2])};
        return true;
    } catch (...) {
        return false;
    }
}

void ParseClusterJobHints(const std::map<std::string, std::string>& raw,
                          Parameters* params) {
    if (params == nullptr) {
        return;
    }
    auto& hints = params->cluster;
    TryBool(raw,
            {"cluster.package", "cluster.package_requested", "cluster.enabled",
             "hpc.package", "hpc.enabled"},
            &hints.packageRequested);
    TryString(raw,
              {"cluster.package_path", "cluster.output", "cluster.package_file",
               "hpc.package_path", "hpc.output", "job_package"},
              &hints.packagePath);
    TryString(raw, {"cluster.payload_mode", "hpc.payload_mode"}, &hints.payloadMode);
    TryString(raw, {"cluster.control_format", "hpc.control_format"}, &hints.controlFormat);
    TryString(raw,
              {"cluster.partition_mode", "hpc.partition_mode", "partition_mode"},
              &hints.partitionMode);
    TryString(raw, {"cluster.scheduler", "hpc.scheduler"}, &hints.scheduler);
    TryString(raw, {"cluster.queue", "hpc.queue"}, &hints.queue);
    TryString(raw, {"cluster.walltime", "hpc.walltime"}, &hints.walltime);
    TryInt(raw, {"cluster.nodes", "hpc.nodes"}, &hints.nodes);
    TryInt(raw, {"cluster.ranks", "hpc.ranks", "mpi_ranks"}, &hints.ranks);
    TryInt(raw, {"cluster.gpus", "hpc.gpus"}, &hints.gpus);
    TryInt(raw, {"cluster.gpus_per_node", "hpc.gpus_per_node"}, &hints.gpusPerNode);
    TryInt(raw,
           {"cluster.cpu_threads_per_rank", "hpc.cpu_threads_per_rank"},
           &hints.cpuThreadsPerRank);
    TryDouble(raw, {"cluster.memory_gb", "hpc.memory_gb"}, &hints.memoryGb);
    if (!hints.packagePath.empty()) {
        hints.packageRequested = true;
    }
}

bool IsAutoText(const std::string& rawValue) {
    const std::string value = ToLower(StripQuotes(rawValue));
    return value.empty() || value == "all" || value == "auto" || value == "default" ||
           value == "none";
}

bool IsAutoLowerBound(const std::string& rawValue) {
    const std::string value = ToLower(StripQuotes(rawValue));
    return value == "all" || value == "auto" || value == "default" || value == "none" ||
           value == "min" || value == "minimum" || value == "start" || value == "first";
}

bool IsAutoUpperBound(const std::string& rawValue) {
    const std::string value = ToLower(StripQuotes(rawValue));
    return value == "all" || value == "auto" || value == "default" || value == "none" ||
           value == "max" || value == "maximum" || value == "end" || value == "last";
}

bool TryAxisRange(const std::map<std::string, std::string>& raw,
                  const std::initializer_list<const char*> keys,
                  AxisRange* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    if (IsAutoText(value)) {
        *out = AxisRange{};
        return true;
    }
    const std::vector<std::string> tokens = SplitList(value);
    if (tokens.size() != 2) {
        throw std::runtime_error("Coordinate crop ranges must contain two values, or be empty/all/auto.");
    }

    AxisRange parsed;
    if (!IsAutoLowerBound(tokens[0])) {
        parsed.min = std::stod(tokens[0]);
    }
    if (!IsAutoUpperBound(tokens[1])) {
        parsed.max = std::stod(tokens[1]);
    }
    if (parsed.min.has_value() && parsed.max.has_value() && *parsed.max < *parsed.min) {
        throw std::runtime_error("Coordinate crop range max must be >= min.");
    }
    *out = parsed;
    return true;
}

bool TryFrameRange(const std::map<std::string, std::string>& raw,
                   const std::initializer_list<const char*> keys,
                   FrameRange* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    if (IsAutoText(value)) {
        *out = FrameRange{};
        return true;
    }
    const std::vector<std::string> tokens = SplitList(value);
    if (tokens.size() != 2) {
        throw std::runtime_error("Frame range must contain two frame indices, or be empty/all/auto.");
    }

    FrameRange parsed;
    if (!IsAutoLowerBound(tokens[0])) {
        parsed.start = static_cast<int>(std::nearbyint(std::stod(tokens[0])));
    }
    if (!IsAutoUpperBound(tokens[1])) {
        parsed.end = static_cast<int>(std::nearbyint(std::stod(tokens[1])));
    }
    if (parsed.start.has_value() && parsed.end.has_value() && *parsed.end < *parsed.start) {
        throw std::runtime_error("Frame range end must be >= start.");
    }
    *out = parsed;
    return true;
}

bool TryVector3dOrScalar(const std::map<std::string, std::string>& raw,
                         const std::initializer_list<const char*> keys,
                         std::array<double, 3>* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    if (IsAutoText(value)) {
        *out = {1.0, 1.0, 1.0};
        return true;
    }
    const std::vector<std::string> tokens = SplitList(value);
    if (tokens.size() == 1) {
        const double scale = std::stod(tokens[0]);
        *out = {scale, scale, scale};
        return true;
    }
    if (tokens.size() == 3) {
        *out = {std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2])};
        return true;
    }
    throw std::runtime_error("Input scale must be a scalar or a three-value vector, or be empty/all/auto.");
}

Backend ParseBackend(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "cpu" || value == "host") {
        return Backend::Cpu;
    }
    if (value == "gpu" || value == "cuda" || value == "device") {
        return Backend::Gpu;
    }
    throw std::runtime_error("Unsupported resample backend: " + rawValue);
}

InterpolationMethod ParseInterpolationMethod(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "linear" || value == "probe" || value == "trilinear") {
        return InterpolationMethod::Linear;
    }
    if (value == "idw" || value == "inverse_distance" ||
        value == "inverse_distance_weighting" || value == "shepard") {
        return InterpolationMethod::InverseDistanceWeighting;
    }
    if (value == "kriging" || value == "ordinary_kriging") {
        return InterpolationMethod::Kriging;
    }
    if (value == "vtk" || value == "vtk_resample" ||
        value == "vtk_resample_to_image" || value == "resample_to_image") {
        return InterpolationMethod::VtkResampleToImage;
    }
    throw std::runtime_error("Unsupported resample interpolation method: " + rawValue);
}

NeighborMode ParseNeighborMode(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "radius" || value == "within_radius") {
        return NeighborMode::Radius;
    }
    if (value == "nearest" || value == "n_closest" || value == "nclosest" ||
        value == "neighbors") {
        return NeighborMode::Nearest;
    }
    throw std::runtime_error("Unsupported neighbor mode: " + rawValue);
}

NullStrategy ParseNullStrategy(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "mask" || value == "mask_points" || value == "hidden") {
        return NullStrategy::Mask;
    }
    if (value == "closest" || value == "closest_point") {
        return NullStrategy::Closest;
    }
    if (value == "null" || value == "null_value" || value == "value") {
        return NullStrategy::NullValue;
    }
    if (value == "error" || value == "fail") {
        return NullStrategy::Error;
    }
    throw std::runtime_error("Unsupported null strategy: " + rawValue);
}

KrigingVariogram ParseKrigingVariogram(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "spherical") {
        return KrigingVariogram::Spherical;
    }
    if (value == "exponential") {
        return KrigingVariogram::Exponential;
    }
    if (value == "gaussian") {
        return KrigingVariogram::Gaussian;
    }
    if (value == "linear") {
        return KrigingVariogram::Linear;
    }
    throw std::runtime_error("Unsupported kriging variogram model: " + rawValue);
}

KrigingFallback ParseKrigingFallback(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "idw" || value == "inverse_distance_weighting") {
        return KrigingFallback::Idw;
    }
    if (value == "closest" || value == "closest_point") {
        return KrigingFallback::Closest;
    }
    if (value == "error" || value == "fail") {
        return KrigingFallback::Error;
    }
    throw std::runtime_error("Unsupported kriging fallback: " + rawValue);
}

BoundsMode ParseBoundsMode(const std::string& rawValue) {
    std::string value = ToLower(StripQuotes(rawValue));
    std::replace(value.begin(), value.end(), '-', '_');
    if (value == "global" || value == "series" || value == "time_series") {
        return BoundsMode::Global;
    }
    if (value == "per_frame" || value == "frame" || value == "local") {
        return BoundsMode::PerFrame;
    }
    throw std::runtime_error("Unsupported bounds mode: " + rawValue);
}

std::string DefaultBasename(const std::filesystem::path& deckFile) {
    std::string name = deckFile.filename().string();
    const std::string suffix = ".resample.scpo";
    const std::string lower = ToLower(name);
    if (lower.size() >= suffix.size() &&
        lower.compare(lower.size() - suffix.size(), suffix.size(), suffix) == 0) {
        name.resize(name.size() - suffix.size());
    } else {
        name = deckFile.stem().string();
    }
    return name.empty() ? std::string("resampled") : name;
}

void ParseEnums(const std::map<std::string, std::string>& raw, Parameters* params) {
    std::string value;
    if (TryString(raw, {"backend", "compute", "compute_backend", "implementation"}, &value)) {
        params->backend = ParseBackend(value);
    }
    if (TryString(raw,
                  {"interpolation.method", "interpolation_method", "method", "algorithm"},
                  &value)) {
        params->method = ParseInterpolationMethod(value);
    }
    if (TryString(raw, {"bounds_mode", "output.bounds_mode"}, &value)) {
        params->boundsMode = ParseBoundsMode(value);
    }
    if (TryString(raw, {"interpolation.idw.neighbor_mode", "idw.neighbor_mode"}, &value)) {
        params->idw.neighborMode = ParseNeighborMode(value);
    }
    if (TryString(raw, {"interpolation.idw.null_strategy", "idw.null_strategy"}, &value)) {
        params->idw.nullStrategy = ParseNullStrategy(value);
    }
    if (TryString(raw, {"interpolation.kriging.neighbor_mode", "kriging.neighbor_mode"}, &value)) {
        params->kriging.neighborMode = ParseNeighborMode(value);
    }
    if (TryString(raw,
                  {"interpolation.kriging.variogram", "interpolation.kriging.variogram_model",
                   "kriging.variogram", "kriging.variogram_model"},
                  &value)) {
        params->kriging.variogram = ParseKrigingVariogram(value);
    }
    if (TryString(raw, {"interpolation.kriging.fallback", "kriging.fallback"}, &value)) {
        params->kriging.fallback = ParseKrigingFallback(value);
    }
}

void Validate(Parameters* params) {
    if (params == nullptr) {
        return;
    }
    const std::string status = ToLower(params->solverStatus);
    if (status != "solved" && status != "unsolved" && status != "error") {
        throw std::runtime_error("solver_status must be solved, unsolved, or error.");
    }
    params->solverStatus = status;
    if (params->inputPath.empty()) {
        throw std::runtime_error("Missing required input path. Use input_path or vtk_folder.");
    }
    if (params->outputDir.empty()) {
        params->outputDir = std::filesystem::path("./results/resample_results");
    }
    for (int dim : params->dimensions) {
        if (dim < 1) {
            throw std::runtime_error("resample dimensions/Nx Ny Nz must all be >= 1.");
        }
    }
    if (params->physicalDt <= 0.0 || !std::isfinite(params->physicalDt)) {
        throw std::runtime_error("physical_dt/dt must be positive when provided.");
    }
    if (params->cluster.nodes < 0 || params->cluster.ranks < 0 ||
        params->cluster.gpus < 0 || params->cluster.gpusPerNode < 0 ||
        params->cluster.cpuThreadsPerRank < 0 || params->cluster.memoryGb < 0.0) {
        throw std::runtime_error("cluster resource hints must be non-negative.");
    }
    if (params->frameRange.start.has_value() && *params->frameRange.start < 0) {
        throw std::runtime_error("frame_range/T_range start must be non-negative.");
    }
    if (params->frameRange.start.has_value() && params->frameRange.end.has_value() &&
        *params->frameRange.end < *params->frameRange.start) {
        throw std::runtime_error("frame_range/T_range end must be >= start.");
    }
    for (const AxisRange& range : params->cropRanges) {
        if ((range.min.has_value() && !std::isfinite(*range.min)) ||
            (range.max.has_value() && !std::isfinite(*range.max))) {
            throw std::runtime_error("Coordinate crop ranges must be finite when provided.");
        }
        if (range.min.has_value() && range.max.has_value() && *range.max < *range.min) {
            throw std::runtime_error("Coordinate crop range max must be >= min.");
        }
    }
    for (double scale : params->inputScale) {
        if (scale <= 0.0 || !std::isfinite(scale)) {
            throw std::runtime_error("Input scale values must be positive finite numbers.");
        }
    }
    params->idw.neighbors = std::max(1, params->idw.neighbors);
    params->kriging.neighbors = std::max(1, params->kriging.neighbors);
    params->kriging.minNeighbors = std::max(1, params->kriging.minNeighbors);
    if (params->kriging.minNeighbors > params->kriging.neighbors &&
        params->kriging.neighborMode == NeighborMode::Nearest) {
        params->kriging.minNeighbors = params->kriging.neighbors;
    }
    if (params->idw.power <= 0.0 || !std::isfinite(params->idw.power)) {
        throw std::runtime_error("IDW power must be positive.");
    }
    if (params->idw.radius < 0.0 || params->kriging.radius < 0.0 ||
        params->kriging.range < 0.0) {
        throw std::runtime_error("Interpolation radii/range must be non-negative; use 0/auto for automatic values.");
    }
    if (params->kriging.sill <= 0.0 || !std::isfinite(params->kriging.sill)) {
        throw std::runtime_error("Kriging sill must be positive.");
    }
    if (params->kriging.nugget < 0.0 || params->kriging.regularization < 0.0) {
        throw std::runtime_error("Kriging nugget and regularization must be non-negative.");
    }
    if (params->kriging.fallbackIdwPower <= 0.0 || !std::isfinite(params->kriging.fallbackIdwPower)) {
        throw std::runtime_error("Kriging fallback_idw_power must be positive.");
    }
    if (params->outputBasename.empty()) {
        params->outputBasename = DefaultBasename(params->deckFile);
    }
}

std::string FormatBool(bool value) {
    return value ? "Enabled" : "Disabled";
}

std::string FormatDouble(double value) {
    std::ostringstream out;
    out << std::setprecision(15) << value;
    return out.str();
}

std::string FormatAxisRange(const AxisRange& range) {
    if (!range.min.has_value() && !range.max.has_value()) {
        return "auto";
    }
    return "[" + (range.min.has_value() ? FormatDouble(*range.min) : std::string("auto")) + ", " +
           (range.max.has_value() ? FormatDouble(*range.max) : std::string("auto")) + "]";
}

std::string FormatFrameRange(const FrameRange& range) {
    if (!range.start.has_value() && !range.end.has_value()) {
        return "auto";
    }
    return "[" + (range.start.has_value() ? std::to_string(*range.start) : std::string("auto")) + ", " +
           (range.end.has_value() ? std::to_string(*range.end) : std::string("auto")) + "]";
}

std::string FormatVector3d(const std::array<double, 3>& values) {
    return "[" + FormatDouble(values[0]) + ", " + FormatDouble(values[1]) + ", " +
           FormatDouble(values[2]) + "]";
}

}  // namespace

std::string ToString(Backend value) {
    switch (value) {
    case Backend::Cpu:
        return "CPU";
    case Backend::Gpu:
        return "GPU";
    }
    return "unknown";
}

std::string ToString(InterpolationMethod value) {
    switch (value) {
    case InterpolationMethod::Linear:
        return "linear";
    case InterpolationMethod::InverseDistanceWeighting:
        return "inverse_distance_weighting";
    case InterpolationMethod::Kriging:
        return "kriging";
    case InterpolationMethod::VtkResampleToImage:
        return "vtk_resample_to_image";
    }
    return "unknown";
}

std::string ToString(NeighborMode value) {
    switch (value) {
    case NeighborMode::Radius:
        return "radius";
    case NeighborMode::Nearest:
        return "nearest";
    }
    return "unknown";
}

std::string ToString(NullStrategy value) {
    switch (value) {
    case NullStrategy::Mask:
        return "mask";
    case NullStrategy::Closest:
        return "closest";
    case NullStrategy::NullValue:
        return "null_value";
    case NullStrategy::Error:
        return "error";
    }
    return "unknown";
}

std::string ToString(KrigingVariogram value) {
    switch (value) {
    case KrigingVariogram::Spherical:
        return "spherical";
    case KrigingVariogram::Exponential:
        return "exponential";
    case KrigingVariogram::Gaussian:
        return "gaussian";
    case KrigingVariogram::Linear:
        return "linear";
    }
    return "unknown";
}

std::string ToString(KrigingFallback value) {
    switch (value) {
    case KrigingFallback::Idw:
        return "idw";
    case KrigingFallback::Closest:
        return "closest";
    case KrigingFallback::Error:
        return "error";
    }
    return "unknown";
}

std::string ToString(BoundsMode value) {
    switch (value) {
    case BoundsMode::Global:
        return "global";
    case BoundsMode::PerFrame:
        return "per_frame";
    }
    return "unknown";
}

std::string ToString(DataSetFrameKind value) {
    switch (value) {
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
    return "unknown";
}

Parameters LoadParameters(const std::filesystem::path& deckFile) {
    if (!std::filesystem::exists(deckFile)) {
        throw std::runtime_error("Missing resample .resample.scpo parameter file: " + deckFile.string());
    }

    Parameters params;
    params.deckFile = deckFile;
    params.outputBasename = DefaultBasename(deckFile);

    const std::map<std::string, std::string> raw = ParseDeck(deckFile);
    std::cout << "Resample parameter file '" << std::filesystem::absolute(deckFile).string()
              << "' has been successfully loaded.\n";

    TryString(raw, {"solver_status", "resample.solver_status"}, &params.solverStatus);
    TryBool(raw, {"enabled", "resample.enabled"}, &params.enabled);
    if (std::string pathValue; TryString(raw,
                                         {"input_path", "input", "vtk_folder", "vtk_path", "source"},
                                         &pathValue)) {
        params.inputPath = pathValue;
    }
    if (std::string pathValue; TryString(raw, {"output_dir", "output_folder", "output.path"}, &pathValue)) {
        params.outputDir = pathValue;
    }
    TryString(raw, {"output_basename", "basename", "output.basename"}, &params.outputBasename);
    TryVector3i(raw, {"dimensions", "dims", "resolution", "grid_size", "sampling_dimensions"},
                &params.dimensions);
    TryDouble(raw, {"physical_dt", "dt", "time_step"}, &params.physicalDt);
    TryFrameRange(raw, {"frame_range", "t_range", "time_range", "input.frame_range"},
                  &params.frameRange);
    TryAxisRange(raw,
                 {"crop.x_range", "crop_x_range", "x_range", "xrange", "domain.x_range",
                  "bounds.x_range"},
                 &params.cropRanges[0]);
    TryAxisRange(raw,
                 {"crop.y_range", "crop_y_range", "y_range", "yrange", "domain.y_range",
                  "bounds.y_range"},
                 &params.cropRanges[1]);
    TryAxisRange(raw,
                 {"crop.z_range", "crop_z_range", "z_range", "zrange", "domain.z_range",
                  "bounds.z_range"},
                 &params.cropRanges[2]);
    TryVector3dOrScalar(raw,
                        {"scale", "scale_ratio", "input_scale", "input.scale",
                         "coordinate_scale", "coordinates.scale"},
                        &params.inputScale);

    TryBool(raw,
            {"interpolation.linear.compute_tolerance", "linear.compute_tolerance"},
            &params.linear.computeTolerance);
    TryDouble(raw, {"interpolation.linear.tolerance", "linear.tolerance"}, &params.linear.tolerance);
    TryBool(raw,
            {"interpolation.linear.snap_to_cell_with_closest_point", "linear.snap_to_cell_with_closest_point"},
            &params.linear.snapToCellWithClosestPoint);

    TryInt(raw, {"interpolation.idw.neighbors", "idw.neighbors", "idw.k"}, &params.idw.neighbors);
    TryDouble(raw, {"interpolation.idw.radius", "idw.radius"}, &params.idw.radius);
    TryDouble(raw, {"interpolation.idw.power", "idw.power"}, &params.idw.power);
    TryDouble(raw, {"interpolation.idw.null_value", "idw.null_value"}, &params.idw.nullValue);

    TryInt(raw,
           {"interpolation.kriging.neighbors", "kriging.neighbors", "kriging.k"},
           &params.kriging.neighbors);
    TryInt(raw,
           {"interpolation.kriging.min_neighbors", "kriging.min_neighbors"},
           &params.kriging.minNeighbors);
    TryDouble(raw, {"interpolation.kriging.radius", "kriging.radius"}, &params.kriging.radius);
    TryDouble(raw, {"interpolation.kriging.range", "kriging.range"}, &params.kriging.range);
    TryDouble(raw, {"interpolation.kriging.sill", "kriging.sill"}, &params.kriging.sill);
    TryDouble(raw, {"interpolation.kriging.nugget", "kriging.nugget"}, &params.kriging.nugget);
    TryDouble(raw,
              {"interpolation.kriging.regularization", "kriging.regularization"},
              &params.kriging.regularization);
    TryDouble(raw, {"interpolation.kriging.null_value", "kriging.null_value"}, &params.kriging.nullValue);
    TryDouble(raw,
              {"interpolation.kriging.fallback_idw_power", "kriging.fallback_idw_power"},
              &params.kriging.fallbackIdwPower);

    ParseEnums(raw, &params);
    ParseClusterJobHints(raw, &params);
    Validate(&params);
    return params;
}

std::string PrettyPrint(const Parameters& params) {
    std::ostringstream out;
    out << "\n=== File IO ===\n";
    out << "Input path:             " << params.inputPath.string() << "\n";
    out << "Output directory:       " << params.outputDir.string() << "\n";
    out << "Output basename:        " << params.outputBasename << "\n";
    out << "Physical dt:            " << FormatDouble(params.physicalDt) << "\n";

    out << "\n=== Resample Grid ===\n";
    out << "Dimensions:             [" << params.dimensions[0] << ", "
        << params.dimensions[1] << ", " << params.dimensions[2] << "]\n";
    out << "Bounds mode:            " << ToString(params.boundsMode) << "\n";
    out << "Frame range:            " << FormatFrameRange(params.frameRange) << "\n";
    out << "X crop range:           " << FormatAxisRange(params.cropRanges[0]) << "\n";
    out << "Y crop range:           " << FormatAxisRange(params.cropRanges[1]) << "\n";
    out << "Z crop range:           " << FormatAxisRange(params.cropRanges[2]) << "\n";
    out << "Input scale:            " << FormatVector3d(params.inputScale) << "\n";

    out << "\n=== Compute ===\n";
    out << "Enabled:                " << FormatBool(params.enabled) << "\n";
    out << "Backend:                " << ToString(params.backend) << "\n";
    out << "Interpolation:          " << ToString(params.method) << "\n";

    out << "\n=== Linear ===\n";
    out << "Compute tolerance:      " << FormatBool(params.linear.computeTolerance) << "\n";
    out << "Tolerance:              " << FormatDouble(params.linear.tolerance) << "\n";
    out << "Snap closest cell:      " << FormatBool(params.linear.snapToCellWithClosestPoint) << "\n";

    out << "\n=== IDW ===\n";
    out << "Neighbor mode:          " << ToString(params.idw.neighborMode) << "\n";
    out << "Neighbors:              " << params.idw.neighbors << "\n";
    out << "Radius:                 " << (params.idw.radius > 0.0 ? FormatDouble(params.idw.radius) : "auto") << "\n";
    out << "Power:                  " << FormatDouble(params.idw.power) << "\n";
    out << "Null strategy:          " << ToString(params.idw.nullStrategy) << "\n";
    out << "Null value:             " << FormatDouble(params.idw.nullValue) << "\n";

    out << "\n=== Kriging ===\n";
    out << "Variogram:              " << ToString(params.kriging.variogram) << "\n";
    out << "Neighbor mode:          " << ToString(params.kriging.neighborMode) << "\n";
    out << "Neighbors:              " << params.kriging.neighbors << "\n";
    out << "Min neighbors:          " << params.kriging.minNeighbors << "\n";
    out << "Radius:                 " << (params.kriging.radius > 0.0 ? FormatDouble(params.kriging.radius) : "auto") << "\n";
    out << "Range:                  " << (params.kriging.range > 0.0 ? FormatDouble(params.kriging.range) : "auto") << "\n";
    out << "Sill:                   " << FormatDouble(params.kriging.sill) << "\n";
    out << "Nugget:                 " << FormatDouble(params.kriging.nugget) << "\n";
    out << "Regularization:         " << FormatDouble(params.kriging.regularization) << "\n";
    out << "Fallback:               " << ToString(params.kriging.fallback) << "\n";
    out << "Fallback IDW power:     " << FormatDouble(params.kriging.fallbackIdwPower) << "\n";
    if (streamcenterplus::HasClusterJobHints(params.cluster)) {
        out << "\n=== Cluster Package Hints ===\n";
        out << "Package requested:      " << FormatBool(params.cluster.packageRequested) << "\n";
        out << "Package path:           "
            << (params.cluster.packagePath.empty() ? "auto" : params.cluster.packagePath) << "\n";
        out << "Payload mode:           " << params.cluster.payloadMode << "\n";
        out << "Control format:         " << params.cluster.controlFormat << "\n";
        out << "Partition mode:         " << params.cluster.partitionMode << "\n";
        out << "Scheduler:              " << params.cluster.scheduler << "\n";
        out << "Queue:                  "
            << (params.cluster.queue.empty() ? "auto" : params.cluster.queue) << "\n";
        out << "Walltime:               "
            << (params.cluster.walltime.empty() ? "auto" : params.cluster.walltime) << "\n";
        out << "Nodes:                  " << params.cluster.nodes << "\n";
        out << "MPI ranks:              " << params.cluster.ranks << "\n";
        out << "GPUs:                   " << params.cluster.gpus << "\n";
        out << "GPUs per node:          " << params.cluster.gpusPerNode << "\n";
        out << "CPU threads/rank:       " << params.cluster.cpuThreadsPerRank << "\n";
        out << "Memory GB:              " << FormatDouble(params.cluster.memoryGb) << "\n";
    }
    return out.str();
}

}  // namespace streamcenterplus::resample
