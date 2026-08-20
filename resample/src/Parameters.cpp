#include "streamcenterplus_resample/Core.h"
#include "ParameterValidation.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
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
    const bool startsQuote = !value.empty() && (value.front() == '"' || value.front() == '\'');
    const bool endsQuote = !value.empty() && (value.back() == '"' || value.back() == '\'');
    if (startsQuote || endsQuote) {
        if (value.size() < 2 || value.front() != value.back()) {
            throw std::runtime_error("Parameter value has unmatched or mismatched quotes: " + value);
        }
    }
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
    std::size_t lineNumber = 0;
    while (std::getline(input, line)) {
        ++lineNumber;
        line = StripInlineComment(line);
        if (line.empty()) {
            continue;
        }

        while (!line.empty() && line.front() == '}') {
            if (blockStack.empty()) {
                throw std::runtime_error("Unexpected closing block at parameter line " +
                                         std::to_string(lineNumber) + ".");
            }
            blockStack.pop_back();
            line = Trim(line.substr(1));
        }
        if (line.empty()) {
            continue;
        }

        if (line.back() == '{') {
            const std::string blockName = NormalizeKey(line.substr(0, line.size() - 1));
            if (blockName.empty()) {
                throw std::runtime_error("Empty parameter block name at line " +
                                         std::to_string(lineNumber) + ".");
            }
            blockStack.push_back(blockName);
            continue;
        }

        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::runtime_error("Malformed parameter line " + std::to_string(lineNumber) +
                                     ": expected key=value or a block delimiter.");
        }

        std::string key = NormalizeKey(line.substr(0, eq));
        std::string value = Trim(line.substr(eq + 1));
        int bracketBalance = static_cast<int>(std::count(value.begin(), value.end(), '[')) -
                             static_cast<int>(std::count(value.begin(), value.end(), ']'));
        while (bracketBalance > 0 && std::getline(input, line)) {
            ++lineNumber;
            line = StripInlineComment(line);
            value += " " + line;
            bracketBalance += static_cast<int>(std::count(line.begin(), line.end(), '['));
            bracketBalance -= static_cast<int>(std::count(line.begin(), line.end(), ']'));
        }
        if (bracketBalance != 0) {
            throw std::runtime_error("Unbalanced list brackets near parameter line " +
                                     std::to_string(lineNumber) + ".");
        }

        std::string fullKey;
        for (const std::string& block : blockStack) {
            fullKey += block + ".";
        }
        fullKey += key;
        fullKey = NormalizeKey(fullKey);
        if (key.empty()) {
            throw std::runtime_error("Parameter line " + std::to_string(lineNumber) +
                                     " has an empty key.");
        }
        if (!raw.emplace(fullKey, value).second) {
            throw std::runtime_error("Duplicate parameter key '" + fullKey + "'.");
        }
    }
    if (!blockStack.empty()) {
        throw std::runtime_error("Parameter file ended before closing block '" +
                                 blockStack.back() + "'.");
    }
    return raw;
}

bool TryAliases(const std::map<std::string, std::string>& raw,
                const std::initializer_list<const char*> keys,
                std::string* out) {
    std::set<std::string> visited;
    std::optional<std::pair<std::string, std::string>> selected;
    for (const char* key : keys) {
        const std::string normalized = NormalizeKey(key);
        const std::array<std::string, 2> candidates = {
            normalized,
            "resample." + normalized,
        };
        for (const std::string& candidate : candidates) {
            if (!visited.insert(candidate).second) {
                continue;
            }
            const auto found = raw.find(candidate);
            if (found == raw.end()) {
                continue;
            }
            if (selected.has_value()) {
                throw std::runtime_error("Parameter aliases '" + selected->first + "' and '" +
                                         candidate + "' were both provided.");
            }
            selected = std::make_pair(candidate, found->second);
        }
    }
    if (!selected.has_value()) {
        return false;
    }
    if (out != nullptr) {
        *out = selected->second;
    }
    return true;
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
    if (!TryAliases(raw, keys, &value)) {
        return false;
    }
    if (!TryBoolText(value, out)) {
        throw std::runtime_error("Invalid boolean parameter value: '" + value + "'.");
    }
    return true;
}

int ParseIntToken(const std::string& rawValue, const char* description) {
    const std::string value = StripQuotes(rawValue);
    if (value.empty()) {
        throw std::runtime_error(std::string(description) + " must be an integer, not an empty value.");
    }
    try {
        std::size_t consumed = 0;
        const long long parsed = std::stoll(value, &consumed, 10);
        if (consumed != value.size() ||
            parsed < static_cast<long long>(std::numeric_limits<int>::min()) ||
            parsed > static_cast<long long>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("out of range or trailing characters");
        }
        return static_cast<int>(parsed);
    } catch (const std::exception&) {
        throw std::runtime_error(std::string(description) +
                                 " must be a complete in-range integer token; received '" +
                                 rawValue + "'.");
    }
}

double ParseFiniteDoubleToken(const std::string& rawValue, const char* description) {
    const std::string value = StripQuotes(rawValue);
    if (value.empty()) {
        throw std::runtime_error(std::string(description) + " must be a finite number, not an empty value.");
    }
    try {
        std::size_t consumed = 0;
        const double parsed = std::stod(value, &consumed);
        if (consumed != value.size() || !std::isfinite(parsed)) {
            throw std::runtime_error("non-finite or trailing characters");
        }
        return parsed;
    } catch (const std::exception&) {
        throw std::runtime_error(std::string(description) +
                                 " must be a complete finite numeric token; received '" +
                                 rawValue + "'.");
    }
}

bool TryInt(const std::map<std::string, std::string>& raw,
            const std::initializer_list<const char*> keys,
            int* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    *out = ParseIntToken(value, "Integer parameter");
    return true;
}

bool TryDouble(const std::map<std::string, std::string>& raw,
               const std::initializer_list<const char*> keys,
               double* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    *out = ParseFiniteDoubleToken(value, "Numeric parameter");
    return true;
}

bool TryAutoDouble(const std::map<std::string, std::string>& raw,
                   const std::initializer_list<const char*> keys,
                   double* out) {
    std::string value;
    if (!TryAliases(raw, keys, &value) || out == nullptr) {
        return false;
    }
    const std::string normalized = ToLower(StripQuotes(value));
    if (normalized == "auto" || normalized == "default" || normalized == "none") {
        *out = 0.0;
        return true;
    }
    *out = ParseFiniteDoubleToken(value, "Numeric parameter");
    return true;
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
        throw std::runtime_error("Grid dimensions must contain exactly three integer values.");
    }
    *out = {ParseIntToken(tokens[0], "Grid dimension"),
            ParseIntToken(tokens[1], "Grid dimension"),
            ParseIntToken(tokens[2], "Grid dimension")};
    return true;
}

bool IsAutoText(const std::string& rawValue) {
    const std::string value = ToLower(StripQuotes(rawValue));
    return value == "all" || value == "auto" || value == "default" ||
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
        parsed.min = ParseFiniteDoubleToken(tokens[0], "Coordinate crop lower bound");
    }
    if (!IsAutoUpperBound(tokens[1])) {
        parsed.max = ParseFiniteDoubleToken(tokens[1], "Coordinate crop upper bound");
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
        parsed.start = ParseIntToken(tokens[0], "Frame range start");
    }
    if (!IsAutoUpperBound(tokens[1])) {
        parsed.end = ParseIntToken(tokens[1], "Frame range end");
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
        const double scale = ParseFiniteDoubleToken(tokens[0], "Input scale");
        *out = {scale, scale, scale};
        return true;
    }
    if (tokens.size() == 3) {
        *out = {ParseFiniteDoubleToken(tokens[0], "Input scale"),
                ParseFiniteDoubleToken(tokens[1], "Input scale"),
                ParseFiniteDoubleToken(tokens[2], "Input scale")};
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

void ValidateParametersImpl(Parameters* params) {
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
    constexpr std::uint64_t kMaximumOutputPoints = 1000000000ull;
    std::uint64_t pointCount = 1;
    for (int dim : params->dimensions) {
        const std::uint64_t factor = static_cast<std::uint64_t>(dim);
        if (factor > kMaximumOutputPoints || pointCount > kMaximumOutputPoints / factor) {
            throw std::runtime_error(
                "resample dimensions request more than the 1,000,000,000-point safety limit.");
        }
        pointCount *= factor;
    }
    if (params->physicalDt <= 0.0 || !std::isfinite(params->physicalDt)) {
        throw std::runtime_error("physical_dt/dt must be positive when provided.");
    }
    if (params->frameRange.start.has_value() && *params->frameRange.start < 0) {
        throw std::runtime_error("frame_range/T_range start must be non-negative.");
    }
    if (params->frameRange.end.has_value() && *params->frameRange.end < 0) {
        throw std::runtime_error("frame_range/T_range end must be non-negative.");
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
    if (!std::isfinite(params->linear.tolerance) || params->linear.tolerance < 0.0) {
        throw std::runtime_error("Linear tolerance must be a finite non-negative number.");
    }
    if (params->idw.neighbors < 1) {
        throw std::runtime_error("IDW neighbors must be at least 1.");
    }
    if (params->kriging.neighbors < 1 || params->kriging.minNeighbors < 1) {
        throw std::runtime_error("Kriging neighbors and min_neighbors must be at least 1.");
    }
    if (params->kriging.neighbors > kMaximumKrigingNeighbors) {
        throw std::runtime_error(
            "Kriging neighbors cannot exceed " + std::to_string(kMaximumKrigingNeighbors) +
            " because the solver uses a dense local system.");
    }
    if (params->kriging.minNeighbors > params->kriging.neighbors) {
        throw std::runtime_error("Kriging min_neighbors cannot exceed the neighbor limit.");
    }
    if (params->idw.power <= 0.0 || !std::isfinite(params->idw.power)) {
        throw std::runtime_error("IDW power must be positive.");
    }
    if (!std::isfinite(params->idw.radius) || !std::isfinite(params->kriging.radius) ||
        !std::isfinite(params->kriging.range) || params->idw.radius < 0.0 ||
        params->kriging.radius < 0.0 || params->kriging.range < 0.0) {
        throw std::runtime_error("Interpolation radii/range must be non-negative; use 0/auto for automatic values.");
    }
    if (!std::isfinite(params->idw.nullValue)) {
        throw std::runtime_error("IDW null_value must be finite.");
    }
    if (params->kriging.sill <= 0.0 || !std::isfinite(params->kriging.sill)) {
        throw std::runtime_error("Kriging sill must be positive.");
    }
    if (!std::isfinite(params->kriging.nugget) ||
        !std::isfinite(params->kriging.regularization) ||
        params->kriging.nugget < 0.0 || params->kriging.regularization < 0.0) {
        throw std::runtime_error("Kriging nugget and regularization must be non-negative.");
    }
    if (!std::isfinite(params->kriging.nullValue)) {
        throw std::runtime_error("Kriging null_value must be finite.");
    }
    if (params->kriging.fallbackIdwPower <= 0.0 || !std::isfinite(params->kriging.fallbackIdwPower)) {
        throw std::runtime_error("Kriging fallback_idw_power must be positive.");
    }
    if (params->outputBasename.empty()) {
        params->outputBasename = DefaultBasename(params->deckFile);
    }
    const std::filesystem::path basenamePath(params->outputBasename);
    const bool containsSeparator = params->outputBasename.find('/') != std::string::npos ||
                                   params->outputBasename.find('\\') != std::string::npos;
    const bool containsTraversal = params->outputBasename.find("..") != std::string::npos;
    const bool containsPortableInvalid =
        params->outputBasename.find_first_of("<>:\"|?*") != std::string::npos;
    const bool containsControl = std::any_of(
        params->outputBasename.begin(), params->outputBasename.end(), [](unsigned char ch) {
            return ch < 0x20u || ch == 0x7fu;
        });
    if (basenamePath.is_absolute() || basenamePath.has_root_name() ||
        basenamePath.has_root_directory() || basenamePath.has_parent_path() ||
        containsSeparator || containsTraversal || containsPortableInvalid || containsControl ||
        params->outputBasename == ".") {
        throw std::runtime_error(
            "output_basename must be one portable filename component without absolute roots, "
            "path separators, '..', control characters, or reserved filename characters.");
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

void ValidateParameters(Parameters* params) {
    ValidateParametersImpl(params);
}

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
        if (pathValue.empty()) {
            throw std::runtime_error("input_path was provided but is empty.");
        }
        params.inputPath = pathValue;
    }
    if (std::string pathValue; TryString(raw, {"output_dir", "output_folder", "output.path"}, &pathValue)) {
        if (pathValue.empty()) {
            throw std::runtime_error("output_dir was provided but is empty.");
        }
        params.outputDir = pathValue;
    }
    if (std::string basename;
        TryString(raw, {"output_basename", "basename", "output.basename"}, &basename)) {
        if (basename.empty()) {
            throw std::runtime_error("output_basename was provided but is empty.");
        }
        params.outputBasename = basename;
    }
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
    TryAutoDouble(raw, {"interpolation.idw.radius", "idw.radius"}, &params.idw.radius);
    TryDouble(raw, {"interpolation.idw.power", "idw.power"}, &params.idw.power);
    TryDouble(raw, {"interpolation.idw.null_value", "idw.null_value"}, &params.idw.nullValue);

    TryInt(raw,
           {"interpolation.kriging.neighbors", "kriging.neighbors", "kriging.k"},
           &params.kriging.neighbors);
    TryInt(raw,
           {"interpolation.kriging.min_neighbors", "kriging.min_neighbors"},
           &params.kriging.minNeighbors);
    TryAutoDouble(raw, {"interpolation.kriging.radius", "kriging.radius"}, &params.kriging.radius);
    TryAutoDouble(raw, {"interpolation.kriging.range", "kriging.range"}, &params.kriging.range);
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
    ValidateParameters(&params);
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
    return out.str();
}

}  // namespace streamcenterplus::resample
