#include <mpi.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gslib.h>

namespace fs = std::filesystem;

namespace {

constexpr std::size_t kNekHeaderBytes = 132;
constexpr std::size_t kNekEndianTagBytes = 4;
constexpr std::uint32_t kCacheVersion = 1;

#pragma pack(push, 1)
struct BinaryVelocityHeader {
    char magic[8];
    std::uint32_t version;
    std::uint32_t headerBytes;
    std::uint32_t nx;
    std::uint32_t ny;
    std::uint32_t nz;
    std::uint32_t reserved0;
    std::uint64_t pointCount;
    double time;
    double x0;
    double y0;
    double z0;
    double dx;
    double dy;
    double dz;
    std::uint64_t validCount;
    std::uint32_t sourceFrame;
    char reserved[20];
};
#pragma pack(pop)

static_assert(sizeof(BinaryVelocityHeader) == 128, "Unexpected FTLE cache header size");

struct Options {
    fs::path inputDir;
    fs::path outputDir;
    std::string prefix = "hemi0.f";
    int startFrame = 0;
    int endFrame = 999;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    double x0 = 0.0;
    double y0 = 0.0;
    double z0 = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double bboxTolerance = 0.01;
    double newtonTolerance = 1.0e-10;
    unsigned pointsPerBatch = 256;
    bool overwrite = false;
};

struct NekHeader {
    int wordSize = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::uint64_t elementCount = 0;
    double time = 0.0;
    std::string code;
    bool littleEndian = true;

    std::uint64_t pointsPerElement() const {
        return static_cast<std::uint64_t>(nx) * static_cast<std::uint64_t>(ny) *
               static_cast<std::uint64_t>(nz);
    }

    std::uint64_t dataBase() const {
        return kNekHeaderBytes + kNekEndianTagBytes +
               sizeof(std::int32_t) * elementCount;
    }
};

struct QueryPoint {
    double x;
    double y;
    double z;
    double r[3];
    double dist2;
    uint code;
    uint proc;
    uint el;
};

std::string Trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string JsonEscape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        switch (ch) {
        case '\\': out << "\\\\"; break;
        case '"': out << "\\\""; break;
        case '\n': out << "\\n"; break;
        case '\r': out << "\\r"; break;
        case '\t': out << "\\t"; break;
        default: out << ch; break;
        }
    }
    return out.str();
}

[[noreturn]] void Usage(const std::string& message = {}) {
    if (!message.empty()) {
        std::cerr << "Error: " << message << "\n\n";
    }
    std::cerr
        << "Usage: nek_ftle_resampler --input-dir <dir> --output-dir <dir>\n"
        << "       --start <frame> --end <frame> --dims <nx> <ny> <nz>\n"
        << "       --origin <x0> <y0> <z0> --spacing <dx> <dy> <dz>\n"
        << "       [--prefix hemi0.f] [--bbox-tol 0.01] [--newton-tol 1e-10]\n"
        << "       [--points-per-batch 256] [--overwrite]\n";
    std::exit(2);
}

int ParseInt(const char* value, const char* label) {
    try {
        std::size_t used = 0;
        const int parsed = std::stoi(value, &used);
        if (used != std::strlen(value)) {
            Usage(std::string("Invalid integer for ") + label + ": " + value);
        }
        return parsed;
    } catch (...) {
        Usage(std::string("Invalid integer for ") + label + ": " + value);
    }
}

double ParseDouble(const char* value, const char* label) {
    try {
        std::size_t used = 0;
        const double parsed = std::stod(value, &used);
        if (used != std::strlen(value) || !std::isfinite(parsed)) {
            Usage(std::string("Invalid number for ") + label + ": " + value);
        }
        return parsed;
    } catch (...) {
        Usage(std::string("Invalid number for ") + label + ": " + value);
    }
}

Options ParseOptions(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require = [&](int count) {
            if (i + count >= argc) {
                Usage("Missing value after " + arg);
            }
        };
        if (arg == "--input-dir") {
            require(1);
            options.inputDir = argv[++i];
        } else if (arg == "--output-dir") {
            require(1);
            options.outputDir = argv[++i];
        } else if (arg == "--prefix") {
            require(1);
            options.prefix = argv[++i];
        } else if (arg == "--start") {
            require(1);
            options.startFrame = ParseInt(argv[++i], "--start");
        } else if (arg == "--end") {
            require(1);
            options.endFrame = ParseInt(argv[++i], "--end");
        } else if (arg == "--dims") {
            require(3);
            options.nx = ParseInt(argv[++i], "--dims nx");
            options.ny = ParseInt(argv[++i], "--dims ny");
            options.nz = ParseInt(argv[++i], "--dims nz");
        } else if (arg == "--origin") {
            require(3);
            options.x0 = ParseDouble(argv[++i], "--origin x0");
            options.y0 = ParseDouble(argv[++i], "--origin y0");
            options.z0 = ParseDouble(argv[++i], "--origin z0");
        } else if (arg == "--spacing") {
            require(3);
            options.dx = ParseDouble(argv[++i], "--spacing dx");
            options.dy = ParseDouble(argv[++i], "--spacing dy");
            options.dz = ParseDouble(argv[++i], "--spacing dz");
        } else if (arg == "--bbox-tol") {
            require(1);
            options.bboxTolerance = ParseDouble(argv[++i], "--bbox-tol");
        } else if (arg == "--newton-tol") {
            require(1);
            options.newtonTolerance = ParseDouble(argv[++i], "--newton-tol");
        } else if (arg == "--points-per-batch") {
            require(1);
            const int parsed = ParseInt(argv[++i], "--points-per-batch");
            if (parsed <= 0) {
                Usage("--points-per-batch must be positive");
            }
            options.pointsPerBatch = static_cast<unsigned>(parsed);
        } else if (arg == "--overwrite") {
            options.overwrite = true;
        } else if (arg == "--help" || arg == "-h") {
            Usage();
        } else {
            Usage("Unknown option: " + arg);
        }
    }

    if (options.inputDir.empty() || options.outputDir.empty()) {
        Usage("--input-dir and --output-dir are required");
    }
    if (options.startFrame < 0 || options.endFrame < options.startFrame) {
        Usage("Invalid frame range");
    }
    if (options.nx < 2 || options.ny < 2 || options.nz < 2) {
        Usage("All grid dimensions must be at least 2");
    }
    if (!(options.dx > 0.0) || !(options.dy > 0.0) || !(options.dz > 0.0)) {
        Usage("All spacings must be positive");
    }
    return options;
}

fs::path SourceFramePath(const Options& options, int frame) {
    std::ostringstream name;
    name << options.prefix << std::setw(5) << std::setfill('0') << frame;
    return options.inputDir / name.str();
}

fs::path CacheFramePath(const Options& options, int frame) {
    std::ostringstream name;
    name << "frame_" << std::setw(5) << std::setfill('0') << frame << ".ftlebin";
    return options.outputDir / name.str();
}

NekHeader ReadNekHeader(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open Nek field: " + path.string());
    }
    std::array<char, kNekHeaderBytes> text{};
    std::array<char, kNekEndianTagBytes> tag{};
    input.read(text.data(), static_cast<std::streamsize>(text.size()));
    input.read(tag.data(), static_cast<std::streamsize>(tag.size()));
    if (!input) {
        throw std::runtime_error("Truncated Nek header: " + path.string());
    }

    std::istringstream tokens(std::string(text.data(), text.size()));
    std::string marker;
    std::uint64_t fileElements = 0;
    int step = 0;
    int fileId = 0;
    int fileCount = 0;
    NekHeader header;
    tokens >> marker >> header.wordSize >> header.nx >> header.ny >> header.nz >>
        header.elementCount >> fileElements >> header.time >> step >> fileId >> fileCount >> header.code;
    if (!tokens || marker != "#std") {
        throw std::runtime_error("Unsupported Nek field header: " + path.string());
    }
    float littleValue = 0.0f;
    std::memcpy(&littleValue, tag.data(), sizeof(float));
    header.littleEndian = std::abs(static_cast<double>(littleValue) - 6.54321) < 1.0e-4;
    if (header.wordSize != 4) {
        throw std::runtime_error("Only FP32 Nek fields are supported: " + path.string());
    }
    if (!header.littleEndian) {
        throw std::runtime_error("Only little-endian Nek fields are supported: " + path.string());
    }
    if (header.code.find('U') == std::string::npos) {
        throw std::runtime_error("Nek field has no velocity block: " + path.string());
    }
    return header;
}

std::uint64_t PartitionStart(std::uint64_t total, int rank, int ranks) {
    return total * static_cast<std::uint64_t>(rank) / static_cast<std::uint64_t>(ranks);
}

std::uint64_t PartitionEnd(std::uint64_t total, int rank, int ranks) {
    return total * static_cast<std::uint64_t>(rank + 1) / static_cast<std::uint64_t>(ranks);
}

void CheckMpi(int code, const std::string& operation) {
    if (code == MPI_SUCCESS) {
        return;
    }
    std::array<char, MPI_MAX_ERROR_STRING> message{};
    int length = 0;
    MPI_Error_string(code, message.data(), &length);
    throw std::runtime_error(operation + " failed: " + std::string(message.data(), static_cast<std::size_t>(length)));
}

void ReadFloatBlock(const fs::path& path,
                    std::uint64_t byteOffset,
                    std::vector<float>* values,
                    MPI_Comm communicator) {
    if (values == nullptr) {
        throw std::runtime_error("ReadFloatBlock received a null output vector");
    }
    if (values->size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("MPI float read exceeds INT_MAX");
    }
    MPI_File file = MPI_FILE_NULL;
    CheckMpi(MPI_File_open(communicator,
                           const_cast<char*>(path.string().c_str()),
                           MPI_MODE_RDONLY,
                           MPI_INFO_NULL,
                           &file),
             "MPI_File_open(" + path.string() + ")");
    const int result = MPI_File_read_at_all(file,
                                            static_cast<MPI_Offset>(byteOffset),
                                            values->data(),
                                            static_cast<int>(values->size()),
                                            MPI_FLOAT,
                                            MPI_STATUS_IGNORE);
    const int closeResult = MPI_File_close(&file);
    CheckMpi(result, "MPI_File_read_at_all(" + path.string() + ")");
    CheckMpi(closeResult, "MPI_File_close(" + path.string() + ")");
}

std::uint64_t VelocityBlockOffset(const NekHeader& header, std::uint64_t elementStart) {
    const std::uint64_t nxyz = header.pointsPerElement();
    const std::uint64_t coordinateValues = header.code.find('X') == std::string::npos
                                               ? 0
                                               : 3 * header.elementCount * nxyz;
    return header.dataBase() +
           sizeof(float) * (coordinateValues + elementStart * 3 * nxyz);
}

std::uint64_t CoordinateBlockOffset(const NekHeader& header, std::uint64_t elementStart) {
    if (header.code.find('X') == std::string::npos) {
        throw std::runtime_error("Coordinate source field has no X block");
    }
    return header.dataBase() + sizeof(float) * elementStart * 3 * header.pointsPerElement();
}

void WriteMask(const Options& options,
               const std::vector<QueryPoint>& points,
               std::uint64_t pointStart,
               MPI_Comm communicator) {
    std::vector<unsigned char> mask(points.size(), 0);
    for (std::size_t i = 0; i < points.size(); ++i) {
        mask[i] = points[i].code < 2 ? static_cast<unsigned char>(1) : static_cast<unsigned char>(0);
    }
    MPI_File file = MPI_FILE_NULL;
    const fs::path path = options.outputDir / "valid_mask.uint8.bin";
    CheckMpi(MPI_File_open(communicator,
                           const_cast<char*>(path.string().c_str()),
                           MPI_MODE_WRONLY | MPI_MODE_CREATE,
                           MPI_INFO_NULL,
                           &file),
             "MPI_File_open(mask)");
    CheckMpi(MPI_File_write_at_all(file,
                                   static_cast<MPI_Offset>(pointStart),
                                   mask.data(),
                                   static_cast<int>(mask.size()),
                                   MPI_UNSIGNED_CHAR,
                                   MPI_STATUS_IGNORE),
             "MPI_File_write_at_all(mask)");
    CheckMpi(MPI_File_close(&file), "MPI_File_close(mask)");
}

void WriteVelocityFrame(const Options& options,
                        int frame,
                        double time,
                        std::uint64_t totalPoints,
                        std::uint64_t validCount,
                        std::uint64_t pointStart,
                        const std::array<std::vector<float>, 3>& values,
                        MPI_Comm communicator) {
    const fs::path finalPath = CacheFramePath(options, frame);
    const fs::path partialPath = finalPath.string() + ".partial";
    MPI_File file = MPI_FILE_NULL;
    CheckMpi(MPI_File_open(communicator,
                           const_cast<char*>(partialPath.string().c_str()),
                           MPI_MODE_WRONLY | MPI_MODE_CREATE,
                           MPI_INFO_NULL,
                           &file),
             "MPI_File_open(cache frame)");

    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    if (rank == 0) {
        BinaryVelocityHeader header{};
        std::memcpy(header.magic, "FTLEBIN", 7);
        header.version = kCacheVersion;
        header.headerBytes = sizeof(BinaryVelocityHeader);
        header.nx = static_cast<std::uint32_t>(options.nx);
        header.ny = static_cast<std::uint32_t>(options.ny);
        header.nz = static_cast<std::uint32_t>(options.nz);
        header.pointCount = totalPoints;
        header.time = time;
        header.x0 = options.x0;
        header.y0 = options.y0;
        header.z0 = options.z0;
        header.dx = options.dx;
        header.dy = options.dy;
        header.dz = options.dz;
        header.validCount = validCount;
        header.sourceFrame = static_cast<std::uint32_t>(frame);
        CheckMpi(MPI_File_write_at(file,
                                   0,
                                   &header,
                                   static_cast<int>(sizeof(header)),
                                   MPI_BYTE,
                                   MPI_STATUS_IGNORE),
                 "MPI_File_write_at(cache header)");
    }
    CheckMpi(MPI_Barrier(communicator), "MPI_Barrier(cache header)");

    for (int component = 0; component < 3; ++component) {
        const std::uint64_t valueOffset =
            sizeof(BinaryVelocityHeader) +
            sizeof(float) * (static_cast<std::uint64_t>(component) * totalPoints + pointStart);
        CheckMpi(MPI_File_write_at_all(file,
                                       static_cast<MPI_Offset>(valueOffset),
                                       const_cast<float*>(values[static_cast<std::size_t>(component)].data()),
                                       static_cast<int>(values[static_cast<std::size_t>(component)].size()),
                                       MPI_FLOAT,
                                       MPI_STATUS_IGNORE),
                 "MPI_File_write_at_all(cache values)");
    }
    CheckMpi(MPI_File_sync(file), "MPI_File_sync(cache frame)");
    CheckMpi(MPI_File_close(&file), "MPI_File_close(cache frame)");
    CheckMpi(MPI_Barrier(communicator), "MPI_Barrier(cache close)");
    if (rank == 0) {
        std::error_code error;
        fs::remove(finalPath, error);
        error.clear();
        fs::rename(partialPath, finalPath, error);
        if (error) {
            throw std::runtime_error("Cannot finalize cache frame " + finalPath.string() + ": " + error.message());
        }
    }
    CheckMpi(MPI_Barrier(communicator), "MPI_Barrier(cache rename)");
}

void WriteMetadata(const Options& options,
                   const NekHeader& source,
                   std::uint64_t totalPoints,
                   std::uint64_t validCount,
                   int ranks) {
    std::ofstream output(options.outputDir / "grid.json");
    if (!output) {
        throw std::runtime_error("Cannot write grid.json");
    }
    output << std::setprecision(17)
           << "{\n"
           << "  \"format\": \"streamcenterplus_ftle_binary_velocity_v1\",\n"
           << "  \"version\": " << kCacheVersion << ",\n"
           << "  \"source_dir\": \"" << JsonEscape(options.inputDir.string()) << "\",\n"
           << "  \"source_prefix\": \"" << JsonEscape(options.prefix) << "\",\n"
           << "  \"source_frame_range\": [" << options.startFrame << ", " << options.endFrame << "],\n"
           << "  \"source_elements\": " << source.elementCount << ",\n"
           << "  \"source_gll_shape\": [" << source.nx << ", " << source.ny << ", " << source.nz << "],\n"
           << "  \"dimensions\": [" << options.nx << ", " << options.ny << ", " << options.nz << "],\n"
           << "  \"origin\": [" << options.x0 << ", " << options.y0 << ", " << options.z0 << "],\n"
           << "  \"spacing\": [" << options.dx << ", " << options.dy << ", " << options.dz << "],\n"
           << "  \"point_count\": " << totalPoints << ",\n"
           << "  \"valid_point_count\": " << validCount << ",\n"
           << "  \"invalid_point_count\": " << (totalPoints - validCount) << ",\n"
           << "  \"bbox_tolerance\": " << options.bboxTolerance << ",\n"
           << "  \"newton_tolerance\": " << options.newtonTolerance << ",\n"
           << "  \"mpi_ranks\": " << ranks << ",\n"
           << "  \"component_layout\": \"header128,u_float32,v_float32,w_float32\",\n"
           << "  \"mask_file\": \"valid_mask.uint8.bin\"\n"
           << "}\n";
}

}  // namespace

int main(int argc, char** argv) {
    int initialized = 0;
    try {
        CheckMpi(MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &initialized), "MPI_Init_thread");
        int rank = 0;
        int ranks = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &ranks);

        const Options options = ParseOptions(argc, argv);
        if (rank == 0) {
            fs::create_directories(options.outputDir);
        }
        CheckMpi(MPI_Barrier(MPI_COMM_WORLD), "MPI_Barrier(output directory)");

        const fs::path coordinatePath = SourceFramePath(options, 0);
        const NekHeader coordinateHeader = ReadNekHeader(coordinatePath);
        const std::uint64_t sourceStart = PartitionStart(coordinateHeader.elementCount, rank, ranks);
        const std::uint64_t sourceEnd = PartitionEnd(coordinateHeader.elementCount, rank, ranks);
        const std::uint64_t localElements = sourceEnd - sourceStart;
        const std::uint64_t nxyz = coordinateHeader.pointsPerElement();
        const std::uint64_t localSourcePoints = localElements * nxyz;

        if (localElements > std::numeric_limits<uint>::max()) {
            throw std::runtime_error("Local element count exceeds gslib uint capacity");
        }
        if (localSourcePoints * 3 > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Local coordinate read exceeds MPI int capacity");
        }

        std::vector<float> coordinateRaw(static_cast<std::size_t>(localSourcePoints * 3));
        ReadFloatBlock(coordinatePath,
                       CoordinateBlockOffset(coordinateHeader, sourceStart),
                       &coordinateRaw,
                       MPI_COMM_WORLD);
        std::array<std::vector<double>, 3> coordinates;
        for (auto& component : coordinates) {
            component.resize(static_cast<std::size_t>(localSourcePoints));
        }
        for (std::uint64_t element = 0; element < localElements; ++element) {
            for (int component = 0; component < 3; ++component) {
                const std::uint64_t sourceBase = (element * 3 + static_cast<std::uint64_t>(component)) * nxyz;
                const std::uint64_t targetBase = element * nxyz;
                std::copy_n(coordinateRaw.data() + static_cast<std::size_t>(sourceBase),
                            static_cast<std::size_t>(nxyz),
                            coordinates[static_cast<std::size_t>(component)].data() + static_cast<std::size_t>(targetBase));
            }
        }
        coordinateRaw.clear();
        coordinateRaw.shrink_to_fit();

        const std::uint64_t totalPoints = static_cast<std::uint64_t>(options.nx) *
                                          static_cast<std::uint64_t>(options.ny) *
                                          static_cast<std::uint64_t>(options.nz);
        const std::uint64_t pointStart = PartitionStart(totalPoints, rank, ranks);
        const std::uint64_t pointEnd = PartitionEnd(totalPoints, rank, ranks);
        const std::uint64_t localPoints = pointEnd - pointStart;
        if (localPoints > std::numeric_limits<uint>::max()) {
            throw std::runtime_error("Local target point count exceeds gslib uint capacity");
        }

        std::vector<QueryPoint> points(static_cast<std::size_t>(localPoints));
        for (std::uint64_t local = 0; local < localPoints; ++local) {
            const std::uint64_t global = pointStart + local;
            const std::uint64_t yz = static_cast<std::uint64_t>(options.ny) * static_cast<std::uint64_t>(options.nz);
            const std::uint64_t i = global / yz;
            const std::uint64_t remainder = global - i * yz;
            const std::uint64_t j = remainder / static_cast<std::uint64_t>(options.nz);
            const std::uint64_t k = remainder - j * static_cast<std::uint64_t>(options.nz);
            QueryPoint& point = points[static_cast<std::size_t>(local)];
            point.x = options.x0 + static_cast<double>(i) * options.dx;
            point.y = options.y0 + static_cast<double>(j) * options.dy;
            point.z = options.z0 + static_cast<double>(k) * options.dz;
        }

        struct comm gslibCommunicator;
        comm_init(&gslibCommunicator, MPI_COMM_WORLD);
        const double* coordinateBases[3] = {
            coordinates[0].data(), coordinates[1].data(), coordinates[2].data()
        };
        const unsigned shape[3] = {
            static_cast<unsigned>(coordinateHeader.nx),
            static_cast<unsigned>(coordinateHeader.ny),
            static_cast<unsigned>(coordinateHeader.nz)
        };
        const unsigned expandedShape[3] = {
            static_cast<unsigned>(2 * coordinateHeader.nx),
            static_cast<unsigned>(2 * coordinateHeader.ny),
            static_cast<unsigned>(2 * coordinateHeader.nz)
        };
        const std::uint64_t hashCandidate = std::max<std::uint64_t>(localSourcePoints, 1024);
        const uint hashSize = static_cast<uint>(std::min<std::uint64_t>(hashCandidate, std::numeric_limits<uint>::max()));
        struct findpts_data_3* finder = findpts_setup_3(
            &gslibCommunicator,
            coordinateBases,
            shape,
            static_cast<uint>(localElements),
            expandedShape,
            options.bboxTolerance,
            hashSize,
            hashSize,
            options.pointsPerBatch,
            options.newtonTolerance);
        if (finder == nullptr) {
            throw std::runtime_error("findpts_setup_3 returned null");
        }

        const double* queryBases[3] = {
            points.empty() ? nullptr : &points[0].x,
            points.empty() ? nullptr : &points[0].y,
            points.empty() ? nullptr : &points[0].z
        };
        const unsigned queryStrides[3] = {
            sizeof(QueryPoint), sizeof(QueryPoint), sizeof(QueryPoint)
        };
        findpts_3(points.empty() ? nullptr : &points[0].code, sizeof(QueryPoint),
                  points.empty() ? nullptr : &points[0].proc, sizeof(QueryPoint),
                  points.empty() ? nullptr : &points[0].el, sizeof(QueryPoint),
                  points.empty() ? nullptr : points[0].r, sizeof(QueryPoint),
                  points.empty() ? nullptr : &points[0].dist2, sizeof(QueryPoint),
                  queryBases, queryStrides, static_cast<uint>(localPoints), finder);

        std::uint64_t localValid = 0;
        double localMaxDistance = 0.0;
        for (const QueryPoint& point : points) {
            if (point.code < 2) {
                ++localValid;
                localMaxDistance = std::max(localMaxDistance, std::sqrt(std::max(0.0, point.dist2)));
            }
        }
        std::uint64_t validCount = 0;
        double maxDistance = 0.0;
        CheckMpi(MPI_Allreduce(&localValid, &validCount, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD),
                 "MPI_Allreduce(valid points)");
        CheckMpi(MPI_Allreduce(&localMaxDistance, &maxDistance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD),
                 "MPI_Allreduce(max distance)");
        WriteMask(options, points, pointStart, MPI_COMM_WORLD);

        if (rank == 0) {
            WriteMetadata(options, coordinateHeader, totalPoints, validCount, ranks);
            std::cout << "[mapping] source_elements=" << coordinateHeader.elementCount
                      << " target_points=" << totalPoints
                      << " valid=" << validCount
                      << " invalid=" << (totalPoints - validCount)
                      << " max_distance=" << std::setprecision(8) << maxDistance
                      << " mpi_ranks=" << ranks << "\n";
        }

        std::array<std::vector<double>, 3> localVelocity;
        for (auto& component : localVelocity) {
            component.resize(static_cast<std::size_t>(localSourcePoints));
        }
        std::array<std::vector<double>, 3> evaluated;
        std::array<std::vector<float>, 3> cacheValues;
        for (int component = 0; component < 3; ++component) {
            evaluated[static_cast<std::size_t>(component)].resize(static_cast<std::size_t>(localPoints));
            cacheValues[static_cast<std::size_t>(component)].resize(static_cast<std::size_t>(localPoints));
        }
        std::vector<float> velocityRaw(static_cast<std::size_t>(localSourcePoints * 3));
        const std::uint64_t expectedFrameBytes = sizeof(BinaryVelocityHeader) + sizeof(float) * 3 * totalPoints;

        for (int frame = options.startFrame; frame <= options.endFrame; ++frame) {
            const fs::path outputPath = CacheFramePath(options, frame);
            int skip = 0;
            if (rank == 0 && !options.overwrite) {
                std::error_code error;
                skip = fs::exists(outputPath, error) && !error && fs::file_size(outputPath, error) == expectedFrameBytes ? 1 : 0;
            }
            CheckMpi(MPI_Bcast(&skip, 1, MPI_INT, 0, MPI_COMM_WORLD), "MPI_Bcast(skip frame)");
            if (skip != 0) {
                if (rank == 0) {
                    std::cout << "[skip] frame=" << frame << " path=" << outputPath.string() << "\n";
                }
                continue;
            }

            const fs::path sourcePath = SourceFramePath(options, frame);
            const NekHeader header = ReadNekHeader(sourcePath);
            if (header.elementCount != coordinateHeader.elementCount ||
                header.nx != coordinateHeader.nx || header.ny != coordinateHeader.ny || header.nz != coordinateHeader.nz) {
                throw std::runtime_error("Nek mesh dimensions changed at " + sourcePath.string());
            }
            ReadFloatBlock(sourcePath,
                           VelocityBlockOffset(header, sourceStart),
                           &velocityRaw,
                           MPI_COMM_WORLD);
            for (std::uint64_t element = 0; element < localElements; ++element) {
                for (int component = 0; component < 3; ++component) {
                    const std::uint64_t sourceBase = (element * 3 + static_cast<std::uint64_t>(component)) * nxyz;
                    const std::uint64_t targetBase = element * nxyz;
                    std::copy_n(velocityRaw.data() + static_cast<std::size_t>(sourceBase),
                                static_cast<std::size_t>(nxyz),
                                localVelocity[static_cast<std::size_t>(component)].data() + static_cast<std::size_t>(targetBase));
                }
            }

            for (int component = 0; component < 3; ++component) {
                findpts_eval_3(evaluated[static_cast<std::size_t>(component)].data(), sizeof(double),
                               points.empty() ? nullptr : &points[0].code, sizeof(QueryPoint),
                               points.empty() ? nullptr : &points[0].proc, sizeof(QueryPoint),
                               points.empty() ? nullptr : &points[0].el, sizeof(QueryPoint),
                               points.empty() ? nullptr : points[0].r, sizeof(QueryPoint),
                               static_cast<uint>(localPoints),
                               localVelocity[static_cast<std::size_t>(component)].data(), finder);
                for (std::size_t i = 0; i < points.size(); ++i) {
                    cacheValues[static_cast<std::size_t>(component)][i] =
                        points[i].code < 2 ? static_cast<float>(evaluated[static_cast<std::size_t>(component)][i]) : 0.0f;
                }
            }

            WriteVelocityFrame(options,
                               frame,
                               header.time,
                               totalPoints,
                               validCount,
                               pointStart,
                               cacheValues,
                               MPI_COMM_WORLD);
            if (rank == 0) {
                std::cout << "[frame] " << (frame - options.startFrame + 1) << "/"
                          << (options.endFrame - options.startFrame + 1)
                          << " source=" << sourcePath.filename().string()
                          << " t=" << std::setprecision(12) << header.time
                          << " output=" << outputPath.filename().string() << "\n";
            }
        }

        findpts_free_3(finder);
        comm_free(&gslibCommunicator);
        MPI_Finalize();
        return 0;
    } catch (const std::exception& error) {
        int rank = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        std::cerr << "[rank " << rank << "] fatal: " << error.what() << "\n";
        int finalized = 0;
        MPI_Finalized(&finalized);
        if (!finalized) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        return 1;
    }
}
