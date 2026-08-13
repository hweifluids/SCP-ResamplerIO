#include "Integrity.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <system_error>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace streamcenterplus::resample::integrity {

namespace {

constexpr std::array<std::uint32_t, 64> kRoundConstants = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
    0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
    0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
    0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
    0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
    0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u,
};

std::uint32_t RotateRight(std::uint32_t value, unsigned int bits) {
    return (value >> bits) | (value << (32u - bits));
}

unsigned long ProcessId() {
#ifdef _WIN32
    return static_cast<unsigned long>(GetCurrentProcessId());
#else
    return static_cast<unsigned long>(getpid());
#endif
}

int HexDigit(char ch) {
    if (ch >= '0' && ch <= '9') {
        return ch - '0';
    }
    if (ch >= 'a' && ch <= 'f') {
        return ch - 'a' + 10;
    }
    if (ch >= 'A' && ch <= 'F') {
        return ch - 'A' + 10;
    }
    return -1;
}

}  // namespace

Sha256::Sha256()
    : state_{0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
             0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u} {
}

void Sha256::Transform(const unsigned char block[64]) {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16; ++i) {
        const std::size_t offset = i * 4;
        words[i] = (static_cast<std::uint32_t>(block[offset]) << 24u) |
                   (static_cast<std::uint32_t>(block[offset + 1]) << 16u) |
                   (static_cast<std::uint32_t>(block[offset + 2]) << 8u) |
                   static_cast<std::uint32_t>(block[offset + 3]);
    }
    for (std::size_t i = 16; i < words.size(); ++i) {
        const std::uint32_t s0 = RotateRight(words[i - 15], 7) ^
                                 RotateRight(words[i - 15], 18) ^
                                 (words[i - 15] >> 3u);
        const std::uint32_t s1 = RotateRight(words[i - 2], 17) ^
                                 RotateRight(words[i - 2], 19) ^
                                 (words[i - 2] >> 10u);
        words[i] = words[i - 16] + s0 + words[i - 7] + s1;
    }

    std::uint32_t a = state_[0];
    std::uint32_t b = state_[1];
    std::uint32_t c = state_[2];
    std::uint32_t d = state_[3];
    std::uint32_t e = state_[4];
    std::uint32_t f = state_[5];
    std::uint32_t g = state_[6];
    std::uint32_t h = state_[7];
    for (std::size_t i = 0; i < words.size(); ++i) {
        const std::uint32_t sum1 = RotateRight(e, 6) ^ RotateRight(e, 11) ^ RotateRight(e, 25);
        const std::uint32_t choose = (e & f) ^ ((~e) & g);
        const std::uint32_t temp1 = h + sum1 + choose + kRoundConstants[i] + words[i];
        const std::uint32_t sum0 = RotateRight(a, 2) ^ RotateRight(a, 13) ^ RotateRight(a, 22);
        const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
        const std::uint32_t temp2 = sum0 + majority;
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }
    state_[0] += a;
    state_[1] += b;
    state_[2] += c;
    state_[3] += d;
    state_[4] += e;
    state_[5] += f;
    state_[6] += g;
    state_[7] += h;
}

void Sha256::Update(const void* data, std::size_t size) {
    if (finalized_) {
        throw std::runtime_error("Cannot update a finalized SHA-256 digest.");
    }
    if (size == 0) {
        return;
    }
    if (data == nullptr ||
        static_cast<std::uintmax_t>(size) >
            static_cast<std::uintmax_t>(std::numeric_limits<std::uint64_t>::max()) ||
        static_cast<std::uint64_t>(size) >
            std::numeric_limits<std::uint64_t>::max() - totalBytes_) {
        throw std::runtime_error("SHA-256 input is null or exceeds the supported size.");
    }
    const auto* bytes = static_cast<const unsigned char*>(data);
    totalBytes_ += static_cast<std::uint64_t>(size);

    if (bufferedBytes_ > 0) {
        const std::size_t copied = std::min(size, buffer_.size() - bufferedBytes_);
        std::memcpy(buffer_.data() + bufferedBytes_, bytes, copied);
        bufferedBytes_ += copied;
        bytes += copied;
        size -= copied;
        if (bufferedBytes_ == buffer_.size()) {
            Transform(buffer_.data());
            bufferedBytes_ = 0;
        }
    }
    while (size >= buffer_.size()) {
        Transform(bytes);
        bytes += buffer_.size();
        size -= buffer_.size();
    }
    if (size > 0) {
        std::memcpy(buffer_.data(), bytes, size);
        bufferedBytes_ = size;
    }
}

void Sha256::Update(const std::string& value) {
    Update(value.data(), value.size());
}

std::string Sha256::FinalHex() {
    if (finalized_) {
        throw std::runtime_error("SHA-256 digest was already finalized.");
    }
    if (totalBytes_ > std::numeric_limits<std::uint64_t>::max() / 8ull) {
        throw std::runtime_error("SHA-256 bit length overflows uint64.");
    }
    const std::uint64_t bitLength = totalBytes_ * 8ull;
    buffer_[bufferedBytes_++] = 0x80u;
    if (bufferedBytes_ > 56) {
        std::fill(buffer_.begin() + static_cast<std::ptrdiff_t>(bufferedBytes_), buffer_.end(), 0u);
        Transform(buffer_.data());
        bufferedBytes_ = 0;
    }
    std::fill(buffer_.begin() + static_cast<std::ptrdiff_t>(bufferedBytes_), buffer_.begin() + 56, 0u);
    for (int i = 0; i < 8; ++i) {
        buffer_[56 + static_cast<std::size_t>(i)] =
            static_cast<unsigned char>((bitLength >> (56 - 8 * i)) & 0xffu);
    }
    Transform(buffer_.data());
    finalized_ = true;

    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (std::uint32_t value : state_) {
        output << std::setw(8) << value;
    }
    return output.str();
}

std::string Sha256String(const std::string& value) {
    Sha256 digest;
    digest.Update(value);
    return digest.FinalHex();
}

std::string Sha256File(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("Cannot open file for checksum: " + path.string());
    }
    Sha256 digest;
    // Keep the buffer well below the default 1 MiB Windows thread stack;
    // callers may checksum a frame from an already deep VTK call path.
    std::array<char, 64u << 10> buffer{};
    while (input) {
        input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize count = input.gcount();
        if (count > 0) {
            digest.Update(buffer.data(), static_cast<std::size_t>(count));
        }
    }
    if (!input.eof()) {
        throw std::runtime_error("Failed while checksumming file: " + path.string());
    }
    return digest.FinalHex();
}

std::string HexEncode(const std::string& value) {
    static constexpr char kDigits[] = "0123456789abcdef";
    std::string encoded;
    if (value.size() > std::numeric_limits<std::size_t>::max() / 2) {
        throw std::runtime_error("Hexadecimal encoding size overflows size_t.");
    }
    encoded.reserve(value.size() * 2);
    for (unsigned char ch : value) {
        encoded.push_back(kDigits[ch >> 4u]);
        encoded.push_back(kDigits[ch & 0x0fu]);
    }
    return encoded;
}

std::string HexDecode(const std::string& value) {
    if (value.size() % 2 != 0) {
        throw std::runtime_error("Manifest hexadecimal value has odd length.");
    }
    std::string decoded;
    decoded.reserve(value.size() / 2);
    for (std::size_t i = 0; i < value.size(); i += 2) {
        const int high = HexDigit(value[i]);
        const int low = HexDigit(value[i + 1]);
        if (high < 0 || low < 0) {
            throw std::runtime_error("Manifest hexadecimal value contains a non-hexadecimal digit.");
        }
        decoded.push_back(static_cast<char>((high << 4) | low));
    }
    return decoded;
}

std::filesystem::path TemporarySibling(const std::filesystem::path& destination) {
    static std::atomic<unsigned long long> counter{0};
    const auto ticks = static_cast<unsigned long long>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    for (int attempt = 0; attempt < 64; ++attempt) {
        std::ostringstream name;
        name << destination.filename().string() << ".tmp." << ProcessId() << "."
             << ticks << "." << counter.fetch_add(1, std::memory_order_relaxed);
        const std::filesystem::path candidate = destination.parent_path() / name.str();
        std::error_code error;
        if (!std::filesystem::exists(candidate, error) && !error) {
            return candidate;
        }
    }
    throw std::runtime_error("Cannot allocate a unique temporary file beside: " +
                             destination.string());
}

void AtomicReplaceFile(const std::filesystem::path& temporary,
                       const std::filesystem::path& destination) {
    if (!std::filesystem::is_regular_file(temporary)) {
        throw std::runtime_error("Atomic publish source is not a regular file: " + temporary.string());
    }
#ifdef _WIN32
    if (!MoveFileExW(temporary.c_str(), destination.c_str(),
                     MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH)) {
        const DWORD code = GetLastError();
        throw std::runtime_error("Atomic file replacement failed for '" + destination.string() +
                                 "' (Windows error " + std::to_string(code) + ").");
    }
#else
    std::error_code error;
    std::filesystem::rename(temporary, destination, error);
    if (error) {
        throw std::runtime_error("Atomic file replacement failed for '" + destination.string() +
                                 "': " + error.message());
    }
#endif
}

}  // namespace streamcenterplus::resample::integrity
