#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

namespace streamcenterplus::resample::integrity {

class Sha256 {
public:
    Sha256();

    void Update(const void* data, std::size_t size);
    void Update(const std::string& value);
    std::string FinalHex();

private:
    void Transform(const unsigned char block[64]);

    std::array<std::uint32_t, 8> state_{};
    std::array<unsigned char, 64> buffer_{};
    std::uint64_t totalBytes_ = 0;
    std::size_t bufferedBytes_ = 0;
    bool finalized_ = false;
};

std::string Sha256String(const std::string& value);
std::string Sha256File(const std::filesystem::path& path);

std::string HexEncode(const std::string& value);
std::string HexDecode(const std::string& value);

std::filesystem::path TemporarySibling(const std::filesystem::path& destination);
void AtomicReplaceFile(const std::filesystem::path& temporary,
                       const std::filesystem::path& destination);

}  // namespace streamcenterplus::resample::integrity
