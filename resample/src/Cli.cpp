#include "streamcenterplus_resample/Core.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace streamcenterplus::resample {

namespace {

std::filesystem::path ResolveDefaultDeckFile(const char* argv0) {
    std::filesystem::path executablePath = argv0 != nullptr ? std::filesystem::path(argv0) : std::filesystem::path();
    std::filesystem::path executableDir = executablePath.has_parent_path()
                                              ? executablePath.parent_path()
                                              : std::filesystem::current_path();
    const std::vector<std::filesystem::path> candidates = {
        std::filesystem::current_path() / "config" / "default.resample.scpo",
        executableDir / "config" / "default.resample.scpo",
        executableDir / "default.resample.scpo",
        std::filesystem::current_path() / "default.resample.scpo",
    };
    for (const std::filesystem::path& candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }
    return candidates.front();
}

std::string ToLower(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool HasResampleSuffix(const std::filesystem::path& path) {
    const std::string fileName = ToLower(path.filename().string());
    const std::string suffix = ".resample.scpo";
    return fileName.size() >= suffix.size() &&
           fileName.compare(fileName.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string Usage() {
    return "Usage: Resample_solver [--par <path-to-resample.scpo-file>]\n"
           "       Resample_solver --smoke-vtk-smp";
}

}  // namespace

int RunCli(int argc, char** argv) {
    try {
        std::cout.setf(std::ios::unitbuf);
        std::cerr.setf(std::ios::unitbuf);

        std::filesystem::path deckFile = ResolveDefaultDeckFile(argc >= 1 ? argv[0] : nullptr);
        bool explicitDeckFile = false;
        bool runVtkSmpSmoke = false;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i] == nullptr ? std::string() : std::string(argv[i]);
            if (arg == "--par") {
                if (runVtkSmpSmoke) {
                    throw std::runtime_error(Usage());
                }
                if (i + 1 >= argc || argv[i + 1] == nullptr || std::string(argv[i + 1]).empty()) {
                    throw std::runtime_error(Usage());
                }
                deckFile = std::filesystem::path(argv[++i]);
                explicitDeckFile = true;
            } else if (arg == "--smoke-vtk-smp" || arg == "--vtk-smp-smoke") {
                if (explicitDeckFile || runVtkSmpSmoke) {
                    throw std::runtime_error(Usage());
                }
                runVtkSmpSmoke = true;
            } else if (arg == "-h" || arg == "--help") {
                std::cout << Usage() << "\n";
                return 0;
            } else if (!explicitDeckFile && !runVtkSmpSmoke && argc == 2 && !arg.empty() && arg.front() != '-') {
                deckFile = std::filesystem::path(arg);
                explicitDeckFile = true;
            } else {
                throw std::runtime_error(Usage());
            }
        }

        if (runVtkSmpSmoke) {
            return RunVtkSmpSmokeTest();
        }

        if (!HasResampleSuffix(deckFile)) {
            throw std::runtime_error("Resample solver only accepts .resample.scpo operation files. Received: " +
                                     deckFile.string());
        }
        if (explicitDeckFile && !std::filesystem::exists(deckFile)) {
            throw std::runtime_error("Resample operation file does not exist: " + deckFile.string());
        }
        if (!explicitDeckFile) {
            std::cout << "\nWARNING: USING DEFAULT RESAMPLE OPERATION FILE: "
                      << std::filesystem::absolute(deckFile).string() << "\n";
        }

        Parameters params = LoadParameters(deckFile);
        return RunSolver(params);
    } catch (const std::exception& ex) {
        std::cerr << "Fatal error: " << ex.what() << "\n";
        return 1;
    }
}

}  // namespace streamcenterplus::resample
