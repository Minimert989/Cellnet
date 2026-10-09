#pragma once

#include "cellnet/network.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>

namespace cellnet {

class BinaryFormatError : public std::runtime_error {
public:
    explicit BinaryFormatError(const std::string& message) : std::runtime_error(message) {}
};

void save_compiled_network(const CompiledNetwork& network,
                           const std::filesystem::path& path);

[[nodiscard]] CompiledNetwork load_compiled_network(const std::filesystem::path& path);

}  // namespace cellnet

