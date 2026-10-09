#pragma once

#include "cellnet/network.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>

namespace cellnet {

class ParseError : public std::runtime_error {
public:
    explicit ParseError(const std::string& message) : std::runtime_error(message) {}
};

[[nodiscard]] NetworkSpec parse_csv_directory(const std::filesystem::path& directory);

[[nodiscard]] NodeType parse_node_type(const std::string& text);
[[nodiscard]] Compartment parse_compartment(const std::string& text);
[[nodiscard]] Primitive parse_primitive(const std::string& text);
[[nodiscard]] RegulationOp parse_regulation_op(const std::string& text);

}  // namespace cellnet

