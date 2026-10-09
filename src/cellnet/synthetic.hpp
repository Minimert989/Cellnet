#pragma once

#include "cellnet/network.hpp"

#include <cstdint>

namespace cellnet {

struct SyntheticConfig {
    std::uint32_t node_count{400};
    std::uint32_t module_size{40};
    std::uint32_t average_degree{3};
    std::uint64_t seed{42};
};

[[nodiscard]] NetworkSpec generate_synthetic_network(const SyntheticConfig& config);

}  // namespace cellnet

