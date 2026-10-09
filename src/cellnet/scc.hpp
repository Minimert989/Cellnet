#pragma once

#include "cellnet/types.hpp"

#include <cstdint>
#include <vector>

namespace cellnet {

struct SccResult {
    std::vector<SccId> component_of;
    std::vector<std::vector<NodeId>> components;
};

[[nodiscard]] SccResult tarjan_scc(
    const std::vector<std::uint64_t>& offsets,
    const std::vector<NodeId>& targets);

}  // namespace cellnet

