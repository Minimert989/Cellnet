#pragma once

#include "cellnet/types.hpp"

#include <cstdint>
#include <vector>

namespace cellnet {

struct NodeContextMultiplier {
    NodeId node{kInvalidNode};
    float multiplier{1.0F};
};

struct ReactionContextMultiplier {
    ReactionId reaction{kInvalidReaction};
    float multiplier{1.0F};
};

struct Context {
    std::uint64_t context_id{0};
    std::vector<NodeContextMultiplier> node_multipliers;
    std::vector<ReactionContextMultiplier> reaction_multipliers;
};

}  // namespace cellnet

