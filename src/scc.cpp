#include "cellnet/scc.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cellnet {

SccResult tarjan_scc(
    const std::vector<std::uint64_t>& offsets,
    const std::vector<NodeId>& targets) {
    if (offsets.empty()) {
        throw std::invalid_argument("SCC offsets must contain at least one entry");
    }
    const std::size_t node_count = offsets.size() - 1U;
    if (offsets.back() != targets.size()) {
        throw std::invalid_argument("SCC CSR offsets do not match target count");
    }
    if (node_count > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
        throw std::invalid_argument("SCC node count exceeds uint32 capacity");
    }

    struct Frame {
        NodeId node;
        std::uint64_t next_edge;
        NodeId parent;
    };

    constexpr std::uint32_t kUnvisited = std::numeric_limits<std::uint32_t>::max();
    std::vector<std::uint32_t> index(node_count, kUnvisited);
    std::vector<std::uint32_t> lowlink(node_count, 0U);
    std::vector<std::uint8_t> on_stack(node_count, 0U);
    std::vector<NodeId> tarjan_stack;
    std::vector<Frame> dfs_stack;
    tarjan_stack.reserve(node_count);
    dfs_stack.reserve(256U);

    SccResult result;
    result.component_of.resize(node_count, kInvalidScc);
    std::uint32_t next_index = 0U;

    for (std::size_t root_index = 0; root_index < node_count; ++root_index) {
        if (index[root_index] != kUnvisited) {
            continue;
        }
        const NodeId root = static_cast<NodeId>(root_index);
        index[root] = next_index;
        lowlink[root] = next_index;
        ++next_index;
        tarjan_stack.push_back(root);
        on_stack[root] = 1U;
        dfs_stack.push_back(Frame{root, offsets[root], kInvalidNode});

        while (!dfs_stack.empty()) {
            Frame& frame = dfs_stack.back();
            const NodeId node = frame.node;
            const std::uint64_t edge_end = offsets[static_cast<std::size_t>(node) + 1U];

            if (frame.next_edge < edge_end) {
                const NodeId target = targets[static_cast<std::size_t>(frame.next_edge)];
                ++frame.next_edge;
                if (target >= node_count) {
                    throw std::invalid_argument("SCC target node is out of range");
                }
                if (index[target] == kUnvisited) {
                    index[target] = next_index;
                    lowlink[target] = next_index;
                    ++next_index;
                    tarjan_stack.push_back(target);
                    on_stack[target] = 1U;
                    dfs_stack.push_back(Frame{target, offsets[target], node});
                    continue;
                }
                if (on_stack[target] != 0U) {
                    lowlink[node] = std::min(lowlink[node], index[target]);
                }
                continue;
            }

            const NodeId parent = frame.parent;
            dfs_stack.pop_back();
            if (parent != kInvalidNode) {
                lowlink[parent] = std::min(lowlink[parent], lowlink[node]);
            }

            if (lowlink[node] == index[node]) {
                const SccId component_id = static_cast<SccId>(result.components.size());
                result.components.emplace_back();
                auto& component = result.components.back();
                while (true) {
                    if (tarjan_stack.empty()) {
                        throw std::logic_error("Tarjan stack underflow");
                    }
                    const NodeId member = tarjan_stack.back();
                    tarjan_stack.pop_back();
                    on_stack[member] = 0U;
                    result.component_of[member] = component_id;
                    component.push_back(member);
                    if (member == node) {
                        break;
                    }
                }
                std::sort(component.begin(), component.end());
            }
        }
    }

    return result;
}

}  // namespace cellnet

