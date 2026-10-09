#include "cellnet/synthetic.hpp"

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <vector>

namespace cellnet {

NetworkSpec generate_synthetic_network(const SyntheticConfig& config) {
    if (config.node_count < 2U) {
        throw std::invalid_argument("synthetic network requires at least two nodes");
    }
    if (config.module_size < 2U || config.average_degree == 0U) {
        throw std::invalid_argument("module size must be >= 2 and degree must be positive");
    }

    NetworkSpec spec;
    spec.nodes.reserve(config.node_count);
    spec.reactions.reserve(config.node_count);
    std::mt19937_64 random(config.seed);
    std::uniform_real_distribution<float> anchor_value(0.25F, 0.85F);

    for (std::uint32_t node_index = 0; node_index < config.node_count; ++node_index) {
        const std::uint32_t local = node_index % config.module_size;
        NodeSpec node;
        node.external_id = node_index;
        node.name = "N" + std::to_string(node_index);
        node.entity = "E" + std::to_string(node_index);
        node.state = local == 0U ? "signal" : "active";
        node.type = local == 0U ? NodeType::ExternalSignal : NodeType::ProteinState;
        node.compartment = local == 0U ? Compartment::Extracellular : Compartment::Cytosol;
        node.initial_value = local == 0U ? anchor_value(random) : 0.0F;
        spec.nodes.push_back(std::move(node));
    }
    spec.nodes.back().type = NodeType::Phenotype;
    spec.nodes.back().name = "SyntheticPhenotype";

    // One stable relaxation reaction per non-anchor node. Local cycles, hubs, feed-forward
    // connections, and sparse forward cross-module links are encoded as additional inputs.
    for (std::uint32_t target = 0; target < config.node_count; ++target) {
        const std::uint32_t module_begin = (target / config.module_size) * config.module_size;
        const std::uint32_t module_end =
            std::min(config.node_count, module_begin + config.module_size);
        const std::uint32_t local = target - module_begin;
        if (local == 0U) {
            continue;
        }

        ReactionSpec reaction;
        reaction.external_id = spec.reactions.size();
        reaction.primitive = Primitive::Modulate;
        reaction.input_op = RegulationOp::ContinuousOr;
        reaction.k = 1.0F;
        reaction.outputs.push_back(WeightedNodeSpec{target, 1.0F});
        reaction.inputs.push_back(WeightedNodeSpec{target - 1U, 1.0F});

        // Close a local feedback loop without joining separate modules.
        if (local == 1U && module_end - module_begin > 3U) {
            reaction.inputs.push_back(WeightedNodeSpec{module_end - 1U, 0.2F});
        }

        // The first dynamic node acts as a hub inside each module.
        if (local > 2U && reaction.inputs.size() < config.average_degree) {
            reaction.inputs.push_back(WeightedNodeSpec{module_begin + 1U, 0.35F});
        }

        while (reaction.inputs.size() < config.average_degree && local > 1U) {
            std::uniform_int_distribution<std::uint32_t> choose(module_begin, target - 1U);
            const std::uint32_t candidate = choose(random);
            const auto duplicate = std::find_if(
                reaction.inputs.begin(), reaction.inputs.end(),
                [candidate](const WeightedNodeSpec& input) {
                    return input.node_external_id == candidate;
                });
            if (duplicate == reaction.inputs.end()) {
                reaction.inputs.push_back(WeightedNodeSpec{candidate, 0.25F});
            } else if (reaction.inputs.size() >= local) {
                break;
            }
        }

        // Sparse feed-forward cross-module edge.
        if (local == 1U && module_begin >= config.module_size) {
            const std::uint32_t previous_hub = module_begin - config.module_size + 1U;
            reaction.inputs.push_back(WeightedNodeSpec{previous_hub, 0.15F});
        }
        spec.reactions.push_back(std::move(reaction));
    }

    return spec;
}

}  // namespace cellnet

