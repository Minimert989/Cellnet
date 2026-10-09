#include "cellnet/compiler.hpp"

#include "cellnet/scc.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cellnet {
namespace {

template <typename T>
void require_finite(T value, const std::string& label) {
    if (!std::isfinite(value)) {
        throw CompileError(label + " must be finite");
    }
}

std::uint32_t checked_u32(std::size_t value, const char* label) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw CompileError(std::string(label) + " exceeds uint32 capacity");
    }
    return static_cast<std::uint32_t>(value);
}

std::vector<NodeId> affected_nodes(const ReactionSpec& reaction,
                                   const std::unordered_map<std::uint64_t, NodeId>& node_ids) {
    std::vector<NodeId> affected;
    const bool consumes_inputs = reaction.primitive == Primitive::Transfer ||
                                 reaction.primitive == Primitive::Destroy ||
                                 reaction.primitive == Primitive::Bind ||
                                 reaction.primitive == Primitive::Unbind ||
                                 reaction.primitive == Primitive::Move;
    if (consumes_inputs) {
        affected.reserve(reaction.inputs.size() + reaction.outputs.size());
        for (const auto& input : reaction.inputs) {
            affected.push_back(node_ids.at(input.node_external_id));
        }
    } else {
        affected.reserve(reaction.outputs.size());
    }
    if (reaction.primitive != Primitive::Destroy) {
        for (const auto& output : reaction.outputs) {
            affected.push_back(node_ids.at(output.node_external_id));
        }
    }
    std::sort(affected.begin(), affected.end());
    affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
    return affected;
}

template <typename T>
std::size_t vector_bytes(const std::vector<T>& values) {
    return values.capacity() * sizeof(T);
}

void build_execution_plans(CompiledNetwork& network) {
    network.execution_plans.assign(network.sccs.size(), SccExecutionPlan{});
    network.plan_modulate_reactions.clear();
    network.plan_modulate_slots.clear();
    network.plan_modulate_output_locals.clear();
    network.plan_modulate_output_stoich.clear();
    network.plan_transfer_reactions.clear();
    network.plan_transfer_slots.clear();
    network.plan_produce_reactions.clear();
    network.plan_produce_slots.clear();
    network.plan_destroy_reactions.clear();
    network.plan_destroy_slots.clear();
    network.plan_generic_reactions.clear();
    network.plan_generic_slots.clear();

    for (std::size_t scc_index = 0; scc_index < network.sccs.size(); ++scc_index) {
        const auto& block = network.sccs[scc_index];
        auto& plan = network.execution_plans[scc_index];
        plan.modulate_offset = static_cast<std::uint32_t>(network.plan_modulate_reactions.size());
        plan.transfer_offset = static_cast<std::uint32_t>(network.plan_transfer_reactions.size());
        plan.produce_offset = static_cast<std::uint32_t>(network.plan_produce_reactions.size());
        plan.destroy_offset = static_cast<std::uint32_t>(network.plan_destroy_reactions.size());
        plan.generic_offset = static_cast<std::uint32_t>(network.plan_generic_reactions.size());

        for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
            const ReactionId reaction_id = network.scc_reactions[block.reaction_offset + slot];
            const Reaction& reaction = network.reactions[reaction_id];
            const bool simple_modulate =
                (reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) &&
                reaction.output_count == 1U && reaction.catalyst_count == 0U &&
                reaction.pos_reg_count == 0U && reaction.neg_reg_count == 0U &&
                network.node_to_scc[network.reaction_outputs[reaction.output_offset]] == scc_index;
            if (simple_modulate) {
                network.plan_modulate_reactions.push_back(reaction_id);
                network.plan_modulate_slots.push_back(block.reaction_offset + slot);
                const NodeId output = network.reaction_outputs[reaction.output_offset];
                network.plan_modulate_output_locals.push_back(network.node_position_in_scc[output]);
                network.plan_modulate_output_stoich.push_back(
                    network.output_stoichiometry[reaction.output_offset]);
            } else if (reaction.primitive == Primitive::Transfer && reaction.catalyst_count == 0U &&
                       reaction.pos_reg_count == 0U && reaction.neg_reg_count == 0U) {
                network.plan_transfer_reactions.push_back(reaction_id);
                network.plan_transfer_slots.push_back(block.reaction_offset + slot);
            } else if (reaction.primitive == Primitive::Produce && reaction.catalyst_count == 0U &&
                       reaction.pos_reg_count == 0U && reaction.neg_reg_count == 0U) {
                network.plan_produce_reactions.push_back(reaction_id);
                network.plan_produce_slots.push_back(block.reaction_offset + slot);
            } else if (reaction.primitive == Primitive::Destroy && reaction.catalyst_count == 0U &&
                       reaction.pos_reg_count == 0U && reaction.neg_reg_count == 0U) {
                network.plan_destroy_reactions.push_back(reaction_id);
                network.plan_destroy_slots.push_back(block.reaction_offset + slot);
            } else {
                network.plan_generic_reactions.push_back(reaction_id);
                network.plan_generic_slots.push_back(block.reaction_offset + slot);
            }
        }

        plan.modulate_count = static_cast<std::uint32_t>(network.plan_modulate_reactions.size()) -
                              plan.modulate_offset;
        plan.transfer_count = static_cast<std::uint32_t>(network.plan_transfer_reactions.size()) -
                              plan.transfer_offset;
        plan.produce_count = static_cast<std::uint32_t>(network.plan_produce_reactions.size()) -
                             plan.produce_offset;
        plan.destroy_count = static_cast<std::uint32_t>(network.plan_destroy_reactions.size()) -
                             plan.destroy_offset;
        plan.generic_count = static_cast<std::uint32_t>(network.plan_generic_reactions.size()) -
                             plan.generic_offset;
    }

    std::vector<std::pair<std::uint32_t, std::uint32_t>> node_slot_pairs;
    for (std::size_t scc_index = 0; scc_index < network.sccs.size(); ++scc_index) {
        const auto& block = network.sccs[scc_index];
        for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
            const Reaction& reaction =
                network.reactions[network.scc_reactions[block.reaction_offset + slot]];
            std::vector<NodeId> dependencies;
            dependencies.reserve(reaction.input_count + reaction.catalyst_count +
                                 reaction.pos_reg_count + reaction.neg_reg_count +
                                 reaction.output_count);
            for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                dependencies.push_back(network.reaction_inputs[reaction.input_offset + index]);
            }
            for (std::uint32_t index = 0; index < reaction.catalyst_count; ++index) {
                dependencies.push_back(network.catalysts[reaction.catalyst_offset + index]);
            }
            for (std::uint32_t index = 0; index < reaction.pos_reg_count; ++index) {
                dependencies.push_back(network.positive_regulators[reaction.pos_reg_offset + index]);
            }
            for (std::uint32_t index = 0; index < reaction.neg_reg_count; ++index) {
                dependencies.push_back(network.negative_regulators[reaction.neg_reg_offset + index]);
            }
            // A relaxation reaction also depends on its current output state.
            if (reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) {
                for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                    dependencies.push_back(network.reaction_outputs[reaction.output_offset + index]);
                }
            }
            std::sort(dependencies.begin(), dependencies.end());
            dependencies.erase(std::unique(dependencies.begin(), dependencies.end()), dependencies.end());
            for (const NodeId dependency : dependencies) {
                if (network.node_to_scc[dependency] == scc_index) {
                    const std::uint32_t flat_node = block.node_offset +
                                                    network.node_position_in_scc[dependency];
                    node_slot_pairs.emplace_back(flat_node, block.reaction_offset + slot);
                }
            }
        }
    }
    std::sort(node_slot_pairs.begin(), node_slot_pairs.end());
    node_slot_pairs.erase(std::unique(node_slot_pairs.begin(), node_slot_pairs.end()),
                          node_slot_pairs.end());
    network.scc_node_reaction_offsets.assign(network.scc_nodes.size() + 1U, 0U);
    for (const auto& pair : node_slot_pairs) {
        ++network.scc_node_reaction_offsets[static_cast<std::size_t>(pair.first) + 1U];
    }
    for (std::size_t index = 1; index < network.scc_node_reaction_offsets.size(); ++index) {
        network.scc_node_reaction_offsets[index] += network.scc_node_reaction_offsets[index - 1U];
    }
    network.scc_node_reaction_slots.reserve(node_slot_pairs.size());
    for (const auto& pair : node_slot_pairs) {
        network.scc_node_reaction_slots.push_back(pair.second);
    }
}

struct LinearSignalForm {
    bool valid{false};
    float constant{0.0F};
    std::vector<std::pair<NodeId, float>> terms;
};

bool is_simple_primitive(Primitive primitive) {
    return primitive == Primitive::Modulate || primitive == Primitive::Delay ||
           primitive == Primitive::Produce || primitive == Primitive::Destroy ||
           primitive == Primitive::Transfer;
}

LinearSignalForm linear_signal_form(const CompiledNetwork& network, const Reaction& reaction) {
    LinearSignalForm form;
    if (reaction.catalyst_count != 0U || reaction.pos_reg_count != 0U ||
        reaction.neg_reg_count != 0U) {
        return form;
    }
    const float scale = reaction.context_weight * reaction.confidence;
    form.valid = true;
    if (reaction.input_count == 0U) {
        form.constant = scale;
        return form;
    }

    auto add_term = [&form, &network, &reaction, scale](std::uint32_t index, float weight) {
        form.terms.emplace_back(network.reaction_inputs[reaction.input_offset + index],
                                scale * weight);
    };
    switch (reaction.input_op) {
        case RegulationOp::None:
            if (reaction.input_count != 1U) {
                form.valid = false;
                return form;
            }
            if (network.input_stoichiometry[reaction.input_offset] < 0.0F ||
                network.input_stoichiometry[reaction.input_offset] > 1.0F) {
                form.valid = false;
                return form;
            }
            add_term(0U, network.input_stoichiometry[reaction.input_offset]);
            return form;
        case RegulationOp::Product:
        case RegulationOp::ContinuousAnd:
            if (reaction.input_count != 1U ||
                network.input_stoichiometry[reaction.input_offset] != 1.0F) {
                form.valid = false;
                return form;
            }
            add_term(0U, 1.0F);
            return form;
        case RegulationOp::ContinuousOr:
            if (reaction.input_count != 1U ||
                network.input_stoichiometry[reaction.input_offset] < 0.0F ||
                network.input_stoichiometry[reaction.input_offset] > 1.0F) {
                form.valid = false;
                return form;
            }
            add_term(0U, network.input_stoichiometry[reaction.input_offset]);
            return form;
        case RegulationOp::Sum: {
            float total_weight = 0.0F;
            for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                const float weight = network.input_stoichiometry[reaction.input_offset + index];
                if (weight < 0.0F) {
                    form.valid = false;
                    return form;
                }
                total_weight += weight;
                add_term(index, weight);
            }
            if (total_weight > 1.0F) form.valid = false;
            return form;
        }
        default:
            form.valid = false;
            return form;
    }
}

void add_matrix_term(std::vector<float>& matrix,
                     std::vector<float>& affine,
                     const CompiledNetwork& network,
                     SccId scc,
                     std::uint32_t dimension,
                     std::uint32_t row,
                     NodeId column,
                     float coefficient,
                     bool& has_external,
                     bool& has_affine) {
    if (network.node_to_scc[column] == scc) {
        const std::uint32_t local = network.node_position_in_scc[column];
        if (!matrix.empty()) matrix[static_cast<std::size_t>(row) * dimension + local] += coefficient;
    } else {
        has_external = true;
        has_affine = true;
        affine[row] += coefficient * network.nodes[column].initial_value;
    }
}

bool add_linear_reaction(const CompiledNetwork& network,
                         SccId scc,
                         const Reaction& reaction,
                         float dt,
                         std::vector<float>& matrix,
                         std::vector<float>& affine,
                         bool& has_external,
                         bool& has_affine) {
    const SccBlock& block = network.sccs[scc];
    const std::uint32_t dimension = block.node_count;
    const LinearSignalForm signal = linear_signal_form(network, reaction);
    if (!signal.valid) return false;

    auto local_output = [&](NodeId node) -> std::optional<std::uint32_t> {
        if (network.node_to_scc[node] != scc) return std::nullopt;
        return network.node_position_in_scc[node];
    };
    auto add_signal_to_row = [&](std::uint32_t row, float scale) {
        affine[row] += scale * signal.constant;
        if (signal.constant != 0.0F) has_affine = true;
        for (const auto& [node, coefficient] : signal.terms) {
            add_matrix_term(matrix, affine, network, scc, dimension, row, node,
                            scale * coefficient, has_external, has_affine);
        }
    };

    if (reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) {
        if (reaction.output_count != 1U) return false;
        const NodeId output = network.reaction_outputs[reaction.output_offset];
        const auto row = local_output(output);
        if (!row.has_value()) return false;
        const float gain = reaction.k * network.output_stoichiometry[reaction.output_offset];
        add_signal_to_row(*row, gain);
        if (!matrix.empty()) matrix[static_cast<std::size_t>(*row) * dimension + *row] -= gain;
        return true;
    }

    if (reaction.primitive == Primitive::Produce) {
        for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
            const NodeId output = network.reaction_outputs[reaction.output_offset + index];
            const auto row = local_output(output);
            if (row.has_value()) {
                add_signal_to_row(*row,
                                  reaction.k * network.output_stoichiometry[
                                      reaction.output_offset + index]);
            }
        }
        return true;
    }

    if (reaction.primitive == Primitive::Destroy || reaction.primitive == Primitive::Transfer) {
        if (reaction.input_count != 1U || signal.constant != 0.0F || signal.terms.size() != 1U) {
            return false;
        }
        const NodeId input = network.reaction_inputs[reaction.input_offset];
        if (signal.terms.front().first != input) return false;
        const float input_stoich = network.input_stoichiometry[reaction.input_offset];
        const float signal_coefficient = signal.terms.front().second;
        const float rate = std::min(reaction.k * signal_coefficient, 1.0F / (dt * input_stoich));
        const auto input_local = local_output(input);
        if (input_local.has_value() && !matrix.empty()) {
            matrix[static_cast<std::size_t>(*input_local) * dimension + *input_local] -=
                input_stoich * rate;
        }
        if (reaction.primitive == Primitive::Transfer) {
            for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                const NodeId output = network.reaction_outputs[reaction.output_offset + index];
                const auto row = local_output(output);
                if (row.has_value()) {
                    const float coefficient =
                        network.output_stoichiometry[reaction.output_offset + index] * rate;
                    if (input_local.has_value()) {
                        if (!matrix.empty()) matrix[static_cast<std::size_t>(*row) * dimension + *input_local] +=
                            coefficient;
                    } else {
                        has_external = true;
                        has_affine = true;
                        affine[*row] += coefficient * network.nodes[input].initial_value;
                    }
                }
            }
        }
        return true;
    }
    return false;
}

std::vector<float> matrix_multiply(const std::vector<float>& lhs,
                                   const std::vector<float>& rhs,
                                   std::uint32_t dimension) {
    std::vector<float> result(static_cast<std::size_t>(dimension) * dimension, 0.0F);
    for (std::uint32_t row = 0; row < dimension; ++row) {
        for (std::uint32_t middle = 0; middle < dimension; ++middle) {
            const float left = lhs[static_cast<std::size_t>(row) * dimension + middle];
            if (left == 0.0F) continue;
            for (std::uint32_t column = 0; column < dimension; ++column) {
                result[static_cast<std::size_t>(row) * dimension + column] +=
                    left * rhs[static_cast<std::size_t>(middle) * dimension + column];
            }
        }
    }
    return result;
}

std::vector<float> matrix_add(const std::vector<float>& lhs,
                              const std::vector<float>& rhs) {
    std::vector<float> result(lhs.size(), 0.0F);
    for (std::size_t index = 0; index < lhs.size(); ++index) {
        result[index] = lhs[index] + rhs[index];
    }
    return result;
}

void build_temporal_plans(CompiledNetwork& network, const TemporalCompileConfig& config) {
    network.temporal_compile_dt = config.dt;
    network.temporal_compile_steps = config.horizon_steps;
    network.temporal_plans.assign(network.sccs.size(), SccTemporalPlan{});
    network.temporal_propagators.clear();
    network.temporal_integrals.clear();
    network.temporal_affine.clear();

    constexpr std::uint32_t max_matrix_dimension = 64U;
    constexpr float tolerance = 1.0e-5F;
    for (std::size_t scc_index = 0; scc_index < network.sccs.size(); ++scc_index) {
        const SccId scc = static_cast<SccId>(scc_index);
        const SccBlock& block = network.sccs[scc];
        auto& plan = network.temporal_plans[scc_index];
        plan.node_count = block.node_count;

        bool all_zero = true;
        bool all_linear = true;
        bool all_simple = true;
        bool has_external = false;
        bool has_affine = false;
        std::vector<float> matrix(block.node_count <= max_matrix_dimension
                                      ? static_cast<std::size_t>(block.node_count) * block.node_count : 0U,
                                  0.0F);
        std::vector<float> affine(block.node_count, 0.0F);
        for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
            const Reaction& reaction =
                network.reactions[network.scc_reactions[block.reaction_offset + slot]];
            all_zero = all_zero && reaction.k == 0.0F;
            all_simple = all_simple && is_simple_primitive(reaction.primitive);
            if (!add_linear_reaction(network, scc, reaction, config.dt, matrix, affine,
                                     has_external, has_affine)) {
                all_linear = false;
            }
        }

        if (all_zero) {
            plan.classification = SccTemporalClass::Static;
            continue;
        }
        if (all_linear) {
            plan.classification = has_affine ? SccTemporalClass::LinearAffine
                                             : SccTemporalClass::LinearHomogeneous;
        } else if (all_simple) {
            plan.classification = SccTemporalClass::SimpleNonlinear;
        } else {
            plan.classification = SccTemporalClass::GenericNonlinear;
        }
        plan.external_inputs = has_external ? 1U : 0U;
        if (!all_linear || block.node_count == 0U || block.node_count > max_matrix_dimension) {
            continue;
        }

        // A finite-horizon local propagator is exact for this engine's explicit Euler equations.
        // P = I + dt*A and S = dt * sum(P^i), so x_N = P^N*x_0 + S*b.
        const std::uint32_t dimension = block.node_count;
        std::vector<float> step_matrix(matrix.size(), 0.0F);
        for (std::uint32_t index = 0; index < dimension; ++index) {
            step_matrix[static_cast<std::size_t>(index) * dimension + index] = 1.0F;
        }
        for (std::size_t index = 0; index < matrix.size(); ++index) {
            step_matrix[index] += config.dt * matrix[index];
        }
        // Prove invariance at every step, not just at the final horizon.
        bool bounded = true;
        for (std::uint32_t row = 0; row < dimension; ++row) {
            float upper = config.dt * affine[row];
            bounded = bounded && affine[row] >= 0.0F;
            for (std::uint32_t column = 0; column < dimension; ++column) {
                const float value = step_matrix[static_cast<std::size_t>(row) * dimension + column];
                bounded = bounded && value >= 0.0F;
                upper += value;
            }
            bounded = bounded && upper <= 1.0F;
        }
        std::vector<float> accumulated_matrix(step_matrix.size(), 0.0F);
        std::vector<float> accumulated_integral(step_matrix.size(), 0.0F);
        for (std::uint32_t index = 0; index < dimension; ++index) {
            accumulated_matrix[static_cast<std::size_t>(index) * dimension + index] = 1.0F;
        }
        std::vector<float> base_integral(step_matrix.size(), 0.0F);
        for (std::uint32_t index = 0; index < dimension; ++index) {
            base_integral[static_cast<std::size_t>(index) * dimension + index] = config.dt;
        }
        std::uint32_t remaining = config.horizon_steps;
        while (remaining != 0U) {
            if ((remaining & 1U) != 0U) {
                accumulated_integral = matrix_add(
                    matrix_multiply(step_matrix, accumulated_integral, dimension), base_integral);
                accumulated_matrix = matrix_multiply(step_matrix, accumulated_matrix, dimension);
            }
            const std::vector<float> doubled_integral = matrix_add(
                matrix_multiply(step_matrix, base_integral, dimension), base_integral);
            base_integral = doubled_integral;
            step_matrix = matrix_multiply(step_matrix, step_matrix, dimension);
            remaining >>= 1U;
        }

        const std::vector<float> propagated_affine = [&]() {
            std::vector<float> result(dimension, 0.0F);
            for (std::uint32_t row = 0; row < dimension; ++row) {
                for (std::uint32_t column = 0; column < dimension; ++column) {
                    result[row] += accumulated_integral[static_cast<std::size_t>(row) * dimension +
                                                        column] * affine[column];
                }
            }
            return result;
        }();

        for (std::uint32_t row = 0; row < dimension; ++row) {
            float row_sum = 0.0F;
            for (std::uint32_t column = 0; column < dimension; ++column) {
                const float value = accumulated_matrix[static_cast<std::size_t>(row) * dimension +
                                                       column];
                bounded = bounded && value >= -tolerance;
                row_sum += value;
            }
            bounded = bounded && propagated_affine[row] >= -tolerance &&
                      row_sum + propagated_affine[row] <= 1.0F + tolerance;
        }
        plan.bounded = bounded ? 1U : 0U;
        plan.propagator = 1U;
        plan.matrix_offset = checked_u32(network.temporal_propagators.size(),
                                         "temporal propagator offset");
        plan.matrix_count = checked_u32(accumulated_matrix.size(), "temporal propagator count");
        plan.integral_offset = checked_u32(network.temporal_integrals.size(),
                                           "temporal integral offset");
        plan.integral_count = checked_u32(accumulated_integral.size(), "temporal integral count");
        plan.affine_offset = checked_u32(network.temporal_affine.size(), "temporal affine offset");
        plan.affine_count = dimension;
        network.temporal_propagators.insert(network.temporal_propagators.end(),
                                            accumulated_matrix.begin(), accumulated_matrix.end());
        network.temporal_integrals.insert(network.temporal_integrals.end(),
                                          accumulated_integral.begin(), accumulated_integral.end());
        network.temporal_affine.insert(network.temporal_affine.end(), propagated_affine.begin(),
                                       propagated_affine.end());
    }
}

}  // namespace

std::size_t CompiledNetwork::estimated_bytes() const noexcept {
    std::size_t bytes = 0U;
    bytes += vector_bytes(nodes);
    bytes += vector_bytes(reactions);
    bytes += vector_bytes(reaction_inputs) + vector_bytes(input_stoichiometry);
    bytes += vector_bytes(reaction_outputs) + vector_bytes(output_stoichiometry);
    bytes += vector_bytes(catalysts);
    bytes += vector_bytes(positive_regulators) + vector_bytes(positive_weights) +
             vector_bytes(positive_K) + vector_bytes(positive_n);
    bytes += vector_bytes(negative_regulators) + vector_bytes(negative_weights) +
             vector_bytes(negative_K) + vector_bytes(negative_n);
    bytes += vector_bytes(dependency_offsets) + vector_bytes(dependency_targets);
    bytes += vector_bytes(node_to_scc) + vector_bytes(node_position_in_scc) +
             vector_bytes(sccs) + vector_bytes(scc_nodes) +
             vector_bytes(scc_reactions) + vector_bytes(downstream_sccs) +
             vector_bytes(topological_sccs) + vector_bytes(reaction_scc_offsets) +
             vector_bytes(reaction_sccs) + vector_bytes(execution_plans) +
             vector_bytes(plan_modulate_reactions) + vector_bytes(plan_modulate_slots) +
             vector_bytes(plan_modulate_output_locals) + vector_bytes(plan_modulate_output_stoich) +
             vector_bytes(plan_transfer_reactions) + vector_bytes(plan_transfer_slots) +
             vector_bytes(plan_produce_reactions) + vector_bytes(plan_produce_slots) +
             vector_bytes(plan_destroy_reactions) + vector_bytes(plan_destroy_slots) +
             vector_bytes(plan_generic_reactions) + vector_bytes(plan_generic_slots) +
             vector_bytes(scc_node_reaction_offsets) + vector_bytes(scc_node_reaction_slots) +
             vector_bytes(temporal_plans) + vector_bytes(temporal_propagators) +
             vector_bytes(temporal_integrals) + vector_bytes(temporal_affine);
    bytes += vector_bytes(linear_plans) + vector_bytes(linear_entries) + vector_bytes(linear_boundaries);
    bytes += vector_bytes(node_external_ids) + vector_bytes(reaction_external_ids) +
             vector_bytes(phenotype_nodes);
    for (const auto& name : node_names) {
        bytes += name.capacity();
    }
    return bytes;
}

CompiledNetwork Compiler::compile(const NetworkSpec& spec, std::uint64_t network_version,
                                  CompileMetrics* metrics,
                                  TemporalCompileConfig temporal) const {
    const auto compile_start = std::chrono::steady_clock::now();
    if (metrics != nullptr) *metrics = CompileMetrics{};
    if (network_version == 0U) {
        throw CompileError("network version must be positive");
    }
    if (!(temporal.dt > 0.0F) || !std::isfinite(temporal.dt) || temporal.horizon_steps == 0U) {
        throw CompileError("temporal compile dt must be finite and positive with a positive horizon");
    }
    if (spec.nodes.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
        throw CompileError("node count exceeds uint32 capacity");
    }
    if (spec.reactions.size() > static_cast<std::size_t>(std::numeric_limits<ReactionId>::max())) {
        throw CompileError("reaction count exceeds uint32 capacity");
    }

    CompiledNetwork network;
    network.network_version = network_version;
    network.nodes.reserve(spec.nodes.size());
    network.node_names.reserve(spec.nodes.size());
    network.node_external_ids.reserve(spec.nodes.size());

    std::unordered_map<std::uint64_t, NodeId> node_ids;
    node_ids.reserve(spec.nodes.size());
    std::unordered_map<std::string, std::uint32_t> entity_ids;
    std::unordered_map<std::string, std::uint32_t> state_ids;

    auto intern = [](const std::string& value,
                     std::unordered_map<std::string, std::uint32_t>& table) -> std::uint32_t {
        const auto found = table.find(value);
        if (found != table.end()) {
            return found->second;
        }
        const std::uint32_t id = static_cast<std::uint32_t>(table.size());
        table.emplace(value, id);
        return id;
    };

    for (std::size_t index = 0; index < spec.nodes.size(); ++index) {
        const auto& source = spec.nodes[index];
        if (!node_ids.emplace(source.external_id, static_cast<NodeId>(index)).second) {
            throw CompileError("duplicate node id: " + std::to_string(source.external_id));
        }
        if (source.name.empty()) {
            throw CompileError("node name must not be empty for id " +
                               std::to_string(source.external_id));
        }
        require_finite(source.initial_value, "node initial value");
        if (source.initial_value < 0.0F || source.initial_value > 1.0F) {
            throw CompileError("node initial value must be in [0,1] for id " +
                               std::to_string(source.external_id));
        }
        Node node;
        node.id = static_cast<NodeId>(index);
        node.type = source.type;
        node.compartment = source.compartment;
        node.initial_value = source.initial_value;
        node.entity_id = intern(source.entity.empty() ? source.name : source.entity, entity_ids);
        node.state_id = intern(source.state, state_ids);
        network.nodes.push_back(node);
        network.node_names.push_back(source.name);
        network.node_external_ids.push_back(source.external_id);
        if (source.type == NodeType::Phenotype) {
            network.phenotype_nodes.push_back(node.id);
        }
    }

    std::unordered_map<std::uint64_t, ReactionId> reaction_ids;
    reaction_ids.reserve(spec.reactions.size());
    network.reactions.reserve(spec.reactions.size());
    network.reaction_external_ids.reserve(spec.reactions.size());

    auto resolve_node = [&node_ids](std::uint64_t external_id, const char* role) -> NodeId {
        const auto found = node_ids.find(external_id);
        if (found == node_ids.end()) {
            throw CompileError(std::string("unknown ") + role + " node id: " +
                               std::to_string(external_id));
        }
        return found->second;
    };

    for (std::size_t index = 0; index < spec.reactions.size(); ++index) {
        const auto& source = spec.reactions[index];
        if (!reaction_ids.emplace(source.external_id, static_cast<ReactionId>(index)).second) {
            throw CompileError("duplicate reaction id: " + std::to_string(source.external_id));
        }
        require_finite(source.k, "reaction k");
        require_finite(source.Km, "reaction Km");
        require_finite(source.hill_n, "reaction hill_n");
        require_finite(source.confidence, "reaction confidence");
        require_finite(source.context_weight, "reaction context_weight");
        if (source.k < 0.0F || source.Km < 0.0F || source.hill_n <= 0.0F) {
            throw CompileError("reaction kinetic parameters are out of range for id " +
                               std::to_string(source.external_id));
        }
        if (source.confidence < 0.0F || source.confidence > 1.0F ||
            source.context_weight < 0.0F || source.context_weight > 1.0F) {
            throw CompileError("reaction confidence/context weight must be in [0,1] for id " +
                               std::to_string(source.external_id));
        }
        if (source.outputs.empty() && source.primitive != Primitive::Destroy) {
            throw CompileError("reaction requires an output for id " +
                               std::to_string(source.external_id));
        }
        if (source.inputs.empty() &&
            (source.primitive == Primitive::Transfer || source.primitive == Primitive::Destroy ||
             source.primitive == Primitive::Bind || source.primitive == Primitive::Unbind ||
             source.primitive == Primitive::Move)) {
            throw CompileError("reaction requires an input for id " +
                               std::to_string(source.external_id));
        }

        Reaction reaction;
        reaction.id = static_cast<ReactionId>(index);
        reaction.primitive = source.primitive;
        reaction.input_op = source.input_op;
        reaction.pos_reg_op = source.pos_reg_op;
        reaction.neg_reg_op = source.neg_reg_op;
        reaction.k = source.k;
        reaction.Km = source.Km;
        reaction.hill_n = source.hill_n;
        reaction.confidence = source.confidence;
        reaction.context_weight = source.context_weight;

        reaction.input_offset = checked_u32(network.reaction_inputs.size(), "input offset");
        reaction.input_count = checked_u32(source.inputs.size(), "input count");
        for (const auto& input : source.inputs) {
            require_finite(input.coefficient, "input stoichiometry");
            if (input.coefficient <= 0.0F) {
                throw CompileError("input stoichiometry must be positive");
            }
            network.reaction_inputs.push_back(resolve_node(input.node_external_id, "input"));
            network.input_stoichiometry.push_back(input.coefficient);
        }

        reaction.output_offset = checked_u32(network.reaction_outputs.size(), "output offset");
        reaction.output_count = checked_u32(source.outputs.size(), "output count");
        for (const auto& output : source.outputs) {
            require_finite(output.coefficient, "output stoichiometry");
            if (output.coefficient <= 0.0F) {
                throw CompileError("output stoichiometry must be positive");
            }
            network.reaction_outputs.push_back(resolve_node(output.node_external_id, "output"));
            network.output_stoichiometry.push_back(output.coefficient);
        }

        reaction.catalyst_offset = checked_u32(network.catalysts.size(), "catalyst offset");
        reaction.catalyst_count = checked_u32(source.catalysts.size(), "catalyst count");
        for (const auto catalyst : source.catalysts) {
            network.catalysts.push_back(resolve_node(catalyst, "catalyst"));
        }

        reaction.pos_reg_offset = checked_u32(network.positive_regulators.size(), "positive offset");
        reaction.pos_reg_count = checked_u32(source.positive_regulators.size(), "positive count");
        for (const auto& regulator : source.positive_regulators) {
            require_finite(regulator.weight, "positive regulator weight");
            require_finite(regulator.K, "positive regulator K");
            require_finite(regulator.n, "positive regulator n");
            if (regulator.K < 0.0F || regulator.n <= 0.0F) {
                throw CompileError("positive regulator parameters are out of range");
            }
            network.positive_regulators.push_back(
                resolve_node(regulator.node_external_id, "positive regulator"));
            network.positive_weights.push_back(regulator.weight);
            network.positive_K.push_back(regulator.K);
            network.positive_n.push_back(regulator.n);
        }

        reaction.neg_reg_offset = checked_u32(network.negative_regulators.size(), "negative offset");
        reaction.neg_reg_count = checked_u32(source.negative_regulators.size(), "negative count");
        for (const auto& regulator : source.negative_regulators) {
            require_finite(regulator.weight, "negative regulator weight");
            require_finite(regulator.K, "negative regulator K");
            require_finite(regulator.n, "negative regulator n");
            if (regulator.K < 0.0F || regulator.n <= 0.0F) {
                throw CompileError("negative regulator parameters are out of range");
            }
            network.negative_regulators.push_back(
                resolve_node(regulator.node_external_id, "negative regulator"));
            network.negative_weights.push_back(regulator.weight);
            network.negative_K.push_back(regulator.K);
            network.negative_n.push_back(regulator.n);
        }

        network.reactions.push_back(reaction);
        network.reaction_external_ids.push_back(source.external_id);
    }
    const auto reaction_lowering_finish = std::chrono::steady_clock::now();
    if (metrics != nullptr) {
        metrics->reaction_lowering_ms =
            std::chrono::duration<double, std::milli>(reaction_lowering_finish - compile_start).count();
    }

    // Build a sorted edge list first. This is temporary compiler memory, not runtime layout.
    std::vector<std::pair<NodeId, NodeId>> dependency_edges;
    std::vector<std::vector<NodeId>> reaction_affected(spec.reactions.size());
    for (std::size_t reaction_index = 0; reaction_index < spec.reactions.size(); ++reaction_index) {
        const auto& source = spec.reactions[reaction_index];
        auto affected = affected_nodes(source, node_ids);
        reaction_affected[reaction_index] = affected;

        std::vector<NodeId> dependencies;
        dependencies.reserve(source.inputs.size() + source.catalysts.size() +
                             source.positive_regulators.size() +
                             source.negative_regulators.size());
        for (const auto& input : source.inputs) {
            dependencies.push_back(resolve_node(input.node_external_id, "input"));
        }
        for (const auto catalyst : source.catalysts) {
            dependencies.push_back(resolve_node(catalyst, "catalyst"));
        }
        for (const auto& regulator : source.positive_regulators) {
            dependencies.push_back(resolve_node(regulator.node_external_id, "positive regulator"));
        }
        for (const auto& regulator : source.negative_regulators) {
            dependencies.push_back(resolve_node(regulator.node_external_id, "negative regulator"));
        }
        std::sort(dependencies.begin(), dependencies.end());
        dependencies.erase(std::unique(dependencies.begin(), dependencies.end()), dependencies.end());
        for (const NodeId dependency : dependencies) {
            for (const NodeId changed : affected) {
                dependency_edges.emplace_back(dependency, changed);
            }
        }
    }
    std::sort(dependency_edges.begin(), dependency_edges.end());
    dependency_edges.erase(std::unique(dependency_edges.begin(), dependency_edges.end()),
                           dependency_edges.end());

    network.dependency_offsets.assign(network.nodes.size() + 1U, 0U);
    for (const auto& edge : dependency_edges) {
        ++network.dependency_offsets[static_cast<std::size_t>(edge.first) + 1U];
    }
    for (std::size_t index = 1; index < network.dependency_offsets.size(); ++index) {
        network.dependency_offsets[index] += network.dependency_offsets[index - 1U];
    }
    network.dependency_targets.reserve(dependency_edges.size());
    for (const auto& edge : dependency_edges) {
        network.dependency_targets.push_back(edge.second);
    }
    const auto dependency_finish = std::chrono::steady_clock::now();
    if (metrics != nullptr) {
        metrics->dependency_csr_ms =
            std::chrono::duration<double, std::milli>(dependency_finish - reaction_lowering_finish)
                .count();
    }

    const SccResult scc_result = tarjan_scc(network.dependency_offsets, network.dependency_targets);
    network.node_to_scc = scc_result.component_of;
    network.node_position_in_scc.resize(network.nodes.size(), 0U);
    network.sccs.resize(scc_result.components.size());
    network.scc_nodes.reserve(network.nodes.size());

    std::vector<std::uint8_t> self_loop(network.nodes.size(), 0U);
    for (const auto& edge : dependency_edges) {
        if (edge.first == edge.second) {
            self_loop[edge.first] = 1U;
        }
    }

    for (std::size_t scc_index = 0; scc_index < scc_result.components.size(); ++scc_index) {
        auto& block = network.sccs[scc_index];
        block.id = static_cast<SccId>(scc_index);
        block.node_offset = checked_u32(network.scc_nodes.size(), "SCC node offset");
        block.node_count = checked_u32(scc_result.components[scc_index].size(), "SCC node count");
        block.cyclic = block.node_count > 1U ||
                       (block.node_count == 1U &&
                        self_loop[scc_result.components[scc_index].front()] != 0U);
        network.scc_nodes.insert(network.scc_nodes.end(), scc_result.components[scc_index].begin(),
                                 scc_result.components[scc_index].end());
        for (std::size_t local = 0; local < scc_result.components[scc_index].size(); ++local) {
            network.node_position_in_scc[scc_result.components[scc_index][local]] =
                static_cast<std::uint32_t>(local);
        }
    }

    std::vector<std::pair<SccId, SccId>> dag_edges;
    dag_edges.reserve(dependency_edges.size());
    for (const auto& edge : dependency_edges) {
        const SccId from = network.node_to_scc[edge.first];
        const SccId to = network.node_to_scc[edge.second];
        if (from != to) {
            dag_edges.emplace_back(from, to);
        }
    }
    std::sort(dag_edges.begin(), dag_edges.end());
    dag_edges.erase(std::unique(dag_edges.begin(), dag_edges.end()), dag_edges.end());

    std::vector<std::pair<SccId, ReactionId>> scc_reaction_pairs;
    for (std::size_t reaction_index = 0; reaction_index < reaction_affected.size(); ++reaction_index) {
        for (const NodeId node : reaction_affected[reaction_index]) {
            scc_reaction_pairs.emplace_back(network.node_to_scc[node],
                                            static_cast<ReactionId>(reaction_index));
        }
    }
    std::sort(scc_reaction_pairs.begin(), scc_reaction_pairs.end());
    scc_reaction_pairs.erase(std::unique(scc_reaction_pairs.begin(), scc_reaction_pairs.end()),
                             scc_reaction_pairs.end());

    std::vector<std::pair<ReactionId, SccId>> reaction_scc_pairs;
    reaction_scc_pairs.reserve(scc_reaction_pairs.size());
    for (const auto& pair : scc_reaction_pairs) {
        reaction_scc_pairs.emplace_back(pair.second, pair.first);
    }
    std::sort(reaction_scc_pairs.begin(), reaction_scc_pairs.end());
    network.reaction_scc_offsets.assign(network.reactions.size() + 1U, 0U);
    for (const auto& pair : reaction_scc_pairs) {
        ++network.reaction_scc_offsets[static_cast<std::size_t>(pair.first) + 1U];
        network.reaction_sccs.push_back(pair.second);
    }
    for (std::size_t index = 1; index < network.reaction_scc_offsets.size(); ++index) {
        network.reaction_scc_offsets[index] += network.reaction_scc_offsets[index - 1U];
    }

    std::size_t reaction_pair_index = 0U;
    std::size_t dag_edge_index = 0U;
    for (std::size_t scc_index = 0; scc_index < network.sccs.size(); ++scc_index) {
        auto& block = network.sccs[scc_index];
        block.reaction_offset = checked_u32(network.scc_reactions.size(), "SCC reaction offset");
        while (reaction_pair_index < scc_reaction_pairs.size() &&
               scc_reaction_pairs[reaction_pair_index].first == scc_index) {
            network.scc_reactions.push_back(scc_reaction_pairs[reaction_pair_index].second);
            ++reaction_pair_index;
        }
        block.reaction_count = checked_u32(
            network.scc_reactions.size() - block.reaction_offset, "SCC reaction count");

        block.downstream_offset = checked_u32(network.downstream_sccs.size(), "downstream offset");
        while (dag_edge_index < dag_edges.size() && dag_edges[dag_edge_index].first == scc_index) {
            network.downstream_sccs.push_back(dag_edges[dag_edge_index].second);
            ++dag_edge_index;
        }
        block.downstream_count = checked_u32(
            network.downstream_sccs.size() - block.downstream_offset, "downstream count");
    }

    std::vector<std::uint32_t> indegree(network.sccs.size(), 0U);
    for (const auto& edge : dag_edges) {
        ++indegree[edge.second];
    }
    std::priority_queue<SccId, std::vector<SccId>, std::greater<>> ready;
    for (std::size_t index = 0; index < indegree.size(); ++index) {
        if (indegree[index] == 0U) {
            ready.push(static_cast<SccId>(index));
        }
    }
    while (!ready.empty()) {
        const SccId current = ready.top();
        ready.pop();
        network.topological_sccs.push_back(current);
        const auto& block = network.sccs[current];
        for (std::uint32_t offset = 0; offset < block.downstream_count; ++offset) {
            const SccId next = network.downstream_sccs[block.downstream_offset + offset];
            if (--indegree[next] == 0U) {
                ready.push(next);
            }
        }
    }
    if (network.topological_sccs.size() != network.sccs.size()) {
        throw CompileError("internal error: SCC condensation graph is cyclic");
    }

    const auto scc_finish = std::chrono::steady_clock::now();
    if (metrics != nullptr) {
        metrics->scc_condensation_ms =
            std::chrono::duration<double, std::milli>(scc_finish - dependency_finish).count();
    }

    const auto plan_start = std::chrono::steady_clock::now();
    build_execution_plans(network);
    const auto plan_finish = std::chrono::steady_clock::now();
    if (metrics != nullptr) {
        metrics->execution_plan_ms =
            std::chrono::duration<double, std::milli>(plan_finish - plan_start).count();
    }

    const auto temporal_start = std::chrono::steady_clock::now();
    build_temporal_plans(network, temporal);
    compile_linear_execution(network);
    const auto temporal_finish = std::chrono::steady_clock::now();
    if (metrics != nullptr) {
        metrics->temporal_plan_ms =
            std::chrono::duration<double, std::milli>(temporal_finish - temporal_start).count();
        metrics->total_ms =
            std::chrono::duration<double, std::milli>(temporal_finish - compile_start).count();
    }

    return network;
}

}  // namespace cellnet
