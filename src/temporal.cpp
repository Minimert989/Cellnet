#include "cellnet/network.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cellnet {
namespace {
bool linear_signal(const CompiledNetwork& n, const Reaction& r,
                   std::vector<std::pair<NodeId, float>>& terms, float& constant) {
    if (r.catalyst_count || r.pos_reg_count || r.neg_reg_count) return false;
    constant = r.input_count == 0 ? 1.0F : 0.0F;
    if (r.input_count == 0) return true;
    auto add = [&](std::uint32_t i, float w) {
        terms.emplace_back(n.reaction_inputs[r.input_offset + i], w);
    };
    switch (r.input_op) {
        case RegulationOp::None:
        case RegulationOp::ContinuousOr:
            if (r.input_count != 1 || n.input_stoichiometry[r.input_offset] > 1.0F) return false;
            add(0, n.input_stoichiometry[r.input_offset]); return true;
        case RegulationOp::Product:
        case RegulationOp::ContinuousAnd:
            if (r.input_count != 1 || n.input_stoichiometry[r.input_offset] != 1.0F) return false;
            add(0, 1.0F); return true;
        case RegulationOp::Sum:
            for (std::uint32_t i = 0; i < r.input_count; ++i)
                add(i, n.input_stoichiometry[r.input_offset + i]);
            return true;
        default: return false;
    }
}
}

void compile_linear_execution(CompiledNetwork& n) {
    n.linear_plans.assign(n.sccs.size(), LinearExecutionPlan{});
    n.linear_entries.clear(); n.linear_boundaries.clear();
    for (const auto& block : n.sccs) {
        auto& plan = n.linear_plans[block.id];
        plan.entry_offset = static_cast<std::uint32_t>(n.linear_entries.size());
        plan.boundary_offset = static_cast<std::uint32_t>(n.linear_boundaries.size());
        const auto cls = n.temporal_plans[block.id].classification;
        if (cls == SccTemporalClass::Static) { plan.backend = TemporalBackend::StaticDirect; continue; }
        if (cls != SccTemporalClass::LinearHomogeneous && cls != SccTemporalClass::LinearAffine) {
            plan.backend = TemporalBackend::NumericalNonlinear; continue;
        }
        bool supported = true;
        std::vector<NodeId> boundaries;
        for (std::uint32_t i = 0; i < block.reaction_count && supported; ++i) {
            const auto& r = n.reactions[n.scc_reactions[block.reaction_offset + i]];
            std::vector<std::pair<NodeId, float>> signal;
            float constant = 0;
            supported = linear_signal(n, r, signal, constant);
            if (!supported) break;
            const auto signal_scale = r.context_weight * r.confidence;
            if (r.primitive == Primitive::Modulate || r.primitive == Primitive::Delay) {
                if (r.output_count != 1 || n.node_to_scc[n.reaction_outputs[r.output_offset]] != block.id) {
                    supported = false; break;
                }
                const NodeId row = n.reaction_outputs[r.output_offset];
                const float gain = r.k * n.output_stoichiometry[r.output_offset];
                // The reference relaxes by k*(signal - state): an inhibition scales the
                // signal contribution, while the -k*state diagonal remains fixed.
                n.linear_entries.push_back({row, row, r.id, -gain, 1, 1, LinearCoefficientKind::Fixed});
                for (const auto& [node, w] : signal) {
                    n.linear_entries.push_back({row, node, r.id, gain * w, 1, 1, LinearCoefficientKind::Signal});
                    if (n.node_to_scc[node] != block.id) boundaries.push_back(node);
                }
                if (constant != 0) n.linear_entries.push_back({row, kInvalidNode, r.id,
                    gain * constant, 1, 1, LinearCoefficientKind::Signal});
            } else if (r.primitive == Primitive::Produce) {
                for (std::uint32_t j = 0; j < r.output_count; ++j) {
                    const NodeId row = n.reaction_outputs[r.output_offset + j];
                    if (n.node_to_scc[row] != block.id) { supported = false; break; }
                    const float gain = r.k * n.output_stoichiometry[r.output_offset + j];
                    for (const auto& [node, w] : signal) {
                        n.linear_entries.push_back({row, node, r.id, gain * w, 1, 1, LinearCoefficientKind::Signal});
                        if (n.node_to_scc[node] != block.id) boundaries.push_back(node);
                    }
                    if (constant != 0) n.linear_entries.push_back({row, kInvalidNode, r.id,
                        gain * constant, 1, 1, LinearCoefficientKind::Signal});
                }
            } else if ((r.primitive == Primitive::Destroy || r.primitive == Primitive::Transfer) &&
                       r.input_count == 1 && signal.size() == 1 && signal[0].first == n.reaction_inputs[r.input_offset] &&
                       constant == 0) {
                const NodeId input = n.reaction_inputs[r.input_offset];
                const float in_st = n.input_stoichiometry[r.input_offset];
                const float w = signal[0].second;
                // A transfer reaction can affect the reactant and product in different
                // SCCs. Emit only rows owned by this SCC; the shared reaction is lowered
                // again in each affected SCC so the fused-region compiler can join the
                // pieces without placing an external row in the wrong region.
                if (n.node_to_scc[input] == block.id) {
                    n.linear_entries.push_back({input, input, r.id, -r.k * w * in_st,
                                                1.0F / n.temporal_compile_dt, 1,
                                                LinearCoefficientKind::Capped});
                }
                if (r.primitive == Primitive::Transfer) for (std::uint32_t j = 0; j < r.output_count; ++j) {
                    const NodeId row = n.reaction_outputs[r.output_offset + j];
                    if (n.node_to_scc[row] != block.id) continue;
                    n.linear_entries.push_back({row, input, r.id,
                        r.k * w * n.output_stoichiometry[r.output_offset + j],
                        n.output_stoichiometry[r.output_offset + j] /
                            (n.temporal_compile_dt * in_st), 1,
                        LinearCoefficientKind::Capped});
                }
            } else supported = false;
            for (const auto& [node, _] : signal)
                if (n.node_to_scc[node] != block.id) boundaries.push_back(node);
            (void)signal_scale;
        }
        if (!supported) {
            n.linear_entries.resize(plan.entry_offset);
            plan.backend = TemporalBackend::NumericalNonlinear;
            continue;
        }
        std::sort(boundaries.begin(), boundaries.end());
        boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());
        n.linear_boundaries.insert(n.linear_boundaries.end(), boundaries.begin(), boundaries.end());
        plan.entry_count = static_cast<std::uint32_t>(n.linear_entries.size() - plan.entry_offset);
        plan.boundary_count = static_cast<std::uint32_t>(boundaries.size());
        const auto pairs = static_cast<std::uint64_t>(plan.entry_count);
        plan.density = block.node_count == 0 ? 0.0 :
            std::min(1.0, static_cast<double>(pairs) / (static_cast<double>(block.node_count) * block.node_count));
        if (!boundaries.empty()) plan.backend = TemporalBackend::BoundaryDrivenLinear;
        else if (block.node_count == 1) plan.backend = TemporalBackend::ScalarClosedForm;
        else if (block.node_count <= 64) plan.backend = TemporalBackend::PrecomputedDensePropagator;
        else if (block.node_count <= 512 && plan.density >= 0.15)
            plan.backend = TemporalBackend::RuntimeDensePropagator;
        else plan.backend = TemporalBackend::SparseLinear;
    }
}
} // namespace cellnet
