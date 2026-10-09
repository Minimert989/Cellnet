#include "cellnet/engine.hpp"

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cellnet {
namespace {

float clamp01(float value) {
    return std::clamp(value, 0.0F, 1.0F);
}

float hill_positive(float value, float K, float n) {
    const float x = std::max(value, 0.0F);
    const float safe_K = std::max(K, 1.0e-12F);
    const float xn = std::pow(x, n);
    const float kn = std::pow(safe_K, n);
    const float denominator = kn + xn;
    return denominator > 0.0F ? xn / denominator : 0.0F;
}

float hill_negative(float value, float K, float n) {
    return 1.0F - hill_positive(value, K, n);
}

bool consumes_inputs(Primitive primitive) {
    return primitive == Primitive::Transfer || primitive == Primitive::Destroy ||
           primitive == Primitive::Bind || primitive == Primitive::Unbind ||
           primitive == Primitive::Move;
}

bool changes_outputs(Primitive primitive) {
    return primitive != Primitive::Destroy;
}

template <class T>
std::size_t vector_bytes(const std::vector<T>& values) noexcept {
    return values.capacity() * sizeof(T);
}

}  // namespace

Engine::Engine(CompiledNetwork network, SolverConfig config)
    : network_(std::move(network)), config_(config) {
    if (!(config_.dt > 0.0F) || !std::isfinite(config_.dt)) {
        throw std::invalid_argument("solver dt must be finite and positive");
    }
    if (!(config_.epsilon > 0.0F) || !std::isfinite(config_.epsilon)) {
        throw std::invalid_argument("solver epsilon must be finite and positive");
    }
    if (!(config_.result_epsilon >= config_.epsilon) || !std::isfinite(config_.result_epsilon)) {
        throw std::invalid_argument("result epsilon must be finite and >= solver epsilon");
    }
    if (config_.baseline_max_steps == 0U || config_.perturbation_steps == 0U) {
        throw std::invalid_argument("solver step limits must be positive");
    }

    state_.resize(network_.nodes.size(), 0.0F);
    protocol_saved_state_.resize(network_.nodes.size(), 0.0F);
    protocol_saved_baseline_.resize(network_.nodes.size(), 0.0F);
    reaction_multiplier_.resize(network_.reactions.size(), 1.0F);
    protocol_saved_reaction_multiplier_.resize(network_.reactions.size(), 1.0F);
    context_reaction_multiplier_.resize(network_.reactions.size(), 1.0F);
    context_node_multiplier_.resize(network_.nodes.size(), 1.0F);
    extra_node_decay_.resize(network_.nodes.size(), 0.0F);
    protocol_saved_extra_node_decay_.resize(network_.nodes.size(), 0.0F);
    clamp_value_.resize(network_.nodes.size(), 0.0F);
    protocol_saved_clamp_value_.resize(network_.nodes.size(), 0.0F);
    clamp_active_.resize(network_.nodes.size(), 0U);
    protocol_saved_clamp_active_.resize(network_.nodes.size(), 0U);
    journal_generation_.resize(network_.nodes.size(), 0U);
    queue_generation_.resize(network_.sccs.size(), 0U);
    dirty_counts_.resize(network_.sccs.size(), 0U);
    dirty_count_generation_.resize(network_.sccs.size(), 0U);
    dirty_reaction_generation_.resize(network_.scc_reactions.size(), 0U);
    scc_visit_generation_.resize(network_.sccs.size(), 0U);
    analytic_visit_generation_.resize(network_.sccs.size(), 0U);
    numerical_visit_generation_.resize(network_.sccs.size(), 0U);
    temporal_completed_generation_.resize(network_.sccs.size(), 0U);
    initialize_linear_workspace();

    std::size_t largest_scc = 0U;
    for (const auto& block : network_.sccs) {
        largest_scc = std::max(largest_scc, static_cast<std::size_t>(block.node_count));
    }
    scratch_values_.resize(largest_scc, 0.0F);
    scratch_flags_.resize(largest_scc, 0U);
    journal_.reserve(network_.nodes.size());
    pending_.reserve(network_.nodes.size());
    changed_nodes_this_step_.reserve(network_.nodes.size());
    current_frontier_.reserve(network_.sccs.size());
    next_frontier_.reserve(network_.sccs.size());
    seed_sccs_.reserve(network_.sccs.size());
    protocol_seeds_.reserve(network_.sccs.size());
    protocol_causal_queue_.reserve(network_.sccs.size());
    protocol_topological_rank_.resize(network_.sccs.size(), 0U);
    protocol_seed_generation_.resize(network_.sccs.size(), 0U);
    for (std::uint32_t rank = 0; rank < network_.topological_sccs.size(); ++rank)
        protocol_topological_rank_[network_.topological_sccs[rank]] = rank;
    active_deltas_.reserve(network_.sccs.size());
}

EngineMemoryBreakdown Engine::memory_breakdown() const noexcept {
    EngineMemoryBreakdown result;
    const auto& n = network_;
    result.base_graph_bytes = vector_bytes(n.nodes) + vector_bytes(n.reactions) +
        vector_bytes(n.reaction_inputs) + vector_bytes(n.input_stoichiometry) +
        vector_bytes(n.reaction_outputs) + vector_bytes(n.output_stoichiometry) +
        vector_bytes(n.catalysts) + vector_bytes(n.positive_regulators) +
        vector_bytes(n.positive_weights) + vector_bytes(n.positive_K) + vector_bytes(n.positive_n) +
        vector_bytes(n.negative_regulators) + vector_bytes(n.negative_weights) +
        vector_bytes(n.negative_K) + vector_bytes(n.negative_n) +
        vector_bytes(n.node_external_ids) + vector_bytes(n.reaction_external_ids) +
        vector_bytes(n.phenotype_nodes);
    for (const auto& name : n.node_names) result.base_graph_bytes += name.capacity();
    result.base_graph_bytes += vector_bytes(n.node_names);

    result.scc_metadata_bytes = vector_bytes(n.dependency_offsets) + vector_bytes(n.dependency_targets) +
        vector_bytes(n.node_to_scc) + vector_bytes(n.node_position_in_scc) + vector_bytes(n.sccs) +
        vector_bytes(n.scc_nodes) + vector_bytes(n.scc_reactions) + vector_bytes(n.downstream_sccs) +
        vector_bytes(n.topological_sccs) + vector_bytes(n.reaction_scc_offsets) +
        vector_bytes(n.reaction_sccs) + vector_bytes(n.execution_plans) +
        vector_bytes(n.plan_modulate_reactions) + vector_bytes(n.plan_modulate_slots) +
        vector_bytes(n.plan_modulate_output_locals) + vector_bytes(n.plan_modulate_output_stoich) +
        vector_bytes(n.plan_transfer_reactions) + vector_bytes(n.plan_transfer_slots) +
        vector_bytes(n.plan_produce_reactions) + vector_bytes(n.plan_produce_slots) +
        vector_bytes(n.plan_destroy_reactions) + vector_bytes(n.plan_destroy_slots) +
        vector_bytes(n.plan_generic_reactions) + vector_bytes(n.plan_generic_slots) +
        vector_bytes(n.scc_node_reaction_offsets) + vector_bytes(n.scc_node_reaction_slots);
    result.temporal_plan_bytes = vector_bytes(n.temporal_plans) + vector_bytes(n.temporal_propagators) +
        vector_bytes(n.temporal_integrals) + vector_bytes(n.temporal_affine) +
        vector_bytes(n.linear_plans) + vector_bytes(n.linear_entries) + vector_bytes(n.linear_boundaries);
    result.fused_region_bytes = vector_bytes(linear_regions_) + vector_bytes(linear_region_sccs_) +
        vector_bytes(linear_region_nodes_) + vector_bytes(linear_region_boundaries_) +
        vector_bytes(linear_region_entries_) + vector_bytes(linear_region_of_scc_) +
        vector_bytes(linear_region_node_slot_) + vector_bytes(linear_region_dependency_offsets_) +
        vector_bytes(linear_region_dependency_targets_) + vector_bytes(linear_region_entry_row_offsets_) +
        vector_bytes(linear_region_entry_row_indices_) + vector_bytes(linear_region_baseline_stationary_) +
        vector_bytes(cached_causal_cone_nodes_) +
        vector_bytes(cached_causal_cone_entries_) + vector_bytes(cached_causal_boundary_nodes_);
    result.response_workspace_bytes = vector_bytes(boundary_previous_) +
        vector_bytes(boundary_deferred_scc_) + vector_bytes(boundary_trace_nodes_) +
        vector_bytes(boundary_trace_index_) + vector_bytes(boundary_trajectory_) +
        vector_bytes(boundary_saved_state_) + vector_bytes(boundary_segment_start_) +
        vector_bytes(region_accumulator_) + vector_bytes(region_initial_state_) +
        vector_bytes(region_segment_response_) + vector_bytes(causal_cone_generation_) +
        vector_bytes(causal_cone_depth_) + vector_bytes(causal_cone_queue_) +
        vector_bytes(active_region_node_slots_) + vector_bytes(active_region_entry_slots_) +
        vector_bytes(active_boundary_nodes_) + vector_bytes(causal_cone_scc_generation_) +
        vector_bytes(boundary_previous_generation_) +
        vector_bytes(active_region_sccs_);
    result.scratch_workspace_bytes = vector_bytes(state_) + vector_bytes(baseline_) +
        vector_bytes(protocol_saved_state_) + vector_bytes(protocol_saved_baseline_) +
        vector_bytes(protocol_saved_reaction_multiplier_) + vector_bytes(protocol_saved_extra_node_decay_) +
        vector_bytes(protocol_saved_clamp_value_) + vector_bytes(protocol_saved_clamp_active_) +
        vector_bytes(protocol_seeds_) + vector_bytes(protocol_causal_queue_) +
        vector_bytes(protocol_topological_rank_) + vector_bytes(protocol_seed_generation_) +
        vector_bytes(reaction_multiplier_) + vector_bytes(context_reaction_multiplier_) +
        vector_bytes(context_node_multiplier_) + vector_bytes(extra_node_decay_) +
        vector_bytes(clamp_value_) + vector_bytes(clamp_active_) + vector_bytes(journal_generation_) +
        vector_bytes(queue_generation_) + vector_bytes(journal_) + vector_bytes(pending_) +
        vector_bytes(current_frontier_) + vector_bytes(next_frontier_) + vector_bytes(scratch_values_) +
        vector_bytes(scratch_flags_) + vector_bytes(seed_sccs_) + vector_bytes(active_deltas_) +
        vector_bytes(dirty_counts_) + vector_bytes(dirty_reaction_generation_) +
        vector_bytes(scc_visit_generation_) + vector_bytes(analytic_visit_generation_) +
        vector_bytes(numerical_visit_generation_) + vector_bytes(temporal_completed_generation_) +
        vector_bytes(dirty_count_generation_) +
        vector_bytes(changed_nodes_this_step_) + vector_bytes(linear_coefficients_) +
        vector_bytes(linear_scratch_) + vector_bytes(power_action_a_) + vector_bytes(power_action_b_) +
        vector_bytes(power_action_result_) + vector_bytes(power_action_diagonal_) +
        vector_bytes(power_action_offdiagonal_) + vector_bytes(power_action_rowsum_) +
        vector_bytes(linear_generation_) + vector_bytes(boundary_generation_) +
        vector_bytes(dense_base_) + vector_bytes(dense_power_) + vector_bytes(dense_product_) +
        vector_bytes(dense_free_index_);
    result.cache_bytes = vector_bytes(propagator_cache_) + vector_bytes(power_action_cache_) +
        hot_cache_.bucket_count() * sizeof(void*);
    for (const auto& slot : propagator_cache_) result.cache_bytes += vector_bytes(slot.matrix);
    for (const auto& slot : power_action_cache_) result.cache_bytes += vector_bytes(slot.values);
    for (const auto& [key, value] : hot_cache_) {
        (void)key;
        result.cache_bytes += sizeof(CacheKey) + sizeof(SimulationResult) +
            vector_bytes(value.molecular_changes) + vector_bytes(value.phenotype_changes);
    }
    return result;
}

void Engine::initialize_linear_workspace() {
    linear_scratch_.resize(network_.nodes.size(), 0.0);
    power_action_a_.resize(network_.nodes.size() + 1U, 0.0);
    power_action_b_.resize(network_.nodes.size() + 1U, 0.0);
    power_action_result_.resize(network_.nodes.size() + 1U, 0.0);
    power_action_diagonal_.resize(network_.nodes.size(), 0.0);
    power_action_offdiagonal_.resize(network_.nodes.size(), 0.0);
    power_action_rowsum_.resize(network_.nodes.size(), 0.0);
    region_accumulator_.resize(network_.nodes.size(), 0.0);
    region_initial_state_.resize(network_.nodes.size(), 0.0);
    region_segment_response_.resize(network_.nodes.size(), 0.0);
    boundary_segment_start_.resize(config_.perturbation_steps, 0U);
    boundary_previous_.resize(network_.nodes.size(), 0.0F);
    boundary_previous_generation_.resize(network_.nodes.size(), 0U);
    std::size_t largest = 1U;
    for (const auto& block : network_.sccs) largest = std::max(largest, static_cast<std::size_t>(block.node_count));
    dense_free_index_.resize(largest, std::numeric_limits<std::uint32_t>::max());
    constexpr std::size_t max_dense_dimension = 513U;
    dense_base_.resize(max_dense_dimension * max_dense_dimension);
    dense_power_.resize(max_dense_dimension * max_dense_dimension);
    dense_product_.resize(max_dense_dimension * max_dense_dimension);
    propagator_cache_.resize(4U);
    for (auto& slot : propagator_cache_) slot.matrix.resize(max_dense_dimension * max_dense_dimension);
    power_action_cache_.resize(2U);
    for (auto& slot : power_action_cache_) slot.values.resize(largest);
    linear_coefficients_.resize(network_.reactions.size(), 1.0);

    boundary_deferred_scc_.assign(network_.sccs.size(), 0U);
    boundary_trace_index_.assign(network_.nodes.size(), std::numeric_limits<std::uint32_t>::max());
    linear_region_of_scc_.assign(network_.sccs.size(), std::numeric_limits<std::uint32_t>::max());
    if (config_.mode != SimulationMode::Continuous ||
        config_.dt != network_.temporal_compile_dt ||
        config_.perturbation_steps != network_.temporal_compile_steps ||
        network_.linear_plans.empty()) return;
    const auto region_construction_start = std::chrono::steady_clock::now();
    const std::size_t scc_count = network_.sccs.size();
    std::vector<SccId> parent(scc_count);
    std::vector<std::uint32_t> component_size(scc_count, 1U);
    for (SccId scc = 0; scc < scc_count; ++scc) parent[scc] = scc;
    const auto is_linear = [&](SccId scc) {
        if (network_.linear_plans[scc].backend == TemporalBackend::NumericalNonlinear) return false;
        const auto cls = network_.temporal_plans[scc].classification;
        return cls == SccTemporalClass::LinearAffine || cls == SccTemporalClass::LinearHomogeneous;
    };
    const auto root = [&](SccId scc) {
        SccId representative = scc;
        while (parent[representative] != representative)
            representative = parent[representative];
        while (parent[scc] != scc) {
            const SccId next = parent[scc];
            parent[scc] = representative;
            scc = next;
        }
        return representative;
    };
    for (SccId scc = 0; scc < scc_count; ++scc) {
        if (!is_linear(scc)) continue;
        const auto& block = network_.sccs[scc];
        for (std::uint32_t i = 0; i < block.downstream_count; ++i) {
            const SccId downstream = network_.downstream_sccs[block.downstream_offset + i];
            if (!is_linear(downstream)) continue;
            SccId a = root(scc), b = root(downstream);
            if (a == b) continue;
            if (component_size[a] < component_size[b]) std::swap(a, b);
            parent[b] = a;
            component_size[a] += component_size[b];
        }
    }
    std::vector<SccId> root_to_region(scc_count, kInvalidScc);
    for (SccId scc = 0; scc < scc_count; ++scc) {
        if (!is_linear(scc)) continue;
        const SccId r = root(scc);
        if (root_to_region[r] == kInvalidScc) {
            root_to_region[r] = static_cast<SccId>(linear_regions_.size());
            linear_regions_.push_back({});
        }
        const std::uint32_t region = root_to_region[r];
        linear_region_of_scc_[scc] = region;
        auto& descriptor = linear_regions_[region];
        ++descriptor.scc_count;
        descriptor.node_count += network_.sccs[scc].node_count;
        descriptor.entry_count += network_.linear_plans[scc].entry_count;
    }
    std::size_t total_region_sccs = 0U, total_region_nodes = 0U, total_region_entries = 0U;
    std::vector<std::uint32_t> scc_cursor(linear_regions_.size());
    std::vector<std::uint32_t> node_cursor(linear_regions_.size());
    std::vector<std::uint32_t> entry_cursor(linear_regions_.size());
    for (std::size_t region_id = 0; region_id < linear_regions_.size(); ++region_id) {
        auto& descriptor = linear_regions_[region_id];
        descriptor.scc_offset = static_cast<std::uint32_t>(total_region_sccs);
        descriptor.node_offset = static_cast<std::uint32_t>(total_region_nodes);
        descriptor.entry_offset = static_cast<std::uint32_t>(total_region_entries);
        descriptor.boundary_offset = 0U;
        scc_cursor[region_id] = descriptor.scc_offset;
        node_cursor[region_id] = descriptor.node_offset;
        entry_cursor[region_id] = descriptor.entry_offset;
        total_region_sccs += descriptor.scc_count;
        total_region_nodes += descriptor.node_count;
        total_region_entries += descriptor.entry_count;
    }
    linear_region_sccs_.resize(total_region_sccs);
    linear_region_nodes_.resize(total_region_nodes);
    linear_region_entries_.resize(total_region_entries);
    // Gather each SCC once into pre-counted contiguous region ranges. This avoids rescanning
    // the complete condensation graph once per region (which is quadratic for many regions).
    for (SccId scc = 0; scc < scc_count; ++scc) {
        const std::uint32_t region = linear_region_of_scc_[scc];
        if (region == std::numeric_limits<std::uint32_t>::max()) continue;
        const auto& block = network_.sccs[scc];
        const auto& plan = network_.linear_plans[scc];
        linear_region_sccs_[scc_cursor[region]++] = scc;
        for (std::uint32_t i = 0; i < block.node_count; ++i)
            linear_region_nodes_[node_cursor[region]++] = network_.scc_nodes[block.node_offset + i];
        for (std::uint32_t i = 0; i < plan.entry_count; ++i)
            linear_region_entries_[entry_cursor[region]++] =
                network_.linear_entries[plan.entry_offset + i];
    }
    linear_region_node_slot_.assign(network_.nodes.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::uint32_t slot = 0; slot < linear_region_nodes_.size(); ++slot)
        linear_region_node_slot_[linear_region_nodes_[slot]] = slot;
    linear_region_dependency_offsets_.assign(linear_region_nodes_.size() + 1U, 0U);
    for (const auto& entry : linear_region_entries_) {
        if (entry.column == kInvalidNode) continue;
        const auto source_slot = linear_region_node_slot_[entry.column];
        const auto row_slot = linear_region_node_slot_[entry.row];
        if (source_slot == std::numeric_limits<std::uint32_t>::max() ||
            row_slot == std::numeric_limits<std::uint32_t>::max() ||
            linear_region_of_scc_[network_.node_to_scc[entry.column]] !=
                linear_region_of_scc_[network_.node_to_scc[entry.row]]) continue;
        ++linear_region_dependency_offsets_[source_slot + 1U];
    }
    for (std::size_t slot = 1; slot < linear_region_dependency_offsets_.size(); ++slot)
        linear_region_dependency_offsets_[slot] += linear_region_dependency_offsets_[slot - 1U];
    linear_region_dependency_targets_.resize(linear_region_dependency_offsets_.back());
    std::vector<std::uint64_t> dependency_cursor = linear_region_dependency_offsets_;
    for (const auto& entry : linear_region_entries_) {
        if (entry.column == kInvalidNode) continue;
        const auto source_slot = linear_region_node_slot_[entry.column];
        const auto row_slot = linear_region_node_slot_[entry.row];
        if (source_slot == std::numeric_limits<std::uint32_t>::max() ||
            row_slot == std::numeric_limits<std::uint32_t>::max() ||
            linear_region_of_scc_[network_.node_to_scc[entry.column]] !=
                linear_region_of_scc_[network_.node_to_scc[entry.row]]) continue;
        linear_region_dependency_targets_[dependency_cursor[source_slot]++] = row_slot;
    }
    linear_region_entry_row_offsets_.assign(linear_region_nodes_.size() + 1U, 0U);
    for (std::uint32_t entry_slot = 0; entry_slot < linear_region_entries_.size(); ++entry_slot) {
        const NodeId row = linear_region_entries_[entry_slot].row;
        const std::uint32_t row_slot = linear_region_node_slot_[row];
        if (row_slot == std::numeric_limits<std::uint32_t>::max())
            throw std::logic_error("linear region entry row is missing from its region node list");
        ++linear_region_entry_row_offsets_[row_slot + 1U];
    }
    for (std::size_t slot = 1; slot < linear_region_entry_row_offsets_.size(); ++slot)
        linear_region_entry_row_offsets_[slot] += linear_region_entry_row_offsets_[slot - 1U];
    linear_region_entry_row_indices_.resize(linear_region_entries_.size());
    std::vector<std::uint64_t> row_entry_cursor = linear_region_entry_row_offsets_;
    for (std::uint32_t entry_slot = 0; entry_slot < linear_region_entries_.size(); ++entry_slot) {
        const std::uint32_t row_slot = linear_region_node_slot_[linear_region_entries_[entry_slot].row];
        if (row_slot == std::numeric_limits<std::uint32_t>::max())
            throw std::logic_error("linear region entry row is missing while building entry-row CSR");
        linear_region_entry_row_indices_[row_entry_cursor[row_slot]++] = entry_slot;
    }
    causal_cone_generation_.assign(linear_region_nodes_.size(), 0U);
    causal_cone_depth_.resize(linear_region_nodes_.size(), 0U);
    causal_cone_queue_.reserve(linear_region_nodes_.size());
    active_region_node_slots_.reserve(linear_region_nodes_.size());
    active_region_entry_slots_.reserve(linear_region_entries_.size());
    causal_cone_scc_generation_.assign(network_.sccs.size(), 0U);
    active_region_sccs_.reserve(network_.sccs.size());
    std::vector<std::uint8_t> needed(network_.nodes.size(), 0U);
    std::vector<std::uint32_t> boundary_seen(network_.nodes.size(), std::numeric_limits<std::uint32_t>::max());
    for (std::uint32_t region_id = 0; region_id < linear_regions_.size(); ++region_id) {
        auto& region = linear_regions_[region_id];
        region.boundary_offset = static_cast<std::uint32_t>(linear_region_boundaries_.size());
        bool has_nonlinear_boundary = false;
        bool terminal = true;
        for (std::uint32_t i = 0; i < region.node_count; ++i) {
            const NodeId row = linear_region_nodes_[region.node_offset + i];
            const SccId scc = network_.node_to_scc[row];
            const auto& block = network_.sccs[scc];
            for (std::uint32_t j = 0; j < block.downstream_count; ++j) {
                const SccId downstream = network_.downstream_sccs[block.downstream_offset + j];
                if (linear_region_of_scc_[downstream] != region_id &&
                    network_.temporal_plans[downstream].classification != SccTemporalClass::Static)
                    terminal = false;
            }
        }
        for (std::uint32_t i = 0; i < region.entry_count; ++i) {
            const auto& entry = linear_region_entries_[region.entry_offset + i];
            if (entry.column == kInvalidNode || linear_region_of_scc_[network_.node_to_scc[entry.column]] == region_id)
                continue;
            if (boundary_seen[entry.column] != region_id) {
                boundary_seen[entry.column] = region_id;
                linear_region_boundaries_.push_back(entry.column);
                ++region.boundary_count;
            }
            const auto cls = network_.temporal_plans[network_.node_to_scc[entry.column]].classification;
            has_nonlinear_boundary = has_nonlinear_boundary ||
                cls == SccTemporalClass::SimpleNonlinear || cls == SccTemporalClass::GenericNonlinear;
        }
        // Fuse connected linear SCCs as one temporal state vector. A singleton is
        // deferred only when it has a changing nonlinear input trajectory.
        if (!terminal || (region.scc_count == 1U && !has_nonlinear_boundary)) continue;
        region.deferred = true;
        for (std::uint32_t i = 0; i < region.scc_count; ++i)
            boundary_deferred_scc_[linear_region_sccs_[region.scc_offset + i]] = 1U;
        for (std::uint32_t i = 0; i < region.entry_count; ++i) {
            const auto& entry = linear_region_entries_[region.entry_offset + i];
            if (entry.column != kInvalidNode && linear_region_of_scc_[network_.node_to_scc[entry.column]] != region_id)
                needed[entry.column] = 1U;
        }
    }
    boundary_trace_horizon_ = config_.perturbation_steps;
    for (NodeId node = 0; node < needed.size(); ++node) {
        if (!needed[node]) continue;
        boundary_trace_index_[node] = static_cast<std::uint32_t>(boundary_trace_nodes_.size());
        boundary_trace_nodes_.push_back(node);
    }
    constexpr std::size_t max_trajectory_bytes = 64U * 1024U * 1024U;
    if (static_cast<std::uint64_t>(boundary_trace_nodes_.size()) * boundary_trace_horizon_ * sizeof(float) > max_trajectory_bytes) {
        std::fill(boundary_deferred_scc_.begin(), boundary_deferred_scc_.end(), 0U);
        boundary_trace_nodes_.clear();
        boundary_trace_horizon_ = 0U;
        for (auto& region : linear_regions_) region.deferred = false;
        std::fill(boundary_deferred_scc_.begin(), boundary_deferred_scc_.end(), 0U);
        region_construction_time_us_ = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - region_construction_start).count();
        return;
    }
    boundary_trajectory_.resize(boundary_trace_nodes_.size() * boundary_trace_horizon_);
    boundary_saved_state_.resize(boundary_trace_nodes_.size());
    active_boundary_nodes_.reserve(linear_region_boundaries_.size());
    cached_causal_cone_nodes_.reserve(linear_region_nodes_.size());
    cached_causal_cone_entries_.reserve(linear_region_entries_.size());
    cached_causal_boundary_nodes_.reserve(linear_region_boundaries_.size());
    region_construction_time_us_ = std::chrono::duration<double, std::micro>(
        std::chrono::steady_clock::now() - region_construction_start).count();
}

void Engine::initialize_baseline(const BaselineStepObserver& observer) {
    for (std::size_t index = 0; index < network_.nodes.size(); ++index) {
        state_[index] = clamp01(network_.nodes[index].initial_value * context_node_multiplier_[index]);
    }
    std::fill(reaction_multiplier_.begin(), reaction_multiplier_.end(), 1.0F);
    std::fill(extra_node_decay_.begin(), extra_node_decay_.end(), 0.0F);
    std::fill(clamp_active_.begin(), clamp_active_.end(), 0U);
    journal_.clear();
    journaling_ = false;

    const float baseline_threshold = config_.epsilon * std::min(config_.dt, 1.0F);
    const ExecutionStats stats = run_full(config_.baseline_max_steps, baseline_threshold, observer);
    if (stats.numerical_error) {
        throw std::runtime_error("baseline construction produced NaN or Inf");
    }
    if (!stats.converged) {
        throw std::runtime_error("baseline did not converge within the configured step limit");
    }
    baseline_ = state_;
    linear_region_baseline_stationary_.assign(linear_regions_.size(), 1U);
    bool has_deferred_region = false;
    for (const auto& region : linear_regions_) has_deferred_region = has_deferred_region || region.deferred;
    if (has_deferred_region && config_.perturbation_steps > 0U) {
        // The causal cone contains perturbation influence, not spontaneous baseline drift.
        // Measure one unchanged reference timestep at initialization and use a conservative
        // finite-horizon bound; nonstationary baselines retain full-region analytic execution.
        const ExecutionStats drift = run_full(1U, 0.0F);
        if (drift.numerical_error) throw std::runtime_error("baseline drift check produced NaN or Inf");
        const double allowed_step_delta = 5.0e-6 / config_.perturbation_steps;
        for (std::uint32_t region_id = 0; region_id < linear_regions_.size(); ++region_id) {
            const auto& region = linear_regions_[region_id];
            if (!region.deferred) continue;
            double maximum_delta = 0.0;
            for (std::uint32_t i = 0; i < region.node_count; ++i) {
                const NodeId node = linear_region_nodes_[region.node_offset + i];
                maximum_delta = std::max(maximum_delta,
                    static_cast<double>(std::abs(state_[node] - baseline_[node])));
            }
            double maximum_boundary_delta = 0.0;
            for (std::uint32_t i = 0; i < region.boundary_count; ++i) {
                const NodeId node = linear_region_boundaries_[region.boundary_offset + i];
                maximum_boundary_delta = std::max(maximum_boundary_delta,
                    static_cast<double>(std::abs(state_[node] - baseline_[node])));
            }
            bool contractive = true;
            state_ = baseline_;
            for (std::uint32_t i = 0; i < region.scc_count; ++i)
                prepare_linear_coefficients(linear_region_sccs_[region.scc_offset + i]);
            for (std::uint32_t i = 0; i < region.node_count; ++i) {
                const NodeId node = linear_region_nodes_[region.node_offset + i];
                linear_scratch_[node] = 1.0 - config_.dt * extra_node_decay_[node];
                region_accumulator_[node] = 0.0;
            }
            for (std::uint32_t i = 0; i < region.entry_count; ++i) {
                const auto& entry = linear_region_entries_[region.entry_offset + i];
                if (entry.column == kInvalidNode ||
                    linear_region_of_scc_[network_.node_to_scc[entry.column]] != region_id) continue;
                double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
                    ? 1.0 : linear_coefficients_[entry.reaction]);
                if (entry.kind == LinearCoefficientKind::Capped)
                    coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
                const double scaled = config_.dt * coefficient;
                if (entry.column == entry.row) linear_scratch_[entry.row] += scaled;
                else region_accumulator_[entry.row] += std::abs(scaled);
            }
            for (std::uint32_t i = 0; i < region.node_count; ++i) {
                const NodeId node = linear_region_nodes_[region.node_offset + i];
                const double row_norm = std::abs(linear_scratch_[node]) + region_accumulator_[node];
                if (row_norm > 1.0 + 1.0e-9) contractive = false;
            }
            linear_region_baseline_stationary_[region_id] =
                maximum_delta <= allowed_step_delta && maximum_boundary_delta <= 1.0e-9 &&
                contractive ? 1U : 0U;
        }
        state_ = baseline_;
    }
    ++baseline_state_id_;
    clear_cache();
}

void Engine::set_context(const Context& context) {
    std::fill(context_node_multiplier_.begin(), context_node_multiplier_.end(), 1.0F);
    std::fill(context_reaction_multiplier_.begin(), context_reaction_multiplier_.end(), 1.0F);
    for (const auto& entry : context.node_multipliers) {
        if (entry.node >= network_.nodes.size() || !std::isfinite(entry.multiplier) ||
            entry.multiplier < 0.0F || entry.multiplier > 1.0F) {
            throw std::invalid_argument("context node multiplier is invalid");
        }
        context_node_multiplier_[entry.node] = entry.multiplier;
    }
    for (const auto& entry : context.reaction_multipliers) {
        if (entry.reaction >= network_.reactions.size() || !std::isfinite(entry.multiplier) ||
            entry.multiplier < 0.0F || entry.multiplier > 1.0F) {
            throw std::invalid_argument("context reaction multiplier is invalid");
        }
        context_reaction_multiplier_[entry.reaction] = entry.multiplier;
    }
    context_id_ = context.context_id;
    baseline_.clear();
    clear_cache();
}

void Engine::clear_cache() {
    hot_cache_.clear();
}

float Engine::combine_values(
    const NodeId* nodes,
    const float* weights,
    const float* K,
    const float* n,
    std::uint32_t count,
    RegulationOp op,
    float default_K,
    float default_n) const {
    if (count == 0U) {
        return 1.0F;
    }

    auto weight_at = [weights](std::uint32_t index) {
        return weights == nullptr ? 1.0F : weights[index];
    };
    auto K_at = [K, default_K](std::uint32_t index) {
        return K == nullptr ? default_K : K[index];
    };
    auto n_at = [n, default_n](std::uint32_t index) {
        return n == nullptr ? default_n : n[index];
    };

    switch (op) {
        case RegulationOp::None:
            return clamp01(state_[nodes[0]] * weight_at(0U));
        case RegulationOp::Sum: {
            float sum = 0.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                sum += state_[nodes[index]] * weight_at(index);
            }
            return clamp01(sum);
        }
        case RegulationOp::Min: {
            float value = 1.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                value = std::min(value, clamp01(state_[nodes[index]] * weight_at(index)));
            }
            return value;
        }
        case RegulationOp::Max: {
            float value = 0.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                value = std::max(value, clamp01(state_[nodes[index]] * weight_at(index)));
            }
            return value;
        }
        case RegulationOp::ContinuousOr: {
            float inverse_product = 1.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                inverse_product *=
                    1.0F - clamp01(state_[nodes[index]] * weight_at(index));
            }
            return clamp01(1.0F - inverse_product);
        }
        case RegulationOp::HillPositive: {
            float product = 1.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                const float response = hill_positive(state_[nodes[index]], K_at(index), n_at(index));
                product *= std::pow(clamp01(response), std::max(weight_at(index), 0.0F));
            }
            return clamp01(product);
        }
        case RegulationOp::HillNegative: {
            float product = 1.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                const float response = hill_negative(state_[nodes[index]], K_at(index), n_at(index));
                product *= std::pow(clamp01(response), std::max(weight_at(index), 0.0F));
            }
            return clamp01(product);
        }
        case RegulationOp::Threshold: {
            for (std::uint32_t index = 0; index < count; ++index) {
                if (state_[nodes[index]] * weight_at(index) < K_at(index)) {
                    return 0.0F;
                }
            }
            return 1.0F;
        }
        case RegulationOp::Xor: {
            bool parity = false;
            for (std::uint32_t index = 0; index < count; ++index) {
                parity = parity != (state_[nodes[index]] * weight_at(index) >= K_at(index));
            }
            return parity ? 1.0F : 0.0F;
        }
        case RegulationOp::Product:
        case RegulationOp::ContinuousAnd: {
            float product = 1.0F;
            for (std::uint32_t index = 0; index < count; ++index) {
                const float exponent = std::max(weight_at(index), 0.0F);
                product *= std::pow(clamp01(state_[nodes[index]]), exponent);
            }
            return clamp01(product);
        }
    }
    return 0.0F;
}

float Engine::reaction_signal(const Reaction& reaction) const {
    float signal = 1.0F;
    const NodeId* inputs = reaction.input_count == 0U
                               ? nullptr
                               : network_.reaction_inputs.data() + reaction.input_offset;
    const float* input_weights = reaction.input_count == 0U
                                     ? nullptr
                                     : network_.input_stoichiometry.data() + reaction.input_offset;
    if (reaction.input_count > 0U && reaction.input_op == RegulationOp::ContinuousOr) {
        float inverse_product = 1.0F;
        for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
            inverse_product *= 1.0F -
                               clamp01(state_[inputs[index]] * input_weights[index]);
        }
        signal = clamp01(1.0F - inverse_product);
    } else if (reaction.input_count > 0U &&
               (reaction.input_op == RegulationOp::Product ||
                reaction.input_op == RegulationOp::ContinuousAnd)) {
        float product = 1.0F;
        for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
            product *= std::pow(clamp01(state_[inputs[index]]),
                                std::max(input_weights[index], 0.0F));
        }
        signal = clamp01(product);
    } else if (reaction.input_count > 0U) {
        signal = combine_values(inputs, input_weights, nullptr, nullptr, reaction.input_count,
                                reaction.input_op, reaction.Km, reaction.hill_n);
    }

    if (reaction.catalyst_count > 0U) {
        const NodeId* catalyst_nodes = network_.catalysts.data() + reaction.catalyst_offset;
        signal *= combine_values(catalyst_nodes, nullptr, nullptr, nullptr, reaction.catalyst_count,
                                 RegulationOp::Product, reaction.Km, reaction.hill_n);
    }

    if (reaction.pos_reg_count > 0U) {
        const NodeId* positive_nodes =
            network_.positive_regulators.data() + reaction.pos_reg_offset;
        const float* positive_weights = network_.positive_weights.data() + reaction.pos_reg_offset;
        const float* positive_K = network_.positive_K.data() + reaction.pos_reg_offset;
        const float* positive_n = network_.positive_n.data() + reaction.pos_reg_offset;
        signal *= combine_values(positive_nodes, positive_weights, positive_K, positive_n,
                                 reaction.pos_reg_count, reaction.pos_reg_op, reaction.Km,
                                 reaction.hill_n);
    }

    if (reaction.catalyst_count == 0U && reaction.pos_reg_count == 0U &&
        reaction.neg_reg_count == 0U && reaction.context_weight == 1.0F &&
        reaction.confidence == 1.0F && reaction_multiplier_[reaction.id] == 1.0F &&
        context_reaction_multiplier_[reaction.id] == 1.0F) {
        return signal;
    }

    if (reaction.neg_reg_count > 0U) {
        const NodeId* negative_nodes =
            network_.negative_regulators.data() + reaction.neg_reg_offset;
        const float* negative_weights = network_.negative_weights.data() + reaction.neg_reg_offset;
        const float* negative_K = network_.negative_K.data() + reaction.neg_reg_offset;
        const float* negative_n = network_.negative_n.data() + reaction.neg_reg_offset;
        float negative_factor = combine_values(
            negative_nodes, negative_weights, negative_K, negative_n, reaction.neg_reg_count,
            reaction.neg_reg_op, reaction.Km, reaction.hill_n);
        if (reaction.neg_reg_op != RegulationOp::HillNegative) {
            negative_factor = 1.0F - clamp01(negative_factor);
        }
        signal *= negative_factor;
    }

    signal *= reaction.context_weight;
    signal *= context_reaction_multiplier_[reaction.id];
    if (config_.apply_confidence) {
        signal *= reaction.confidence;
    }
    signal *= reaction_multiplier_[reaction.id];
    return clamp01(signal);
}

float Engine::reaction_flux(const Reaction& reaction) const {
    return reaction.k * reaction_signal(reaction);
}

float Engine::evaluate_scc(SccId scc, ExecutionStats& stats) {
    const SccBlock& block = network_.sccs[scc];
    std::fill_n(scratch_values_.begin(), block.node_count, 0.0F);
    std::fill_n(scratch_flags_.begin(), block.node_count, 0U);

    auto local_index = [this, scc](NodeId node) -> std::optional<std::uint32_t> {
        if (network_.node_to_scc[node] != scc) {
            return std::nullopt;
        }
        return network_.node_position_in_scc[node];
    };

    for (std::uint32_t reaction_offset = 0; reaction_offset < block.reaction_count;
         ++reaction_offset) {
        const Reaction& reaction =
            network_.reactions[network_.scc_reactions[block.reaction_offset + reaction_offset]];
        record_reaction_evaluation(scc, stats);

        if (config_.mode == SimulationMode::Boolean) {
            const float signal = reaction_signal(reaction) >= 0.5F ? 1.0F : 0.0F;
            auto assign_boolean = [&](NodeId node, float value) {
                const auto local = local_index(node);
                if (!local.has_value()) {
                    return;
                }
                if (scratch_flags_[*local] == 0U) {
                    scratch_values_[*local] = value;
                    scratch_flags_[*local] = 1U;
                } else {
                    scratch_values_[*local] = std::max(scratch_values_[*local], value);
                }
            };

            if (changes_outputs(reaction.primitive)) {
                for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                    assign_boolean(network_.reaction_outputs[reaction.output_offset + index], signal);
                }
            }
            if (consumes_inputs(reaction.primitive) && signal > 0.5F) {
                for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                    assign_boolean(network_.reaction_inputs[reaction.input_offset + index], 0.0F);
                }
            }
            continue;
        }

        if ((reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) &&
            reaction.output_count == 1U) {
            const float target = reaction_signal(reaction);
            const NodeId node = network_.reaction_outputs[reaction.output_offset];
            if (network_.node_to_scc[node] == scc) {
                const std::uint32_t local = network_.node_position_in_scc[node];
                const float stoich = network_.output_stoichiometry[reaction.output_offset];
                scratch_values_[local] += reaction.k * stoich * (target - state_[node]);
            }
            continue;
        }

        float flux = reaction_flux(reaction);
        if (consumes_inputs(reaction.primitive)) {
            for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + index];
                const float stoich = network_.input_stoichiometry[reaction.input_offset + index];
                flux = std::min(flux, state_[input] / (config_.dt * stoich));
            }
        }

        if (consumes_inputs(reaction.primitive)) {
            for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + index];
                const auto local = local_index(input);
                if (local.has_value()) {
                    const float stoich = network_.input_stoichiometry[reaction.input_offset + index];
                    scratch_values_[*local] -= stoich * flux;
                }
            }
        }
        if (changes_outputs(reaction.primitive)) {
            for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                const NodeId output = network_.reaction_outputs[reaction.output_offset + index];
                const auto local = local_index(output);
                if (local.has_value()) {
                    const float stoich = network_.output_stoichiometry[reaction.output_offset + index];
                    scratch_values_[*local] += stoich * flux;
                }
            }
        }
    }

    float maximum_delta = 0.0F;
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        float next_value = state_[node];
        if (config_.mode == SimulationMode::Boolean) {
            if (scratch_flags_[local] != 0U) {
                next_value = scratch_values_[local];
            }
            if (extra_node_decay_[node] > 0.0F) {
                next_value = 0.0F;
            }
        } else {
            const float derivative = scratch_values_[local] - extra_node_decay_[node] * state_[node];
            next_value += config_.dt * derivative;
            if (config_.clamp_normalized_state) {
                next_value = clamp01(next_value);
            }
        }
        if (clamp_active_[node] != 0U) {
            next_value = clamp_value_[node];
        }
        const float delta = std::abs(next_value - state_[node]);
        maximum_delta = std::max(maximum_delta, delta);
        if (delta > 0.0F || !std::isfinite(next_value)) {
            pending_.push_back(PendingValue{node, next_value, scc, delta});
        }
    }
    ++stats.scc_evaluations;
    ++stats.numerical_steps;
    record_scc_visit(scc, stats);
    record_scc_class(scc, false, stats);
    return maximum_delta;
}

float Engine::evaluate_scc_dirty(SccId scc, ExecutionStats& stats) {
    const SccBlock& block = network_.sccs[scc];
    const SccExecutionPlan& plan = network_.execution_plans[scc];
    std::fill_n(scratch_values_.begin(), block.node_count, 0.0F);
    std::fill_n(scratch_flags_.begin(), block.node_count, 0U);
    dirty_counts_[scc] = 0U;
    dirty_count_generation_[scc] = current_dirty_generation_;

    // Dirty scheduling is safe at the output-node level, not the reaction level:
    // every Modulate term targeting one node contributes to the same derivative
    // (including its -k*x relaxation term). If any term for a node is dirty, sum
    // all terms for that output before committing the next state. Rare primitive
    // mixtures retain the full SCC path until they have equivalent output plans.
    const bool full_scc_evaluation = config_.mode == SimulationMode::Boolean ||
        plan.transfer_count != 0U || plan.produce_count != 0U ||
        plan.destroy_count != 0U || plan.generic_count != 0U;

    auto local_index = [this, scc](NodeId node) -> std::optional<std::uint32_t> {
        if (network_.node_to_scc[node] != scc) return std::nullopt;
        return network_.node_position_in_scc[node];
    };
    auto dirty = [this, full_scc_evaluation](std::uint32_t flat_slot) {
        return full_scc_evaluation ||
               dirty_reaction_generation_[flat_slot] == current_dirty_generation_;
    };
    auto record_reaction = [this, scc, &stats]() { record_reaction_evaluation(scc, stats); };

    if (config_.mode == SimulationMode::Boolean) {
        for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
            const std::uint32_t flat_slot = block.reaction_offset + slot;
            if (!dirty(flat_slot)) continue;
            const Reaction& reaction = network_.reactions[network_.scc_reactions[flat_slot]];
            record_reaction();
            const float signal = reaction_signal(reaction) >= 0.5F ? 1.0F : 0.0F;
            auto assign_boolean = [&](NodeId node, float value) {
                const auto local = local_index(node);
                if (!local.has_value()) return;
                if (scratch_flags_[*local] == 0U) {
                    scratch_values_[*local] = value;
                    scratch_flags_[*local] = 1U;
                } else {
                    scratch_values_[*local] = std::max(scratch_values_[*local], value);
                }
            };
            if (changes_outputs(reaction.primitive)) {
                for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                    assign_boolean(network_.reaction_outputs[reaction.output_offset + index], signal);
                }
            }
            if (consumes_inputs(reaction.primitive) && signal > 0.5F) {
                for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
                    assign_boolean(network_.reaction_inputs[reaction.input_offset + index], 0.0F);
                }
            }
        }
    } else {
        if (!full_scc_evaluation) {
            for (std::uint32_t index = 0; index < plan.modulate_count; ++index) {
                const std::uint32_t plan_index = plan.modulate_offset + index;
                const std::uint32_t flat_slot = network_.plan_modulate_slots[plan_index];
                if (dirty(flat_slot)) {
                    scratch_flags_[network_.plan_modulate_output_locals[plan_index]] = 1U;
                }
            }
        }
        for (std::uint32_t index = 0; index < plan.modulate_count; ++index) {
            const std::uint32_t plan_index = plan.modulate_offset + index;
            const std::uint32_t flat_slot = network_.plan_modulate_slots[plan_index];
            const std::uint32_t local = network_.plan_modulate_output_locals[plan_index];
            if (!dirty(flat_slot) && scratch_flags_[local] == 0U) continue;
            const Reaction& reaction =
                network_.reactions[network_.plan_modulate_reactions[plan_index]];
            record_reaction();
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            const float target = reaction_signal(reaction);
            scratch_values_[local] += reaction.k *
                                      network_.plan_modulate_output_stoich[plan_index] *
                                      (target - state_[node]);
        }

        for (std::uint32_t index = 0; index < plan.transfer_count; ++index) {
            const std::uint32_t plan_index = plan.transfer_offset + index;
            const std::uint32_t flat_slot = network_.plan_transfer_slots[plan_index];
            if (!dirty(flat_slot)) continue;
            const Reaction& reaction = network_.reactions[network_.plan_transfer_reactions[plan_index]];
            record_reaction();
            float flux = reaction_flux(reaction);
            for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                const float stoich = network_.input_stoichiometry[reaction.input_offset + input_index];
                flux = std::min(flux, state_[input] / (config_.dt * stoich));
            }
            for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                const auto local = local_index(input);
                if (local.has_value()) {
                    scratch_values_[*local] -=
                        network_.input_stoichiometry[reaction.input_offset + input_index] * flux;
                }
            }
            for (std::uint32_t output_index = 0; output_index < reaction.output_count; ++output_index) {
                const NodeId output = network_.reaction_outputs[reaction.output_offset + output_index];
                const auto local = local_index(output);
                if (local.has_value()) {
                    scratch_values_[*local] +=
                        network_.output_stoichiometry[reaction.output_offset + output_index] * flux;
                }
            }
        }

        for (std::uint32_t index = 0; index < plan.produce_count; ++index) {
            const std::uint32_t plan_index = plan.produce_offset + index;
            const std::uint32_t flat_slot = network_.plan_produce_slots[plan_index];
            if (!dirty(flat_slot)) continue;
            const Reaction& reaction = network_.reactions[network_.plan_produce_reactions[plan_index]];
            record_reaction();
            const float flux = reaction_flux(reaction);
            for (std::uint32_t output_index = 0; output_index < reaction.output_count; ++output_index) {
                const NodeId output = network_.reaction_outputs[reaction.output_offset + output_index];
                const auto local = local_index(output);
                if (local.has_value()) {
                    scratch_values_[*local] +=
                        network_.output_stoichiometry[reaction.output_offset + output_index] * flux;
                }
            }
        }

        for (std::uint32_t index = 0; index < plan.destroy_count; ++index) {
            const std::uint32_t plan_index = plan.destroy_offset + index;
            const std::uint32_t flat_slot = network_.plan_destroy_slots[plan_index];
            if (!dirty(flat_slot)) continue;
            const Reaction& reaction = network_.reactions[network_.plan_destroy_reactions[plan_index]];
            record_reaction();
            float flux = reaction_flux(reaction);
            for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                const float stoich = network_.input_stoichiometry[reaction.input_offset + input_index];
                flux = std::min(flux, state_[input] / (config_.dt * stoich));
            }
            for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                const auto local = local_index(input);
                if (local.has_value()) {
                    scratch_values_[*local] -=
                        network_.input_stoichiometry[reaction.input_offset + input_index] * flux;
                }
            }
        }

        // Rare/complex forms retain the original primitive dispatch, isolated from common plan blocks.
        for (std::uint32_t index = 0; index < plan.generic_count; ++index) {
            const std::uint32_t plan_index = plan.generic_offset + index;
            const std::uint32_t flat_slot = network_.plan_generic_slots[plan_index];
            if (!dirty(flat_slot)) continue;
            const Reaction& reaction = network_.reactions[network_.plan_generic_reactions[plan_index]];
            record_reaction();
            if ((reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) &&
                reaction.output_count == 1U) {
                const NodeId node = network_.reaction_outputs[reaction.output_offset];
                const auto local = local_index(node);
                if (local.has_value()) {
                    scratch_values_[*local] += reaction.k *
                        network_.output_stoichiometry[reaction.output_offset] *
                        (reaction_signal(reaction) - state_[node]);
                }
                continue;
            }
            float flux = reaction_flux(reaction);
            if (consumes_inputs(reaction.primitive)) {
                for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                    const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                    flux = std::min(flux, state_[input] /
                                           (config_.dt * network_.input_stoichiometry[
                                               reaction.input_offset + input_index]));
                }
                for (std::uint32_t input_index = 0; input_index < reaction.input_count; ++input_index) {
                    const NodeId input = network_.reaction_inputs[reaction.input_offset + input_index];
                    const auto local = local_index(input);
                    if (local.has_value()) {
                        scratch_values_[*local] -=
                            network_.input_stoichiometry[reaction.input_offset + input_index] * flux;
                    }
                }
            }
            if (changes_outputs(reaction.primitive)) {
                for (std::uint32_t output_index = 0; output_index < reaction.output_count; ++output_index) {
                    const NodeId output = network_.reaction_outputs[reaction.output_offset + output_index];
                    const auto local = local_index(output);
                    if (local.has_value()) {
                        scratch_values_[*local] +=
                            network_.output_stoichiometry[reaction.output_offset + output_index] * flux;
                    }
                }
            }
        }
    }

    float maximum_delta = 0.0F;
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        float next_value = state_[node];
        if (config_.mode == SimulationMode::Boolean) {
            if (scratch_flags_[local] != 0U) next_value = scratch_values_[local];
            if (extra_node_decay_[node] > 0.0F) next_value = 0.0F;
        } else {
            const float derivative = scratch_values_[local] - extra_node_decay_[node] * state_[node];
            next_value += config_.dt * derivative;
            if (config_.clamp_normalized_state) next_value = clamp01(next_value);
        }
        if (clamp_active_[node] != 0U) next_value = clamp_value_[node];
        const float delta = std::abs(next_value - state_[node]);
        maximum_delta = std::max(maximum_delta, delta);
        if (delta > 0.0F || !std::isfinite(next_value)) {
            pending_.push_back(PendingValue{node, next_value, scc, delta});
        }
    }
    ++stats.scc_evaluations;
    ++stats.numerical_steps;
    record_scc_visit(scc, stats);
    record_scc_class(scc, false, stats);
    return maximum_delta;
}

bool Engine::scalar_temporal_coefficients(SccId scc, double& p, double& q) const {
    const auto& block = network_.sccs[scc];
    if (config_.mode != SimulationMode::Continuous || block.node_count != 1U ||
        block.downstream_count != 0U) return false;
    const NodeId target = network_.scc_nodes[block.node_offset];
    if (clamp_active_[target] != 0U) {
        p = 0.0; q = clamp_value_[target]; return true;
    }
    double a = -extra_node_decay_[target];
    double b = 0.0;
    for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
        const auto& reaction = network_.reactions[network_.scc_reactions[block.reaction_offset + slot]];
        if (reaction.catalyst_count || reaction.pos_reg_count || reaction.neg_reg_count) return false;
        double constant = reaction.input_count == 0U ? 1.0 : 0.0;
        double coefficient = 0.0;
        for (std::uint32_t index = 0; index < reaction.input_count; ++index) {
            const NodeId input = network_.reaction_inputs[reaction.input_offset + index];
            const double weight = network_.input_stoichiometry[reaction.input_offset + index];
            if (input == target) {
                coefficient += weight;
            } else {
                const auto& upstream = network_.temporal_plans[network_.node_to_scc[input]];
                if (upstream.classification != SccTemporalClass::Static ||
                    (extra_node_decay_[input] != 0.0F && clamp_active_[input] == 0U)) return false;
                constant += weight * state_[input];
            }
        }
        const double scale = static_cast<double>(reaction.context_weight) *
            context_reaction_multiplier_[reaction.id] * reaction_multiplier_[reaction.id] *
            (config_.apply_confidence ? reaction.confidence : 1.0F);
        constant *= scale; coefficient *= scale;
        if (coefficient == 0.0) constant = std::clamp(constant, 0.0, 1.0);
        else if (constant < 0.0 || coefficient < 0.0 || constant + coefficient > 1.0) return false;
        if (reaction.primitive == Primitive::Modulate || reaction.primitive == Primitive::Delay) {
            if (reaction.output_count != 1U || network_.reaction_outputs[reaction.output_offset] != target)
                return false;
            const double gain = reaction.k * network_.output_stoichiometry[reaction.output_offset];
            // Modulate relaxes toward reaction_signal: inhibition/context scales the
            // signal term, not the intrinsic -k*x relaxation term.
            a -= gain;
            a += gain * coefficient;
            b += gain * constant;
        } else if (reaction.primitive == Primitive::Produce) {
            for (std::uint32_t index = 0; index < reaction.output_count; ++index) {
                if (network_.reaction_outputs[reaction.output_offset + index] != target) continue;
                const double gain = reaction.k * network_.output_stoichiometry[reaction.output_offset + index];
                a += gain * coefficient; b += gain * constant;
            }
        } else if (reaction.primitive == Primitive::Destroy && reaction.input_count == 1U &&
                   network_.reaction_inputs[reaction.input_offset] == target && constant == 0.0) {
            const double stoich = network_.input_stoichiometry[reaction.input_offset];
            a -= stoich * std::min(reaction.k * coefficient, 1.0 / (config_.dt * stoich));
        } else return false;
    }
    p = 1.0 + config_.dt * a; q = config_.dt * b;
    return std::isfinite(p) && std::isfinite(q) && p >= 0.0 && q >= 0.0;
}

bool Engine::can_evaluate_scc_temporal(SccId scc,
                                       std::uint32_t remaining_steps,
                                       std::uint32_t max_steps) const {
    const SccTemporalPlan& plan = network_.temporal_plans[scc];
    const SccBlock& block = network_.sccs[scc];
    if (plan.classification == SccTemporalClass::Static) {
        // A decaying source must retain its trajectory for downstream forcing.
        // Boolean evaluation also ignores k, so k==0 does not imply static there.
        if (config_.mode == SimulationMode::Boolean && block.reaction_count != 0U) return false;
        if (config_.mode == SimulationMode::Continuous) {
            for (std::uint32_t local = 0; local < block.node_count; ++local) {
                const NodeId node = network_.scc_nodes[block.node_offset + local];
                if (clamp_active_[node] == 0U && extra_node_decay_[node] > 0.0F &&
                    (block.downstream_count != 0U ||
                     1.0F - config_.dt * extra_node_decay_[node] < 0.0F)) {
                    return false;
                }
            }
        }
        return true;
    }
    if (plan.classification != SccTemporalClass::LinearHomogeneous &&
        plan.classification != SccTemporalClass::LinearAffine) {
        return false;
    }
    if (config_.mode != SimulationMode::Continuous || config_.dt != network_.temporal_compile_dt ||
        remaining_steps != max_steps || max_steps != network_.temporal_compile_steps ||
        network_.linear_plans.empty() ||
        network_.linear_plans[scc].backend == TemporalBackend::NumericalNonlinear) return false;
    double scalar_p = 0.0, scalar_q = 0.0;
    if (scalar_temporal_coefficients(scc, scalar_p, scalar_q)) return true;
    const auto& execution = network_.linear_plans[scc];
    if (execution.boundary_count != 0U) return false;
    for (std::uint32_t index = 0; index < block.downstream_count; ++index) {
        const SccId d = network_.downstream_sccs[block.downstream_offset + index];
        if (network_.temporal_plans[d].classification != SccTemporalClass::Static) return false;
    }
    // Small/medium SCCs use dense powers; larger sparse SCCs use matrix-free
    // uniformization/power action, including runtime clamps and rate modifiers.
    return true;
}

void Engine::prepare_linear_coefficients(SccId scc) {
    const auto& block = network_.sccs[scc];
    for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
        const auto id = network_.scc_reactions[block.reaction_offset + slot];
        const auto& r = network_.reactions[id];
        linear_coefficients_[id] = static_cast<double>(r.context_weight) *
            (config_.apply_confidence ? r.confidence : 1.0F) *
            context_reaction_multiplier_[id] * reaction_multiplier_[id];
    }
}

void Engine::record_linear_execution(SccId scc, ExecutionStats& stats) {
    record_scc_visit(scc, stats); record_scc_class(scc, true, stats);
}

bool Engine::evaluate_dense_linear(SccId scc, std::uint32_t horizon, ExecutionStats& stats) {
    const auto& block = network_.sccs[scc];
    const auto& plan = network_.linear_plans[scc];
    const auto& temporal = network_.temporal_plans[scc];
    if (temporal.propagator != 0U && horizon == network_.temporal_compile_steps &&
        (temporal.bounded != 0U || !config_.clamp_normalized_state) && std::all_of(
            network_.scc_reactions.begin() + block.reaction_offset,
            network_.scc_reactions.begin() + block.reaction_offset + block.reaction_count,
            [&](ReactionId id) {
                return reaction_multiplier_[id] == 1.0F &&
                       context_reaction_multiplier_[id] == 1.0F &&
                       network_.reactions[id].confidence == 1.0F;
            })) {
        bool unmodified = true;
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            unmodified = unmodified && clamp_active_[node] == 0U && extra_node_decay_[node] == 0.0F;
        }
        if (unmodified) {
            const auto n = block.node_count;
            for (std::uint32_t row = 0; row < n; ++row) {
                const NodeId node = network_.scc_nodes[block.node_offset + row];
                double value = network_.temporal_affine[temporal.affine_offset + row];
                for (std::uint32_t column = 0; column < n; ++column) {
                    value += network_.temporal_propagators[
                        temporal.matrix_offset + static_cast<std::size_t>(row) * n + column] *
                        state_[network_.scc_nodes[block.node_offset + column]];
                }
                const float next = config_.clamp_normalized_state
                    ? clamp01(static_cast<float>(value)) : static_cast<float>(value);
                const float delta = std::abs(next - state_[node]);
                if (delta > 0.0F || !std::isfinite(next)) pending_.push_back({node, next, scc, delta});
            }
            ++stats.precomputed_propagator_executions;
            stats.propagator_memory_bytes += static_cast<std::uint64_t>(n) * n * sizeof(float);
            return true;
        }
    }
    std::uint32_t free_count = 0;
    for (std::uint32_t i = 0; i < block.node_count; ++i) {
        const NodeId node = network_.scc_nodes[block.node_offset + i];
        dense_free_index_[i] = clamp_active_[node] ? std::numeric_limits<std::uint32_t>::max() : free_count++;
    }
    const std::uint32_t dim = free_count + 1U;
    if (dim > 513U) return false;
    const std::size_t cells = static_cast<std::size_t>(dim) * dim;
    prepare_linear_coefficients(scc);
    std::uint64_t signature = 1469598103934665603ULL;
    const auto mix = [&signature](std::uint64_t value) {
        signature ^= value;
        signature *= 1099511628211ULL;
    };
    mix(scc); mix(horizon); mix(std::bit_cast<std::uint32_t>(config_.dt));
    mix(config_.clamp_normalized_state ? 1U : 0U);
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        mix(clamp_active_[node]);
        mix(std::bit_cast<std::uint32_t>(clamp_value_[node]));
        mix(std::bit_cast<std::uint32_t>(extra_node_decay_[node]));
    }
    for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
        const ReactionId id = network_.scc_reactions[block.reaction_offset + slot];
        mix(std::bit_cast<std::uint64_t>(linear_coefficients_[id]));
    }
    const auto apply_power = [&](const std::vector<double>& matrix) {
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            const auto row = dense_free_index_[local];
            if (row == std::numeric_limits<std::uint32_t>::max()) continue;
            double value = matrix[static_cast<std::size_t>(row) * dim + free_count];
            for (std::uint32_t column_local = 0; column_local < block.node_count; ++column_local) {
                const auto column = dense_free_index_[column_local];
                if (column == std::numeric_limits<std::uint32_t>::max()) continue;
                const NodeId column_node = network_.scc_nodes[block.node_offset + column_local];
                value += matrix[static_cast<std::size_t>(row) * dim + column] * state_[column_node];
            }
            const float next = config_.clamp_normalized_state
                ? clamp01(static_cast<float>(value)) : static_cast<float>(value);
            const float delta = std::abs(next - state_[node]);
            if (delta > 0.0F || !std::isfinite(next)) pending_.push_back({node, next, scc, delta});
        }
    };
    for (auto& slot : propagator_cache_) {
        if (slot.valid && slot.scc == scc && slot.dimension == dim && slot.signature == signature) {
            apply_power(slot.matrix);
            ++stats.runtime_propagator_cache_hits;
            stats.propagator_memory_bytes += cells * sizeof(double);
            return true;
        }
    }
    std::fill_n(dense_base_.begin(), cells, 0.0);
    for (std::uint32_t i = 0; i < free_count; ++i) dense_base_[static_cast<std::size_t>(i) * dim + i] = 1.0;
    dense_base_[static_cast<std::size_t>(free_count) * dim + free_count] = 1.0;
    for (std::uint32_t i = 0; i < plan.entry_count; ++i) {
        const auto& e = network_.linear_entries[plan.entry_offset + i];
        const auto row_local = network_.node_position_in_scc[e.row];
        const auto row = dense_free_index_[row_local];
        if (row == std::numeric_limits<std::uint32_t>::max()) continue;
        double coeff = e.gain * (e.kind == LinearCoefficientKind::Fixed ? 1.0 :
                                 linear_coefficients_[e.reaction]);
        if (e.kind == LinearCoefficientKind::Capped) coeff = std::min(coeff, static_cast<double>(e.stoichiometry));
        const double scaled = config_.dt * coeff;
        if (e.column == kInvalidNode) dense_base_[static_cast<std::size_t>(row) * dim + free_count] += scaled;
        else {
            const auto col_local = network_.node_position_in_scc[e.column];
            const auto col = dense_free_index_[col_local];
            if (col == std::numeric_limits<std::uint32_t>::max())
                dense_base_[static_cast<std::size_t>(row) * dim + free_count] += scaled * clamp_value_[e.column];
            else dense_base_[static_cast<std::size_t>(row) * dim + col] += scaled;
        }
    }
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        const auto row = dense_free_index_[local];
        if (row != std::numeric_limits<std::uint32_t>::max() && extra_node_decay_[node] != 0.0F)
            dense_base_[static_cast<std::size_t>(row) * dim + row] -= config_.dt * extra_node_decay_[node];
    }
    if (config_.clamp_normalized_state) for (std::uint32_t row = 0; row < free_count; ++row) {
        double sum = 0;
        for (std::uint32_t col = 0; col < dim; ++col) {
            const double x = dense_base_[static_cast<std::size_t>(row) * dim + col];
            if (x < -1.0e-12) return false;
            sum += x;
        }
        if (sum > 1.0 + 1.0e-12) return false;
    }
    std::fill_n(dense_power_.begin(), cells, 0.0);
    for (std::uint32_t i = 0; i < dim; ++i) dense_power_[static_cast<std::size_t>(i) * dim + i] = 1;
    std::uint32_t exponent = horizon;
    auto multiply = [&](const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& out) {
        std::fill_n(out.begin(), cells, 0.0);
        for (std::uint32_t i = 0; i < dim; ++i) for (std::uint32_t k = 0; k < dim; ++k) {
            const double x = a[static_cast<std::size_t>(i) * dim + k]; if (x == 0) continue;
            for (std::uint32_t j = 0; j < dim; ++j)
                out[static_cast<std::size_t>(i) * dim + j] += x * b[static_cast<std::size_t>(k) * dim + j];
        }
    };
    while (exponent) {
        if (exponent & 1U) { multiply(dense_base_, dense_power_, dense_product_); std::copy_n(dense_product_.begin(), cells, dense_power_.begin()); }
        exponent >>= 1U;
        if (exponent) { multiply(dense_base_, dense_base_, dense_product_); std::copy_n(dense_product_.begin(), cells, dense_base_.begin()); }
    }
    apply_power(dense_power_);
    auto& cache_slot = propagator_cache_[next_propagator_slot_++ % propagator_cache_.size()];
    std::copy_n(dense_power_.begin(), cells, cache_slot.matrix.begin());
    cache_slot.valid = true;
    cache_slot.signature = signature;
    cache_slot.scc = scc;
    cache_slot.dimension = dim;
    ++stats.runtime_propagator_builds;
    stats.propagator_memory_bytes += cells * sizeof(double);
    return true;
}

bool Engine::evaluate_sparse_power_action(SccId scc, std::uint32_t horizon,
                                          ExecutionStats& stats) {
    const auto& block = network_.sccs[scc];
    const auto& plan = network_.linear_plans[scc];
    if (horizon == 0U) return false;
    prepare_linear_coefficients(scc);
    std::uint64_t signature = 1469598103934665603ULL;
    const auto mix = [&signature](std::uint64_t value) {
        signature ^= value;
        signature *= 1099511628211ULL;
    };
    mix(scc); mix(horizon); mix(std::bit_cast<std::uint32_t>(config_.dt));
    mix(config_.clamp_normalized_state ? 1U : 0U);
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        mix(std::bit_cast<std::uint32_t>(state_[node]));
        mix(clamp_active_[node]);
        mix(std::bit_cast<std::uint32_t>(clamp_value_[node]));
        mix(std::bit_cast<std::uint32_t>(extra_node_decay_[node]));
    }
    for (std::uint32_t i = 0; i < plan.boundary_count; ++i) {
        const NodeId node = network_.linear_boundaries[plan.boundary_offset + i];
        mix(std::bit_cast<std::uint32_t>(state_[node]));
    }
    for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
        const ReactionId id = network_.scc_reactions[block.reaction_offset + slot];
        mix(std::bit_cast<std::uint64_t>(linear_coefficients_[id]));
    }
    for (auto& cached : power_action_cache_) {
        if (!cached.valid || cached.scc != scc || cached.horizon != horizon ||
            cached.signature != signature) continue;
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            if (clamp_active_[node]) continue;
            const float next = config_.clamp_normalized_state
                ? clamp01(static_cast<float>(cached.values[local]))
                : static_cast<float>(cached.values[local]);
            const float delta = std::abs(next - state_[node]);
            if (delta > 0.0F || !std::isfinite(next)) pending_.push_back({node, next, scc, delta});
        }
        ++stats.runtime_propagator_cache_hits;
        ++stats.sparse_power_action_executions;
        stats.power_action_terms += cached.terms;
        stats.linear_operator_applications += cached.terms > 0U ? cached.terms - 1U : 0U;
        stats.avoided_integration_steps += horizon - 1U;
        record_linear_execution(scc, stats);
        return true;
    }
    std::uint32_t free_count = 0;
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        dense_free_index_[local] = clamp_active_[node]
            ? std::numeric_limits<std::uint32_t>::max() : free_count++;
        power_action_diagonal_[node] = 1.0 - config_.dt * extra_node_decay_[node];
        power_action_offdiagonal_[node] = 0.0;
        power_action_rowsum_[node] = power_action_diagonal_[node];
    }
    bool nonnegative_entries = true;
    for (std::uint32_t index = 0; index < plan.entry_count; ++index) {
        const auto& entry = network_.linear_entries[plan.entry_offset + index];
        double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
            ? 1.0 : linear_coefficients_[entry.reaction]);
        if (entry.kind == LinearCoefficientKind::Capped)
            coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
        const double scaled = config_.dt * coefficient;
        const NodeId row = entry.row;
        if (clamp_active_[row]) continue;
        double effective = scaled;
        if (entry.column != kInvalidNode && clamp_active_[entry.column])
            effective *= clamp_value_[entry.column];
        else if (entry.column != kInvalidNode && network_.node_to_scc[entry.column] != scc)
            effective *= state_[entry.column];
        if (entry.column != row && effective < -1.0e-12) nonnegative_entries = false;
        if (entry.column == row) power_action_diagonal_[row] += scaled;
        else power_action_offdiagonal_[row] += effective;
        power_action_rowsum_[row] += effective;
    }

    // Uniformization writes the affine map as F = (1-alpha)I + alpha P.
    // When F is nonnegative and substochastic, omitted binomial tail mass is a
    // rigorous infinity-norm error bound. Otherwise alpha=1 computes F^H exactly.
    bool substochastic = nonnegative_entries;
    double alpha = 0.0;
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        if (clamp_active_[node]) continue;
        substochastic = substochastic && power_action_diagonal_[node] >= -1.0e-12 &&
            power_action_offdiagonal_[node] >= -1.0e-12 &&
            power_action_rowsum_[node] >= -1.0e-12 &&
            power_action_rowsum_[node] <= 1.0 + 1.0e-12;
        alpha = std::max(alpha, power_action_offdiagonal_[node]);
    }
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        if (!clamp_active_[node] && power_action_diagonal_[node] < 1.0 - alpha - 1.0e-12)
            substochastic = false;
    }
    if (!substochastic) alpha = 1.0;
    if (alpha <= 1.0e-15) {
        // The operator is diagonal; use a per-coordinate closed form.
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            if (clamp_active_[node]) continue;
            const double diagonal = power_action_diagonal_[node];
            const double constant = power_action_rowsum_[node] - diagonal;
            const double powered = std::pow(diagonal, static_cast<double>(horizon));
            const double geometric = std::abs(diagonal - 1.0) < 1.0e-14
                ? static_cast<double>(horizon) : (powered - 1.0) / (diagonal - 1.0);
            const float next = config_.clamp_normalized_state
                ? clamp01(static_cast<float>(powered * state_[node] + constant * geometric))
                : static_cast<float>(powered * state_[node] + constant * geometric);
            power_action_result_[node] = next;
            const float delta = std::abs(next - state_[node]);
            if (delta > 0.0F || !std::isfinite(next)) pending_.push_back({node, next, scc, delta});
        }
        auto& cached = power_action_cache_[next_power_action_slot_++ % power_action_cache_.size()];
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            cached.values[local] = clamp_active_[node] ? state_[node] : power_action_result_[node];
        }
        cached.signature = signature;
        cached.scc = scc;
        cached.horizon = horizon;
        cached.terms = 1U;
        cached.valid = true;
        ++stats.sparse_power_action_executions;
        ++stats.sparse_linear_executions;
        ++stats.linear_operator_applications;
        stats.power_action_terms += 1U;
        stats.avoided_integration_steps += horizon - 1U;
        record_linear_execution(scc, stats);
        return true;
    }

    const std::size_t constant_index = network_.nodes.size();
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        power_action_a_[node] = clamp_active_[node] ? 0.0 : state_[node];
        power_action_result_[node] = 0.0;
    }
    power_action_a_[constant_index] = 1.0;
    power_action_result_[constant_index] = 0.0;

    const auto apply_p = [&]() {
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            if (clamp_active_[node]) continue;
            power_action_b_[node] = power_action_a_[node];
        }
        for (std::uint32_t index = 0; index < plan.entry_count; ++index) {
            const auto& entry = network_.linear_entries[plan.entry_offset + index];
            if (clamp_active_[entry.row]) continue;
            double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
                ? 1.0 : linear_coefficients_[entry.reaction]);
            if (entry.kind == LinearCoefficientKind::Capped)
                coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
            const double value = entry.column == kInvalidNode ? 1.0 :
                (clamp_active_[entry.column] ? clamp_value_[entry.column] :
                 (network_.node_to_scc[entry.column] != scc ? state_[entry.column] :
                                                              power_action_a_[entry.column]));
            power_action_b_[entry.row] += config_.dt * coefficient * value;
        }
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            if (clamp_active_[node]) continue;
            power_action_b_[node] -= config_.dt * extra_node_decay_[node] * power_action_a_[node];
            power_action_b_[node] = (power_action_b_[node] - (1.0 - alpha) * power_action_a_[node]) / alpha;
        }
        power_action_b_[constant_index] = 1.0;
        power_action_a_.swap(power_action_b_);
    };

    long double weight = std::pow(1.0L - static_cast<long double>(alpha),
                                  static_cast<long double>(horizon));
    long double cumulative = 0.0L;
    std::uint32_t terms = 0U;
    const std::uint32_t max_term = substochastic && config_.clamp_normalized_state
        ? horizon : horizon;
    for (std::uint32_t power = 0; power <= max_term; ++power) {
        if (power > 0U) {
            if (alpha >= 1.0 - 1.0e-15)
                weight = power == horizon ? 1.0L : 0.0L;
            else weight *= static_cast<long double>(horizon - power + 1U) /
                static_cast<long double>(power) * static_cast<long double>(alpha) /
                (1.0L - static_cast<long double>(alpha));
        }
        if (weight != 0.0L) {
            for (std::uint32_t local = 0; local < block.node_count; ++local) {
                const NodeId node = network_.scc_nodes[block.node_offset + local];
                if (!clamp_active_[node]) power_action_result_[node] +=
                    static_cast<double>(weight) * power_action_a_[node];
            }
            cumulative += weight;
        }
        ++terms;
        if (power == max_term) break;
        if (power >= horizon) break;
        if (power != 0U && substochastic && config_.clamp_normalized_state &&
            1.0L - cumulative < 2.0e-7L) break;
        apply_p();
    }

    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        if (clamp_active_[node]) continue;
        const float raw = static_cast<float>(power_action_result_[node]);
        const float next = config_.clamp_normalized_state ? clamp01(raw) : raw;
        const float delta = std::abs(next - state_[node]);
        if (delta > 0.0F || !std::isfinite(next)) pending_.push_back({node, next, scc, delta});
    }
    auto& cached = power_action_cache_[next_power_action_slot_++ % power_action_cache_.size()];
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        cached.values[local] = clamp_active_[node] ? state_[node] : power_action_result_[node];
    }
    cached.signature = signature;
    cached.scc = scc;
    cached.horizon = horizon;
    cached.terms = terms;
    cached.valid = true;
    ++stats.sparse_power_action_executions;
    ++stats.sparse_linear_executions;
    stats.linear_operator_applications += terms > 0U ? terms - 1U : 0U;
    stats.power_action_terms += terms;
    stats.avoided_integration_steps += horizon - 1U;
    record_linear_execution(scc, stats);
    return true;
}

void Engine::evaluate_sparse_linear(SccId scc, ExecutionStats& stats) {
    const auto& block = network_.sccs[scc]; const auto& plan = network_.linear_plans[scc];
    prepare_linear_coefficients(scc);
    for (std::uint32_t local = 0; local < block.node_count; ++local)
        linear_scratch_[network_.scc_nodes[block.node_offset + local]] = 0;
    for (std::uint32_t i = 0; i < plan.entry_count; ++i) {
        const auto& e = network_.linear_entries[plan.entry_offset + i];
        double coeff = e.gain * (e.kind == LinearCoefficientKind::Fixed ? 1.0 :
                                 linear_coefficients_[e.reaction]);
        if (e.kind == LinearCoefficientKind::Capped) coeff = std::min(coeff, static_cast<double>(e.stoichiometry));
        const double value = e.column == kInvalidNode ? 1.0 : state_[e.column];
        linear_scratch_[e.row] += coeff * value;
    }
    if (plan.boundary_count) {
        stats.boundary_input_steps += plan.boundary_count;
        for (std::uint32_t i = 0; i < plan.boundary_count; ++i) {
            const NodeId node = network_.linear_boundaries[plan.boundary_offset + i];
            float& previous = boundary_previous_[node];
            if (boundary_previous_generation_[node] != current_boundary_generation_) {
                boundary_previous_generation_[node] = current_boundary_generation_;
                previous = std::numeric_limits<float>::quiet_NaN();
            }
            if (!std::isfinite(previous) || std::abs(previous - state_[node]) > config_.epsilon) {
                ++stats.boundary_segments; previous = state_[node];
            }
        }
        ++stats.boundary_driven_linear_executions;
        stats.trajectory_compression_ratio = static_cast<double>(stats.boundary_input_steps) /
            std::max<std::uint64_t>(1, stats.boundary_segments);
    } else ++stats.sparse_linear_executions;
    for (std::uint32_t local = 0; local < block.node_count; ++local) {
        const NodeId node = network_.scc_nodes[block.node_offset + local];
        float value = state_[node] + config_.dt * static_cast<float>(linear_scratch_[node] - extra_node_decay_[node] * state_[node]);
        if (config_.clamp_normalized_state) value = clamp01(value);
        if (clamp_active_[node]) value = clamp_value_[node];
        const float delta = std::abs(value - state_[node]);
        if (delta > 0 || !std::isfinite(value)) pending_.push_back({node, value, scc, delta});
    }
    ++stats.linear_operator_applications;
    ++stats.scc_evaluations;
    ++stats.analytic_updates;
    record_linear_execution(scc, stats);
}

bool Engine::evaluate_sparse_power_action_region(std::uint32_t region_id,
                                                 std::uint32_t horizon,
                                                 ExecutionStats& stats,
                                                 bool boundary_only,
                                                 bool homogeneous_only,
                                                 bool suppress_boundary,
                                                 std::span<const std::uint32_t> active_nodes,
                                                 std::span<const std::uint32_t> active_entries) {
    if (horizon == 0U || region_id >= linear_regions_.size()) return false;
    const auto& region = linear_regions_[region_id];
    const std::size_t constant = network_.nodes.size();
    const bool restricted = !active_nodes.empty();
    const std::uint64_t visited_node_count = restricted ? active_nodes.size() : region.node_count;
    const std::uint64_t visited_entry_count = restricted ? active_entries.size() : region.entry_count;
    const auto visit_nodes = [&](auto&& function) {
        if (restricted) {
            for (const std::uint32_t slot : active_nodes)
                function(linear_region_nodes_[slot]);
        } else {
            for (std::uint32_t i = 0; i < region.node_count; ++i)
                function(linear_region_nodes_[region.node_offset + i]);
        }
    };
    const auto visit_entries = [&](auto&& function) {
        if (restricted) {
            for (const std::uint32_t slot : active_entries)
                function(linear_region_entries_[slot]);
        } else {
            for (std::uint32_t i = 0; i < region.entry_count; ++i)
                function(linear_region_entries_[region.entry_offset + i]);
        }
    };
    const auto active_internal_node = [&](NodeId node) {
        if (node == kInvalidNode ||
            linear_region_of_scc_[network_.node_to_scc[node]] != region_id) return false;
        if (!restricted) return true;
        const std::uint32_t slot = linear_region_node_slot_[node];
        return slot != std::numeric_limits<std::uint32_t>::max() &&
            causal_cone_generation_[slot] == current_causal_cone_generation_;
    };
    if (restricted) {
        for (const std::uint32_t slot : active_nodes)
            prepare_linear_coefficients(network_.node_to_scc[linear_region_nodes_[slot]]);
    } else {
        for (std::uint32_t i = 0; i < region.scc_count; ++i)
            prepare_linear_coefficients(linear_region_sccs_[region.scc_offset + i]);
    }
    visit_nodes([&](NodeId node) {
        power_action_diagonal_[node] = 1.0 - config_.dt * extra_node_decay_[node];
        power_action_offdiagonal_[node] = 0.0;
        power_action_rowsum_[node] = power_action_diagonal_[node];
        power_action_a_[node] = clamp_active_[node] ? 0.0 : state_[node];
        power_action_result_[node] = 0.0;
    });
    bool nonnegative = true;
    double alpha = 0.0;
    visit_entries([&](const LinearEntry& entry) {
        if (clamp_active_[entry.row]) return;
        const bool internal = active_internal_node(entry.column);
        if (suppress_boundary && entry.column != kInvalidNode && !internal) return;
        if (homogeneous_only &&
            (entry.column == kInvalidNode || !internal || clamp_active_[entry.column])) return;
        if (boundary_only && (entry.column == kInvalidNode ||
            (internal && clamp_active_[entry.column]))) return;
        double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
            ? 1.0 : linear_coefficients_[entry.reaction]);
        if (entry.kind == LinearCoefficientKind::Capped)
            coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
        const double scaled = config_.dt * coefficient;
        const double factor = entry.column == kInvalidNode ? 1.0 :
            (internal && !clamp_active_[entry.column] ? 1.0 :
             (clamp_active_[entry.column] ? clamp_value_[entry.column] : state_[entry.column]));
        const double effective = internal && !clamp_active_[entry.column] ? scaled : scaled * factor;
        if (entry.column == entry.row) power_action_diagonal_[entry.row] += scaled;
        else power_action_offdiagonal_[entry.row] += effective;
        power_action_rowsum_[entry.row] += effective;
        if (entry.column != entry.row && effective < -1.0e-12) nonnegative = false;
        alpha = std::max(alpha, std::abs(effective));
    });
    bool substochastic = nonnegative;
    visit_nodes([&](NodeId node) {
        if (clamp_active_[node]) return;
        substochastic = substochastic && power_action_diagonal_[node] >= -1.0e-12 &&
            power_action_offdiagonal_[node] >= -1.0e-12 && power_action_rowsum_[node] >= -1.0e-12 &&
            power_action_rowsum_[node] <= 1.0 + 1.0e-12;
        if (power_action_diagonal_[node] < 1.0 - alpha - 1.0e-12) substochastic = false;
    });
    if (!substochastic) alpha = 1.0;
    if (alpha <= 1.0e-15) {
        visit_nodes([&](NodeId node) {
            if (clamp_active_[node]) return;
            const double diagonal = power_action_diagonal_[node];
            const double forcing = power_action_rowsum_[node] - diagonal;
            const double powered = std::pow(diagonal, static_cast<double>(horizon));
            const double geometric = std::abs(diagonal - 1.0) < 1.0e-14
                ? static_cast<double>(horizon) : (powered - 1.0) / (diagonal - 1.0);
            power_action_result_[node] = powered * state_[node] + forcing * geometric;
        });
        stats.power_action_terms += 1U;
        ++stats.sparse_power_action_executions;
        ++stats.sparse_linear_executions;
        ++stats.linear_operator_applications;
    } else {
        power_action_a_[constant] = 1.0;
        const auto apply_p = [&]() {
            ++stats.sparse_operator_passes;
            stats.region_nodes_touched += visited_node_count;
            stats.region_edges_touched += visited_entry_count;
            stats.spmv_equivalent_operations += visited_node_count + visited_entry_count;
            visit_nodes([&](NodeId node) {
                if (!clamp_active_[node]) power_action_b_[node] = power_action_a_[node];
            });
            visit_entries([&](const LinearEntry& entry) {
                if (clamp_active_[entry.row]) return;
                const bool internal = active_internal_node(entry.column);
                if (suppress_boundary && entry.column != kInvalidNode && !internal) return;
                if (homogeneous_only &&
                    (entry.column == kInvalidNode || !internal || clamp_active_[entry.column])) return;
                if (boundary_only && (entry.column == kInvalidNode ||
                    (internal && clamp_active_[entry.column]))) return;
                double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
                    ? 1.0 : linear_coefficients_[entry.reaction]);
                if (entry.kind == LinearCoefficientKind::Capped)
                    coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
                double value = 1.0;
                if (entry.column != kInvalidNode) {
                    if (clamp_active_[entry.column]) value = clamp_value_[entry.column];
                    else if (internal)
                        value = power_action_a_[entry.column];
                    else value = state_[entry.column];
                }
                power_action_b_[entry.row] += config_.dt * coefficient * value;
            });
            visit_nodes([&](NodeId node) {
                if (clamp_active_[node]) return;
                power_action_b_[node] -= config_.dt * extra_node_decay_[node] * power_action_a_[node];
                power_action_b_[node] = (power_action_b_[node] - (1.0 - alpha) * power_action_a_[node]) / alpha;
            });
            power_action_b_[constant] = 1.0;
            power_action_a_.swap(power_action_b_);
        };
        long double weight = std::pow(1.0L - static_cast<long double>(alpha),
                                      static_cast<long double>(horizon));
        long double cumulative = 0.0L;
        std::uint32_t terms = 0U;
        for (std::uint32_t power = 0U; power <= horizon; ++power) {
            if (power > 0U) {
                if (alpha >= 1.0 - 1.0e-15) weight = power == horizon ? 1.0L : 0.0L;
                else weight *= static_cast<long double>(horizon - power + 1U) /
                    static_cast<long double>(power) * static_cast<long double>(alpha) /
                    (1.0L - static_cast<long double>(alpha));
            }
            if (weight != 0.0L) {
                visit_nodes([&](NodeId node) {
                    if (!clamp_active_[node]) power_action_result_[node] +=
                        static_cast<double>(weight) * power_action_a_[node];
                });
                cumulative += weight;
            }
            ++terms;
            if (power == horizon || (power > 0U && substochastic && config_.clamp_normalized_state &&
                                     1.0L - cumulative < 2.0e-7L)) break;
            apply_p();
        }
        stats.power_action_terms += terms;
        stats.linear_operator_applications += terms > 0U ? terms - 1U : 0U;
        ++stats.sparse_power_action_executions;
        ++stats.sparse_linear_executions;
    }
    visit_nodes([&](NodeId node) {
        if (clamp_active_[node]) return;
        float next = static_cast<float>(power_action_result_[node]);
        if (config_.clamp_normalized_state) next = clamp01(next);
        const float delta = std::abs(next - state_[node]);
        if (delta > 0.0F || !std::isfinite(next))
            pending_.push_back({node, next, network_.node_to_scc[node], delta});
    });
    stats.avoided_integration_steps += horizon > 0U ? horizon - 1U : 0U;
    return true;
}

void Engine::execute_deferred_boundary_sccs(std::uint32_t horizon, ExecutionStats& stats) {
    if (horizon == 0U) return;
    const std::uint32_t trace_width = static_cast<std::uint32_t>(boundary_trace_nodes_.size());
    const float tolerance = config_.epsilon * std::min(config_.dt, 1.0F);
    for (std::uint32_t region_id = 0; region_id < linear_regions_.size(); ++region_id) {
        const auto& region = linear_regions_[region_id];
        if (!region.deferred) continue;
        const auto region_start = std::chrono::steady_clock::now();
        stats.region_nodes_total += region.node_count;
        stats.region_edges_total += region.entry_count;
        stats.fused_region_count += 1U;
        stats.boundary_forcing_aggregations += region.boundary_count > 1U ? 1U : 0U;

        active_boundary_nodes_.clear();
        for (std::uint32_t i = 0; i < region.boundary_count; ++i) {
            const NodeId node = linear_region_boundaries_[region.boundary_offset + i];
            boundary_saved_state_[boundary_trace_index_[node]] = state_[node];
            bool changed = std::abs(state_[node] - baseline_[node]) > tolerance;
            for (std::uint32_t step = 0; !changed && step < horizon; ++step) {
                const float value = boundary_trajectory_[static_cast<std::size_t>(step) * trace_width +
                                                          boundary_trace_index_[node]];
                changed = std::abs(value - baseline_[node]) > tolerance;
            }
            if (changed) active_boundary_nodes_.push_back(node);
        }

        const bool full_region_guard = region_id >= linear_region_baseline_stationary_.size() ||
            linear_region_baseline_stationary_[region_id] == 0U;
        if (++current_causal_cone_generation_ == 0U) {
            std::fill(causal_cone_generation_.begin(), causal_cone_generation_.end(), 0U);
            current_causal_cone_generation_ = 1U;
        }
        if (full_region_guard) {
            ++stats.causal_cone_full_region_guard;
            active_region_node_slots_.clear();
            active_region_entry_slots_.clear();
            active_boundary_nodes_.clear();
            for (std::uint32_t i = 0; i < region.node_count; ++i) {
                const std::uint32_t slot = region.node_offset + i;
                active_region_node_slots_.push_back(slot);
                causal_cone_generation_[slot] = current_causal_cone_generation_;
            }
            for (std::uint32_t i = 0; i < region.entry_count; ++i)
                active_region_entry_slots_.push_back(region.entry_offset + i);
            for (std::uint32_t i = 0; i < region.boundary_count; ++i)
                active_boundary_nodes_.push_back(
                    linear_region_boundaries_[region.boundary_offset + i]);
        } else {
        causal_cone_queue_.clear();
        const auto add_cone_seed = [&](std::uint32_t slot, std::uint32_t depth) {
            if (slot < region.node_offset || slot >= region.node_offset + region.node_count) return;
            if (causal_cone_generation_[slot] == current_causal_cone_generation_) {
                if (depth >= causal_cone_depth_[slot]) return;
                causal_cone_depth_[slot] = depth;
                causal_cone_queue_.push_back(slot);
                return;
            }
            causal_cone_generation_[slot] = current_causal_cone_generation_;
            causal_cone_depth_[slot] = depth;
            causal_cone_queue_.push_back(slot);
        };
        for (const NodeId source : active_boundary_nodes_) {
            for (std::uint64_t dep = network_.dependency_offsets[source];
                 dep < network_.dependency_offsets[source + 1U]; ++dep) {
                const NodeId target = network_.dependency_targets[dep];
                const std::uint32_t slot = linear_region_node_slot_[target];
                if (slot != std::numeric_limits<std::uint32_t>::max() &&
                    linear_region_of_scc_[network_.node_to_scc[target]] == region_id)
                    add_cone_seed(slot, 1U);
            }
        }
        for (const auto& change : journal_) {
            if (change.node >= linear_region_node_slot_.size()) continue;
            const std::uint32_t slot = linear_region_node_slot_[change.node];
            if (slot == std::numeric_limits<std::uint32_t>::max()) continue;
            if (linear_region_of_scc_[network_.node_to_scc[change.node]] == region_id &&
                std::abs(state_[change.node] - baseline_[change.node]) > tolerance)
                add_cone_seed(slot, 0U);
        }
        // Protocols may retain multiple coefficient-changing interventions. Scan the
        // current persistent modifiers, not only the most recently scheduled event.
        // Protein degradation and clamping can change a region row without a journaled
        // state delta; reaction perturbations can change it without changing any node.
        for (std::uint32_t i = 0; i < region.node_count; ++i) {
            const std::uint32_t slot = region.node_offset + i;
            const NodeId node = linear_region_nodes_[slot];
            if (extra_node_decay_[node] > 0.0F || clamp_active_[node] != 0U)
                add_cone_seed(slot, 0U);
        }
        for (std::uint32_t i = 0; i < region.entry_count; ++i) {
            const std::uint32_t entry_slot = region.entry_offset + i;
            const auto& entry = linear_region_entries_[entry_slot];
            if (entry.reaction >= reaction_multiplier_.size() ||
                reaction_multiplier_[entry.reaction] == 1.0F) continue;
            const std::uint32_t slot = linear_region_node_slot_[entry.row];
            if (slot != std::numeric_limits<std::uint32_t>::max() &&
                linear_region_of_scc_[network_.node_to_scc[entry.row]] == region_id)
                add_cone_seed(slot, 0U);
        }
        std::uint64_t signature = network_.network_version;
        const auto mix_signature = [&](std::uint64_t value) {
            signature ^= value + 0x9e3779b97f4a7c15ULL + (signature << 6U) + (signature >> 2U);
        };
        mix_signature(region_id);
        mix_signature(horizon);
        mix_signature(context_id_);
        mix_signature(baseline_state_id_);
        for (const std::uint32_t slot : causal_cone_queue_) {
            mix_signature(slot);
            mix_signature(causal_cone_depth_[slot]);
        }
        const bool cone_cache_hit = causal_cone_cache_valid_ &&
            causal_cone_cache_signature_ == signature &&
            causal_cone_cache_region_ == region_id && causal_cone_cache_horizon_ == horizon &&
            causal_cone_cache_context_id_ == context_id_ &&
            causal_cone_cache_baseline_id_ == baseline_state_id_;
        if (cone_cache_hit) {
            active_region_node_slots_ = cached_causal_cone_nodes_;
            active_region_entry_slots_ = cached_causal_cone_entries_;
            active_boundary_nodes_ = cached_causal_boundary_nodes_;
            for (const std::uint32_t slot : active_region_node_slots_)
                causal_cone_generation_[slot] = current_causal_cone_generation_;
            ++stats.causal_cone_cache_hits;
        } else {
            ++stats.causal_cone_cache_misses;
        }
        if (!cone_cache_hit) for (std::size_t head = 0; head < causal_cone_queue_.size(); ++head) {
            const std::uint32_t source_slot = causal_cone_queue_[head];
            const std::uint32_t next_depth = causal_cone_depth_[source_slot] + 1U;
            if (next_depth > horizon) continue;
            for (std::uint64_t edge = linear_region_dependency_offsets_[source_slot];
                 edge < linear_region_dependency_offsets_[source_slot + 1U]; ++edge) {
                const std::uint32_t target = linear_region_dependency_targets_[edge];
                if (linear_region_of_scc_[network_.node_to_scc[linear_region_nodes_[target]]] != region_id ||
                    causal_cone_generation_[target] == current_causal_cone_generation_) continue;
                causal_cone_generation_[target] = current_causal_cone_generation_;
                causal_cone_depth_[target] = next_depth;
                causal_cone_queue_.push_back(target);
            }
        }
        if (!cone_cache_hit) {
        std::sort(causal_cone_queue_.begin(), causal_cone_queue_.end());
        causal_cone_queue_.erase(std::unique(causal_cone_queue_.begin(), causal_cone_queue_.end()),
                                 causal_cone_queue_.end());
        active_region_node_slots_.assign(causal_cone_queue_.begin(), causal_cone_queue_.end());
        active_boundary_nodes_.clear();
        active_region_entry_slots_.clear();
        for (const std::uint32_t row_slot : active_region_node_slots_) {
            for (std::uint64_t position = linear_region_entry_row_offsets_[row_slot];
                 position < linear_region_entry_row_offsets_[row_slot + 1U]; ++position) {
                const std::uint32_t entry_slot = linear_region_entry_row_indices_[position];
                const auto& entry = linear_region_entries_[entry_slot];
                // Keep incoming entries from same-region nodes outside the perturbation cone.
                // They are fixed baseline forcing, not a causal path to be traversed, but
                // omitting them changes the downstream linear trajectory itself.
                active_region_entry_slots_.push_back(entry_slot);
                if (entry.column != kInvalidNode &&
                    linear_region_of_scc_[network_.node_to_scc[entry.column]] != region_id &&
                    std::find(active_boundary_nodes_.begin(), active_boundary_nodes_.end(), entry.column) ==
                        active_boundary_nodes_.end())
                    active_boundary_nodes_.push_back(entry.column);
            }
        }
            cached_causal_cone_nodes_ = active_region_node_slots_;
            cached_causal_cone_entries_ = active_region_entry_slots_;
            cached_causal_boundary_nodes_ = active_boundary_nodes_;
            causal_cone_cache_signature_ = signature;
            causal_cone_cache_region_ = region_id;
            causal_cone_cache_horizon_ = horizon;
            causal_cone_cache_context_id_ = context_id_;
            causal_cone_cache_baseline_id_ = baseline_state_id_;
            causal_cone_cache_valid_ = true;
        }
        }
        const auto active_nodes = std::span<const std::uint32_t>(active_region_node_slots_);
        const auto active_entries = std::span<const std::uint32_t>(active_region_entry_slots_);
        stats.causal_cone_nodes += active_nodes.size();
        stats.causal_cone_edges += active_entries.size();
        stats.boundary_channels += active_boundary_nodes_.size();
        if (active_nodes.empty()) {
            stats.time_region_execution_us += std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - region_start).count();
            continue;
        }
        for (const std::uint32_t slot : active_nodes) {
            const NodeId node = linear_region_nodes_[slot];
            region_initial_state_[node] = state_[node];
            region_accumulator_[node] = 0.0;
        }

        // Base response holds all upstream boundary values at their final values. The
        // following single sparse recurrence adds the trajectory deviations as one forcing.
        pending_.clear();
        const auto power_start = std::chrono::steady_clock::now();
        ++stats.power_action_calls;
        ++stats.region_passes_per_perturbation;
        if (!evaluate_sparse_power_action_region(region_id, horizon, stats, false, false, false,
                                                  active_nodes, active_entries)) {
            stats.numerical_error = true;
            return;
        }
        stats.time_power_action_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - power_start).count();
        for (const std::uint32_t slot : active_nodes) {
            const NodeId node = linear_region_nodes_[slot];
            if (!clamp_active_[node]) region_accumulator_[node] = power_action_result_[node];
        }

        const auto same_segment = [&](std::uint32_t lhs, std::uint32_t rhs) {
            if (!config_.enable_boundary_compression) return false;
            for (const NodeId node : active_boundary_nodes_) {
                const std::uint32_t trace = boundary_trace_index_[node];
                if (std::abs(boundary_trajectory_[static_cast<std::size_t>(lhs) * trace_width + trace] -
                             boundary_trajectory_[static_cast<std::size_t>(rhs) * trace_width + trace]) > tolerance)
                    return false;
            }
            return true;
        };
        const auto compression_start = std::chrono::steady_clock::now();
        std::uint32_t start = 0U;
        std::uint64_t region_segments = 1U;
        for (std::uint32_t step = 0U; step < horizon; ++step) {
            if (step > 0U && !same_segment(start, step)) {
                start = step;
                ++region_segments;
            }
            boundary_segment_start_[step] = start;
        }
        stats.time_boundary_compression_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - compression_start).count();
        bool has_dynamic_forcing = false;
        const auto recording_start = std::chrono::steady_clock::now();
        for (const NodeId source : active_boundary_nodes_) {
            const std::uint32_t trace = boundary_trace_index_[source];
            for (std::uint32_t step = 0; step < horizon; ++step) {
                const float input = boundary_trajectory_[static_cast<std::size_t>(
                    boundary_segment_start_[step]) * trace_width + trace];
                if (std::abs(input - boundary_saved_state_[trace]) > tolerance) {
                    has_dynamic_forcing = true;
                    break;
                }
            }
            if (has_dynamic_forcing) break;
        }
        stats.time_boundary_recording_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - recording_start).count();
        const std::uint64_t recurrence_cost = static_cast<std::uint64_t>(horizon) *
            (active_nodes.size() + active_entries.size());
        const std::uint64_t independent_boundary_cost = recurrence_cost *
            std::max<std::size_t>(1U, active_boundary_nodes_.size());
        stats.backend_selected_by_cost_model = true;
        if (has_dynamic_forcing && recurrence_cost <= independent_boundary_cost)
            ++stats.selected_sparse_recurrence_runs;
        else
            ++stats.selected_sparse_power_action_runs;
        if (has_dynamic_forcing) {
            const auto response_start = std::chrono::steady_clock::now();
            ++stats.response_kernel_applications;
            ++stats.sparse_power_action_executions;
            ++stats.power_action_calls;
            for (const std::uint32_t slot : active_nodes)
                power_action_a_[linear_region_nodes_[slot]] = 0.0;
            for (std::uint32_t step = 0; step < horizon; ++step) {
                ++stats.sparse_operator_passes;
                ++stats.region_passes_per_perturbation;
                stats.region_nodes_touched += active_nodes.size();
                stats.region_edges_touched += active_entries.size();
                stats.spmv_equivalent_operations += active_nodes.size() + active_entries.size();
                for (const std::uint32_t slot : active_nodes) {
                    const NodeId node = linear_region_nodes_[slot];
                    power_action_b_[node] = clamp_active_[node] ? 0.0 :
                        power_action_a_[node] * (1.0 - config_.dt * extra_node_decay_[node]);
                }
                for (const std::uint32_t entry_slot : active_entries) {
                    const auto& entry = linear_region_entries_[entry_slot];
                    if (clamp_active_[entry.row]) continue;
                    double coefficient = entry.gain * (entry.kind == LinearCoefficientKind::Fixed
                        ? 1.0 : linear_coefficients_[entry.reaction]);
                    if (entry.kind == LinearCoefficientKind::Capped)
                        coefficient = std::min(coefficient, static_cast<double>(entry.stoichiometry));
                    const bool same_region = entry.column != kInvalidNode &&
                        linear_region_of_scc_[network_.node_to_scc[entry.column]] == region_id;
                    const std::uint32_t column_slot = same_region
                        ? linear_region_node_slot_[entry.column]
                        : std::numeric_limits<std::uint32_t>::max();
                    const bool active_column = column_slot != std::numeric_limits<std::uint32_t>::max() &&
                        causal_cone_generation_[column_slot] == current_causal_cone_generation_;
                    if (same_region && active_column) {
                        if (!clamp_active_[entry.column])
                            power_action_b_[entry.row] += config_.dt * coefficient * power_action_a_[entry.column];
                    } else if (!same_region && entry.column != kInvalidNode) {
                        const std::uint32_t trace = boundary_trace_index_[entry.column];
                        const std::uint32_t sample_step = boundary_segment_start_[step];
                        const float input = boundary_trajectory_[static_cast<std::size_t>(sample_step) *
                                                                  trace_width + trace];
                        const double delta = static_cast<double>(input) - boundary_saved_state_[trace];
                        power_action_b_[entry.row] += config_.dt * coefficient * delta;
                    }
                }
                power_action_a_.swap(power_action_b_);
            }
            stats.power_action_terms += horizon;
            stats.linear_operator_applications += horizon;
            for (const std::uint32_t slot : active_nodes) {
                const NodeId node = linear_region_nodes_[slot];
                if (!clamp_active_[node]) region_accumulator_[node] += power_action_a_[node];
            }
            stats.time_power_action_us += std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - response_start).count();
        }
        stats.boundary_segments += region_segments;
        stats.boundary_steps += static_cast<std::uint64_t>(active_boundary_nodes_.size()) * horizon;
        stats.boundary_input_steps += static_cast<std::uint64_t>(active_boundary_nodes_.size()) * horizon;
        if (stats.numerical_error) return;
        pending_.clear();
        active_region_sccs_.clear();
        if (++current_causal_scc_generation_ == 0U) {
            std::fill(causal_cone_scc_generation_.begin(), causal_cone_scc_generation_.end(), 0U);
            current_causal_scc_generation_ = 1U;
        }
        for (const std::uint32_t slot : active_nodes) {
            const NodeId node = linear_region_nodes_[slot];
            state_[node] = static_cast<float>(region_initial_state_[node]);
            const float value = clamp_active_[node] ? clamp_value_[node] :
                static_cast<float>(region_accumulator_[node]);
            const float delta = std::abs(value - state_[node]);
            if (delta > 0.0F || !std::isfinite(value))
                pending_.push_back({node, value, network_.node_to_scc[node], delta});
            const SccId scc = network_.node_to_scc[node];
            if (causal_cone_scc_generation_[scc] != current_causal_scc_generation_) {
                causal_cone_scc_generation_[scc] = current_causal_scc_generation_;
                active_region_sccs_.push_back(scc);
            }
        }
        static_cast<void>(commit_pending(stats));
        for (const SccId scc : active_region_sccs_) {
            ++stats.scc_evaluations;
            ++stats.analytic_updates;
            record_linear_execution(scc, stats);
        }
        if (active_region_sccs_.size() > 1U) ++stats.fused_linear_regions;
        if (stats.boundary_segments > 0U && stats.boundary_steps > 0U)
            stats.trajectory_compression_ratio = static_cast<double>(stats.boundary_steps) /
                static_cast<double>(stats.boundary_segments);
        stats.avoided_linear_integration_steps +=
            static_cast<std::uint64_t>(horizon) * active_region_sccs_.size();
        stats.time_region_execution_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - region_start).count();
        if (stats.numerical_error) return;
    }
}

float Engine::evaluate_scc_temporal(SccId scc,
                                    std::uint32_t remaining_steps,
                                    std::uint32_t max_steps,
                                    ExecutionStats& stats) {
    const SccBlock& block = network_.sccs[scc];
    const SccTemporalPlan& plan = network_.temporal_plans[scc];
    dirty_counts_[scc] = 0U;
    dirty_count_generation_[scc] = current_dirty_generation_;
    float maximum_delta = 0.0F;
    bool completed_horizon = false;

    if (plan.classification == SccTemporalClass::Static) {
        completed_horizon = true;
        ++stats.static_scc_executions;
        ++stats.analytic_updates;
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network_.scc_nodes[block.node_offset + local];
            float next_value = state_[node];
            if (config_.mode == SimulationMode::Boolean) {
                if (extra_node_decay_[node] > 0.0F) next_value = 0.0F;
            } else if (extra_node_decay_[node] > 0.0F && clamp_active_[node] == 0U) {
                const float factor = 1.0F - config_.dt * extra_node_decay_[node];
                next_value = state_[node] *
                             static_cast<float>(std::pow(static_cast<double>(factor),
                                                         static_cast<double>(remaining_steps)));
                if (config_.clamp_normalized_state) next_value = clamp01(next_value);
            }
            if (clamp_active_[node] != 0U) next_value = clamp_value_[node];
            const float delta = std::abs(next_value - state_[node]);
            maximum_delta = std::max(maximum_delta, delta);
            if (delta > 0.0F || !std::isfinite(next_value)) {
                pending_.push_back(PendingValue{node, next_value, scc, delta});
            }
        }
    } else {
        double scalar_p = 0.0, scalar_q = 0.0;
        if (scalar_temporal_coefficients(scc, scalar_p, scalar_q)) {
            const NodeId node = network_.scc_nodes[block.node_offset];
            const double power = std::pow(scalar_p, remaining_steps);
            const double sum = scalar_p == 1.0 ? remaining_steps :
                (scalar_p == 0.0 ? 1.0 : std::expm1(remaining_steps * std::log(scalar_p)) / (scalar_p - 1.0));
            float next_value = static_cast<float>(power * state_[node] + scalar_q * sum);
            if (config_.clamp_normalized_state) next_value = clamp01(next_value);
            const float delta = std::abs(next_value - state_[node]);
            maximum_delta = delta;
            if (delta > 0.0F || !std::isfinite(next_value))
                pending_.push_back(PendingValue{node, next_value, scc, delta});
            ++stats.scalar_closed_form_executions;
            ++stats.analytic_updates;
            completed_horizon = true;
        } else if ((network_.linear_plans[scc].backend == TemporalBackend::PrecomputedDensePropagator ||
                    network_.linear_plans[scc].backend == TemporalBackend::RuntimeDensePropagator) &&
                   evaluate_dense_linear(scc, remaining_steps, stats)) {
            // Dense finite-horizon map has committed a batch of pending node values.
            ++stats.analytic_updates;
            completed_horizon = true;
        } else if (evaluate_sparse_power_action(scc, remaining_steps, stats)) {
            ++stats.analytic_updates;
            completed_horizon = true;
        } else {
            evaluate_sparse_linear(scc, stats);
        }
    }

    ++stats.scc_evaluations;
    if (plan.classification == SccTemporalClass::Static || completed_horizon)
        stats.avoided_integration_steps += remaining_steps > 0U ? remaining_steps - 1U : 0U;
    record_scc_visit(scc, stats);
    record_scc_class(scc, true, stats);
    if (completed_horizon) temporal_completed_generation_[scc] = current_temporal_generation_;
    static_cast<void>(max_steps);
    return maximum_delta;
}

float Engine::commit_pending(ExecutionStats& stats, std::vector<float>* scc_delta) {
    float maximum_delta = 0.0F;
    changed_nodes_this_step_.clear();
    for (const auto& pending : pending_) {
        if (!std::isfinite(pending.value)) {
            stats.numerical_error = true;
            continue;
        }
        const float delta = pending.delta;
        maximum_delta = std::max(maximum_delta, delta);
        if (scc_delta != nullptr) {
            (*scc_delta)[pending.scc] = std::max((*scc_delta)[pending.scc], delta);
        }
        if (delta > 0.0F) {
            set_state(pending.node, pending.value);
            changed_nodes_this_step_.push_back(pending.node);
        }
    }
    stats.maximum_step_delta = std::max(stats.maximum_step_delta, maximum_delta);
    return maximum_delta;
}

void Engine::remember_before_change(NodeId node) {
    if (!journaling_) {
        return;
    }
    if (journal_generation_[node] == current_journal_generation_) {
        return;
    }
    journal_generation_[node] = current_journal_generation_;
    journal_.push_back(JournalEntry{node, state_[node]});
}

void Engine::set_state(NodeId node, float value) {
    remember_before_change(node);
    state_[node] = value;
}

void Engine::apply_perturbation(const Perturbation& perturbation) {
    if (!std::isfinite(perturbation.strength)) {
        throw std::invalid_argument("perturbation strength must be finite");
    }
    seed_sccs_.clear();

    const bool reaction_target = perturbation.type == PerturbationType::ReactionInhibition ||
                                 perturbation.type == PerturbationType::ReactionActivation ||
                                 perturbation.type == PerturbationType::EdgeBlock;
    if (reaction_target) {
        if (perturbation.reaction >= network_.reactions.size()) {
            throw std::out_of_range("perturbation reaction is out of range");
        }
        if ((perturbation.type == PerturbationType::ReactionInhibition &&
             (perturbation.strength < 0.0F || perturbation.strength > 1.0F)) ||
            (perturbation.type == PerturbationType::ReactionActivation &&
             perturbation.strength < 0.0F)) {
            throw std::invalid_argument("reaction perturbation strength is out of range");
        }
        if (perturbation.type == PerturbationType::ReactionInhibition) {
            reaction_multiplier_[perturbation.reaction] = 1.0F - perturbation.strength;
        } else if (perturbation.type == PerturbationType::ReactionActivation) {
            reaction_multiplier_[perturbation.reaction] = 1.0F + perturbation.strength;
        } else {
            reaction_multiplier_[perturbation.reaction] = 0.0F;
        }
        const std::uint64_t begin = network_.reaction_scc_offsets[perturbation.reaction];
        const std::uint64_t end = network_.reaction_scc_offsets[perturbation.reaction + 1U];
        for (std::uint64_t index = begin; index < end; ++index) {
            seed_sccs_.push_back(network_.reaction_sccs[static_cast<std::size_t>(index)]);
        }
        return;
    }

    if (perturbation.node >= network_.nodes.size()) {
        throw std::out_of_range("perturbation node is out of range");
    }
    if (perturbation.strength < 0.0F || perturbation.strength > 1.0F) {
        throw std::invalid_argument("node perturbation strength must be in [0,1]");
    }

    const NodeId node = perturbation.node;
    switch (perturbation.type) {
        case PerturbationType::GeneKnockout:
            clamp_active_[node] = 1U;
            clamp_value_[node] = 0.0F;
            set_state(node, 0.0F);
            break;
        case PerturbationType::GeneKnockdown:
        case PerturbationType::ActivityInhibition:
            clamp_active_[node] = 1U;
            clamp_value_[node] = baseline_[node] * (1.0F - perturbation.strength);
            set_state(node, clamp_value_[node]);
            break;
        case PerturbationType::ActivityActivation:
            clamp_active_[node] = 1U;
            clamp_value_[node] =
                baseline_[node] + perturbation.strength * (1.0F - baseline_[node]);
            set_state(node, clamp_value_[node]);
            break;
        case PerturbationType::ProteinDegradation:
            extra_node_decay_[node] = perturbation.strength;
            break;
        case PerturbationType::Clamp:
            clamp_active_[node] = 1U;
            clamp_value_[node] = perturbation.strength;
            set_state(node, clamp_value_[node]);
            break;
        case PerturbationType::ReactionInhibition:
        case PerturbationType::ReactionActivation:
        case PerturbationType::EdgeBlock:
            break;
    }

    const SccId own_scc = network_.node_to_scc[node];
    seed_sccs_.push_back(own_scc);
    const auto& block = network_.sccs[own_scc];
    for (std::uint32_t index = 0; index < block.downstream_count; ++index) {
        seed_sccs_.push_back(network_.downstream_sccs[block.downstream_offset + index]);
    }
}

void Engine::rollback(
    const Perturbation& perturbation,
    const std::optional<std::pair<ReactionId, float>>& reaction_change) {
    for (auto iterator = journal_.rbegin(); iterator != journal_.rend(); ++iterator) {
        state_[iterator->node] = iterator->old_value;
    }
    journal_.clear();
    if (perturbation.node < network_.nodes.size()) {
        clamp_active_[perturbation.node] = 0U;
        extra_node_decay_[perturbation.node] = 0.0F;
    }
    if (reaction_change.has_value()) {
        reaction_multiplier_[reaction_change->first] = reaction_change->second;
    }
    journaling_ = false;
}

ExecutionStats Engine::run_full(std::uint32_t max_steps, float convergence_threshold,
                                const BaselineStepObserver& observer) {
    ExecutionStats stats;
    std::vector<float> per_scc_delta;
    if (observer) per_scc_delta.resize(network_.sccs.size());
    ++current_scc_visit_generation_;
    if (current_scc_visit_generation_ == 0U) {
        std::fill(scc_visit_generation_.begin(), scc_visit_generation_.end(), 0U);
        std::fill(analytic_visit_generation_.begin(), analytic_visit_generation_.end(), 0U);
        std::fill(numerical_visit_generation_.begin(), numerical_visit_generation_.end(), 0U);
        current_scc_visit_generation_ = 1U;
    }
    for (std::uint32_t step = 0; step < max_steps; ++step) {
        pending_.clear();
        if (observer) std::fill(per_scc_delta.begin(), per_scc_delta.end(), 0.0F);
        for (const SccId scc : network_.topological_sccs) {
            const float delta = evaluate_scc(scc, stats);
            if (observer) per_scc_delta[scc] = std::max(per_scc_delta[scc], delta);
        }
        const float delta = commit_pending(stats);
        ++stats.integration_steps;
        if (observer) observer(step, delta, per_scc_delta);
        if (stats.numerical_error) {
            break;
        }
        if (delta <= convergence_threshold) {
            stats.converged = true;
            break;
        }
    }
    return stats;
}

void Engine::enqueue_unique(std::vector<SccId>& queue, SccId scc) {
    if (queue_generation_[scc] == current_queue_generation_) {
        return;
    }
    queue_generation_[scc] = current_queue_generation_;
    queue.push_back(scc);
}

void Engine::record_scc_visit(SccId scc, ExecutionStats& stats) {
    if (scc_visit_generation_[scc] != current_scc_visit_generation_) {
        scc_visit_generation_[scc] = current_scc_visit_generation_;
        ++stats.visited_unique_scc;
    }
}

void Engine::record_scc_class(SccId scc, bool analytic, ExecutionStats& stats) {
    auto& generations = analytic ? analytic_visit_generation_ : numerical_visit_generation_;
    if (generations[scc] == current_scc_visit_generation_) return;
    generations[scc] = current_scc_visit_generation_;
    if (analytic) {
        ++stats.analytic_scc_count;
    } else {
        ++stats.numerical_scc_count;
        const auto kind = network_.temporal_plans[scc].classification;
        if (kind == SccTemporalClass::LinearHomogeneous || kind == SccTemporalClass::LinearAffine)
            ++stats.linear_fallback_scc_count;
    }
}

void Engine::record_reaction_evaluation(SccId scc, ExecutionStats& stats) {
    ++stats.reaction_evaluations;
    const SccTemporalClass classification = network_.temporal_plans[scc].classification;
    if (classification == SccTemporalClass::SimpleNonlinear ||
        classification == SccTemporalClass::GenericNonlinear) {
        ++stats.nonlinear_reaction_evaluations;
    }
}

void Engine::begin_dirty_generation() {
    ++current_dirty_generation_;
    if (current_dirty_generation_ == 0U) {
        std::fill(dirty_reaction_generation_.begin(), dirty_reaction_generation_.end(), 0U);
        std::fill(dirty_count_generation_.begin(), dirty_count_generation_.end(), 0U);
        current_dirty_generation_ = 1U;
    }
}

void Engine::mark_reaction_slot_dirty(SccId scc, std::uint32_t slot) {
    const SccBlock& block = network_.sccs[scc];
    if (slot >= block.reaction_count) return;
    const std::uint32_t flat_slot = block.reaction_offset + slot;
    if (dirty_reaction_generation_[flat_slot] != current_dirty_generation_) {
        dirty_reaction_generation_[flat_slot] = current_dirty_generation_;
        if (dirty_count_generation_[scc] != current_dirty_generation_) {
            dirty_count_generation_[scc] = current_dirty_generation_;
            dirty_counts_[scc] = 0U;
        }
        ++dirty_counts_[scc];
    }
}

void Engine::mark_scc_all_dirty(SccId scc) {
    const SccBlock& block = network_.sccs[scc];
    for (std::uint32_t slot = 0; slot < block.reaction_count; ++slot) {
        mark_reaction_slot_dirty(scc, slot);
    }
}

void Engine::mark_node_dependents_dirty(NodeId node) {
    const SccId scc = network_.node_to_scc[node];
    const SccBlock& block = network_.sccs[scc];
    const std::uint32_t flat_node = block.node_offset + network_.node_position_in_scc[node];
    const std::uint64_t begin = network_.scc_node_reaction_offsets[flat_node];
    const std::uint64_t end = network_.scc_node_reaction_offsets[flat_node + 1U];
    for (std::uint64_t index = begin; index < end; ++index) {
        const std::uint32_t flat_slot =
            network_.scc_node_reaction_slots[static_cast<std::size_t>(index)];
        mark_reaction_slot_dirty(scc, flat_slot - block.reaction_offset);
    }
}

ExecutionStats Engine::run_frontier(const std::vector<SccId>& seeds, std::uint32_t max_steps) {
    ExecutionStats stats;
    // Boundary traces are reserved for the configured maximum horizon, but a caller may
    // request a shorter finite run (notably a scheduled-protocol segment). Never execute
    // deferred linear regions past that actual segment horizon.
    const std::uint32_t trace_horizon = std::min(boundary_trace_horizon_, max_steps);
    // A per-step delta can accumulate over the remaining finite integration horizon. Scale the
    // scheduling threshold by dt so epsilon remains an approximate state-error tolerance.
    const float frontier_threshold = config_.epsilon * std::min(config_.dt, 1.0F);
    current_frontier_.clear();
    if (++current_temporal_generation_ == 0U) {
        std::fill(temporal_completed_generation_.begin(), temporal_completed_generation_.end(), 0U);
        current_temporal_generation_ = 1U;
    }
    if (++current_boundary_generation_ == 0U) {
        std::fill(boundary_previous_generation_.begin(), boundary_previous_generation_.end(), 0U);
        current_boundary_generation_ = 1U;
    }
    begin_dirty_generation();
    ++current_scc_visit_generation_;
    if (current_scc_visit_generation_ == 0U) {
        std::fill(scc_visit_generation_.begin(), scc_visit_generation_.end(), 0U);
        std::fill(analytic_visit_generation_.begin(), analytic_visit_generation_.end(), 0U);
        std::fill(numerical_visit_generation_.begin(), numerical_visit_generation_.end(), 0U);
        current_scc_visit_generation_ = 1U;
    }
    ++current_queue_generation_;
    if (current_queue_generation_ == 0U) {
        std::fill(queue_generation_.begin(), queue_generation_.end(), 0U);
        current_queue_generation_ = 1U;
    }
    for (const SccId seed : seeds) {
        enqueue_unique(current_frontier_, seed);
        mark_scc_all_dirty(seed);
    }
    stats.queue_pushes += current_frontier_.size();

    std::uint32_t recorded_steps = 0U;
    for (std::uint32_t step = 0; step < max_steps && !current_frontier_.empty(); ++step) {
        if (step < trace_horizon) {
            const auto recording_start = std::chrono::steady_clock::now();
            const std::size_t row = static_cast<std::size_t>(step) * boundary_trace_nodes_.size();
            for (std::uint32_t i = 0; i < boundary_trace_nodes_.size(); ++i)
                boundary_trajectory_[row + i] = state_[boundary_trace_nodes_[i]];
            stats.time_boundary_recording_us += std::chrono::duration<double, std::micro>(
                std::chrono::steady_clock::now() - recording_start).count();
            recorded_steps = step + 1U;
        }
        pending_.clear();
        const std::uint64_t numerical_steps_before = stats.numerical_steps;
        for (const SccId scc : current_frontier_) {
            if (temporal_completed_generation_[scc] == current_temporal_generation_) continue;
            if (scc < boundary_deferred_scc_.size() && boundary_deferred_scc_[scc]) continue;
            const std::uint32_t remaining_steps = max_steps - step;
            if (can_evaluate_scc_temporal(scc, remaining_steps, max_steps)) {
                static_cast<void>(evaluate_scc_temporal(scc, remaining_steps, max_steps, stats));
            } else if (config_.mode == SimulationMode::Continuous &&
                       !network_.linear_plans.empty() &&
                       (network_.temporal_plans[scc].classification == SccTemporalClass::LinearAffine ||
                        network_.temporal_plans[scc].classification == SccTemporalClass::LinearHomogeneous) &&
                       network_.linear_plans[scc].backend != TemporalBackend::NumericalNonlinear) {
                evaluate_sparse_linear(scc, stats);
            } else if ((dirty_count_generation_[scc] == current_dirty_generation_ && dirty_counts_[scc] > 0U) ||
                       network_.sccs[scc].reaction_count == 0U) {
                static_cast<void>(evaluate_scc_dirty(scc, stats));
            }
        }
        static_cast<void>(commit_pending(stats));
        if (stats.numerical_steps != numerical_steps_before) ++stats.integration_steps;
        if (stats.numerical_error) {
            break;
        }

        next_frontier_.clear();
        begin_dirty_generation();
        ++current_queue_generation_;
        if (current_queue_generation_ == 0U) {
            std::fill(queue_generation_.begin(), queue_generation_.end(), 0U);
            current_queue_generation_ = 1U;
        }
        for (const auto& pending : pending_) {
            if (pending.delta <= frontier_threshold) {
                continue;
            }
            const SccId changed = pending.scc;
            if (temporal_completed_generation_[changed] != current_temporal_generation_) {
                enqueue_unique(next_frontier_, changed);
                mark_node_dependents_dirty(pending.node);
                if (extra_node_decay_[pending.node] > 0.0F) {
                    mark_scc_all_dirty(changed);
                }
            }
            const auto& block = network_.sccs[changed];
            for (std::uint32_t downstream = 0; downstream < block.downstream_count; ++downstream) {
                const SccId downstream_scc =
                    network_.downstream_sccs[block.downstream_offset + downstream];
                enqueue_unique(next_frontier_, downstream_scc);
                mark_scc_all_dirty(downstream_scc);
            }
        }
        stats.queue_pushes += next_frontier_.size();
        current_frontier_.swap(next_frontier_);
    }
    if (trace_horizon > 0U && recorded_steps > 0U) {
        const auto recording_start = std::chrono::steady_clock::now();
        for (std::uint32_t step = recorded_steps; step < trace_horizon; ++step) {
            const std::size_t row = static_cast<std::size_t>(step) * boundary_trace_nodes_.size();
            for (std::uint32_t i = 0; i < boundary_trace_nodes_.size(); ++i)
                boundary_trajectory_[row + i] = state_[boundary_trace_nodes_[i]];
        }
        stats.time_boundary_recording_us += std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - recording_start).count();
        execute_deferred_boundary_sccs(trace_horizon, stats);
    }
    stats.converged = current_frontier_.empty() && !stats.numerical_error;
    return stats;
}

void Engine::reserve_result_buffers(SimulationResult& result) const {
    result.molecular_changes.reserve(network_.nodes.size());
    result.phenotype_changes.reserve(network_.phenotype_nodes.size());
}

void Engine::collect_result_into(const ExecutionStats& stats,
                                 double elapsed_us,
                                 SimulationResult& result) const {
    result.execution_time_us = elapsed_us;
    result.cache_hit = false;
    result.stats = stats;
    result.molecular_changes.clear();
    result.phenotype_changes.clear();
    for (const JournalEntry& change : journal_) {
        const std::size_t index = change.node;
        const float delta = state_[index] - baseline_[index];
        if (std::abs(delta) <= config_.result_epsilon) {
            continue;
        }
        const NodeId node = change.node;
        result.molecular_changes.push_back(NodeChange{node, baseline_[index], state_[index], delta});
        if (network_.nodes[index].type == NodeType::Phenotype) {
            result.phenotype_changes.push_back(PhenotypeChange{node, delta});
        }
    }
    std::sort(result.molecular_changes.begin(), result.molecular_changes.end(),
              [](const NodeChange& lhs, const NodeChange& rhs) { return lhs.id < rhs.id; });
    std::sort(result.phenotype_changes.begin(), result.phenotype_changes.end(),
              [](const PhenotypeChange& lhs, const PhenotypeChange& rhs) { return lhs.id < rhs.id; });
    result.stats.changed_nodes = result.molecular_changes.size();
}

SimulationResult Engine::run(const Perturbation& perturbation, ExecutionStrategy strategy) {
    SimulationResult result;
    constexpr std::size_t small_result_reserve = 64U;
    result.molecular_changes.reserve(std::min(network_.nodes.size(), small_result_reserve));
    result.phenotype_changes.reserve(std::min(network_.phenotype_nodes.size(), small_result_reserve));
    run_into(perturbation, result, strategy);
    return result;
}

void Engine::run_into(const Perturbation& perturbation,
                      SimulationResult& output,
                      ExecutionStrategy strategy) {
    if (baseline_.size() != network_.nodes.size()) {
        throw std::logic_error("initialize_baseline() must be called before run()");
    }

    const CacheKey cache_key = make_cache_key(perturbation, strategy);
    if (config_.enable_hot_cache) {
        const auto cache_start = std::chrono::steady_clock::now();
        const auto cached = hot_cache_.find(cache_key);
        if (cached != hot_cache_.end()) {
            output = cached->second;
            output.cache_hit = true;
            output.execution_time_us = std::chrono::duration<double, std::micro>(
                                           std::chrono::steady_clock::now() - cache_start)
                                           .count();
            return;
        }
    }

    ++current_journal_generation_;
    if (current_journal_generation_ == 0U) {
        std::fill(journal_generation_.begin(), journal_generation_.end(), 0U);
        current_journal_generation_ = 1U;
    }
    journal_.clear();
    journaling_ = true;

    std::optional<std::pair<ReactionId, float>> reaction_change;
    if ((perturbation.type == PerturbationType::ReactionInhibition ||
         perturbation.type == PerturbationType::ReactionActivation ||
         perturbation.type == PerturbationType::EdgeBlock) &&
        perturbation.reaction < reaction_multiplier_.size()) {
        reaction_change = std::make_pair(perturbation.reaction,
                                         reaction_multiplier_[perturbation.reaction]);
    }

    try {
        const auto start = std::chrono::steady_clock::now();
        apply_perturbation(perturbation);
        active_perturbed_reaction_ = reaction_change.has_value()
            ? perturbation.reaction : kInvalidReaction;
        active_perturbed_node_ = perturbation.type == PerturbationType::ProteinDegradation
            ? perturbation.node : kInvalidNode;
        ExecutionStats stats = strategy == ExecutionStrategy::Full
                                   ? run_full(config_.perturbation_steps, 0.0F)
                                   : run_frontier(seed_sccs_, config_.perturbation_steps);
        const auto finish = std::chrono::steady_clock::now();
        const double elapsed =
            std::chrono::duration<double, std::micro>(finish - start).count();
        const auto materialization_start = std::chrono::steady_clock::now();
        collect_result_into(stats, elapsed, output);
        output.stats.time_result_materialization_us = std::chrono::duration<double, std::micro>(
            std::chrono::steady_clock::now() - materialization_start).count();
        rollback(perturbation, reaction_change);
        active_perturbed_reaction_ = kInvalidReaction;
        active_perturbed_node_ = kInvalidNode;
        if (config_.enable_hot_cache && config_.hot_cache_max_entries > 0U) {
            if (hot_cache_.size() >= config_.hot_cache_max_entries) {
                hot_cache_.clear();
            }
            hot_cache_.emplace(cache_key, output);
        }
    } catch (...) {
        rollback(perturbation, reaction_change);
        active_perturbed_reaction_ = kInvalidReaction;
        active_perturbed_node_ = kInvalidNode;
        throw;
    }
}

std::vector<SimulationResult> Engine::run_batch(
    const std::vector<Perturbation>& perturbations,
    ExecutionStrategy strategy) {
    std::vector<SimulationResult> results;
    results.reserve(perturbations.size());
    for (const auto& perturbation : perturbations) {
        results.push_back(run(perturbation, strategy));
    }
    return results;
}

EngineStateSnapshot Engine::snapshot_state() const {
    return EngineStateSnapshot{state_, reaction_multiplier_, extra_node_decay_, clamp_value_, clamp_active_};
}

void Engine::restore_state(const EngineStateSnapshot& snapshot) {
    if (snapshot.state.size() != network_.nodes.size() ||
        snapshot.reaction_multipliers.size() != network_.reactions.size() ||
        snapshot.extra_node_decay.size() != network_.nodes.size() ||
        snapshot.clamp_values.size() != network_.nodes.size() ||
        snapshot.clamp_active.size() != network_.nodes.size()) {
        throw std::invalid_argument("engine snapshot dimensions do not match the compiled network");
    }
    const auto finite = [](const std::vector<float>& values) {
        return std::all_of(values.begin(), values.end(), [](float value) { return std::isfinite(value); });
    };
    if (!finite(snapshot.state) || !finite(snapshot.reaction_multipliers) ||
        !finite(snapshot.extra_node_decay) || !finite(snapshot.clamp_values)) {
        throw std::invalid_argument("engine snapshot contains a non-finite value");
    }
    state_ = snapshot.state;
    baseline_ = snapshot.state;
    reaction_multiplier_ = snapshot.reaction_multipliers;
    extra_node_decay_ = snapshot.extra_node_decay;
    clamp_value_ = snapshot.clamp_values;
    clamp_active_ = snapshot.clamp_active;
    journal_.clear();
    journaling_ = false;
    active_perturbed_reaction_ = kInvalidReaction;
    active_perturbed_node_ = kInvalidNode;
    ++baseline_state_id_;
    clear_cache();
}

void Engine::reserve_protocol_result(ProtocolResult& result,
                                     std::size_t sample_count,
                                     std::size_t intervention_count) const {
    result.samples.clear();
    result.samples.resize(sample_count);
    for (auto& sample : result.samples) sample.state.resize(network_.nodes.size());
    result.segments.clear();
    result.segments.reserve(sample_count + intervention_count + 1U);
}

void Engine::run_protocol_into(std::span<const TimedIntervention> interventions,
                               double end_time,
                               std::span<const double> sample_times,
                               ProtocolResult& result,
                               ExecutionStrategy strategy) {
    if (baseline_.size() != network_.nodes.size())
        throw std::logic_error("initialize_baseline() must be called before run_protocol_into()");
    if (!std::isfinite(end_time) || end_time < 0.0)
        throw std::invalid_argument("protocol end_time must be finite and non-negative");
    const double dt = config_.dt;
    const auto time_to_step = [dt](double time) -> std::uint64_t {
        if (!std::isfinite(time) || time < 0.0)
            throw std::invalid_argument("protocol times must be finite and non-negative");
        const double raw = time / dt;
        const double nearest = std::round(raw);
        if (std::abs(raw - nearest) > 1.0e-6 * std::max(1.0, std::abs(raw)))
            throw std::invalid_argument("protocol event/sample times must align with solver dt");
        if (nearest > static_cast<double>(std::numeric_limits<std::uint32_t>::max()))
            throw std::out_of_range("protocol time exceeds the supported step range");
        return static_cast<std::uint64_t>(nearest);
    };
    const std::uint64_t end_step = time_to_step(end_time);
    if (end_step > std::numeric_limits<std::uint32_t>::max())
        throw std::out_of_range("protocol exceeds the supported total step range");
    std::uint64_t previous = 0U;
    for (const auto& event : interventions) {
        const auto step = time_to_step(event.time);
        if (step > end_step || step < previous)
            throw std::invalid_argument("interventions must be ordered and occur no later than end_time");
        previous = step;
    }
    previous = 0U;
    for (const double time : sample_times) {
        const auto step = time_to_step(time);
        if (step > end_step || step < previous)
            throw std::invalid_argument("sample_times must be ordered and no later than end_time");
        previous = step;
    }
    if (result.samples.size() != sample_times.size() ||
        result.segments.capacity() < sample_times.size() + interventions.size() + 1U)
        throw std::invalid_argument("protocol result buffers are not reserved for this schedule");
    for (const auto& sample : result.samples) {
        if (sample.state.size() != network_.nodes.size())
            throw std::invalid_argument("protocol result buffers were reserved for a different network");
    }
    result.segments.clear();
    result.numerical_error = false;
    std::copy(state_.begin(), state_.end(), protocol_saved_state_.begin());
    std::copy(baseline_.begin(), baseline_.end(), protocol_saved_baseline_.begin());
    std::copy(reaction_multiplier_.begin(), reaction_multiplier_.end(),
              protocol_saved_reaction_multiplier_.begin());
    std::copy(extra_node_decay_.begin(), extra_node_decay_.end(),
              protocol_saved_extra_node_decay_.begin());
    std::copy(clamp_value_.begin(), clamp_value_.end(), protocol_saved_clamp_value_.begin());
    std::copy(clamp_active_.begin(), clamp_active_.end(), protocol_saved_clamp_active_.begin());
    std::fill(reaction_multiplier_.begin(), reaction_multiplier_.end(), 1.0F);
    std::fill(extra_node_decay_.begin(), extra_node_decay_.end(), 0.0F);
    std::fill(clamp_active_.begin(), clamp_active_.end(), 0U);
    std::fill(clamp_value_.begin(), clamp_value_.end(), 0.0F);
    std::copy(protocol_saved_baseline_.begin(), protocol_saved_baseline_.end(), state_.begin());
    std::copy(state_.begin(), state_.end(), baseline_.begin());
    journal_.clear();
    journaling_ = false;
    active_perturbed_reaction_ = kInvalidReaction;
    active_perturbed_node_ = kInvalidNode;
    protocol_seeds_.clear();
    protocol_causal_queue_.clear();
    if (++current_protocol_seed_generation_ == 0U) {
        std::fill(protocol_seed_generation_.begin(), protocol_seed_generation_.end(), 0U);
        current_protocol_seed_generation_ = 1U;
    }
    const auto add_protocol_seed = [&](SccId scc) {
        if (protocol_seed_generation_[scc] == current_protocol_seed_generation_) return;
        protocol_seed_generation_[scc] = current_protocol_seed_generation_;
        protocol_seeds_.push_back(scc);
        protocol_causal_queue_.push_back(scc);
    };
    std::size_t event_index = 0U, sample_index = 0U;
    std::uint64_t current_step = 0U;
    const auto capture_samples = [&](std::uint64_t at_step, std::size_t& index) {
        while (index < sample_times.size() && time_to_step(sample_times[index]) == at_step) {
            auto& sample = result.samples[index];
            sample.time = static_cast<double>(at_step) * dt;
            std::copy(state_.begin(), state_.end(), sample.state.begin());
            ++index;
        }
    };
    const auto apply_events = [&](std::uint64_t at_step, std::size_t& index) {
        while (index < interventions.size() && time_to_step(interventions[index].time) == at_step) {
            baseline_ = state_;
            seed_sccs_.clear();
            apply_perturbation(interventions[index].action);
            for (const SccId seed : seed_sccs_) add_protocol_seed(seed);
            // Sample boundaries restart the finite Frontier scheduler. Re-seed the full
            // causal closure so an SCC that slept in the previous segment can continue
            // its transient from the checkpointed state.
            for (std::size_t q = 0; q < protocol_causal_queue_.size(); ++q) {
                const SccId scc = protocol_causal_queue_[q];
                const auto& block = network_.sccs[scc];
                for (std::uint32_t child = 0; child < block.downstream_count; ++child) {
                    const SccId downstream = network_.downstream_sccs[block.downstream_offset + child];
                    add_protocol_seed(downstream);
                }
            }
            protocol_causal_queue_.clear();
            std::sort(protocol_seeds_.begin(), protocol_seeds_.end(), [&](SccId left, SccId right) {
                return protocol_topological_rank_[left] < protocol_topological_rank_[right];
            });
            const auto type = interventions[index].action.type;
            if (type == PerturbationType::ReactionInhibition ||
                type == PerturbationType::ReactionActivation || type == PerturbationType::EdgeBlock)
                active_perturbed_reaction_ = interventions[index].action.reaction;
            if (type == PerturbationType::ProteinDegradation)
                active_perturbed_node_ = interventions[index].action.node;
            ++index;
        }
    };
    try {
        apply_events(0U, event_index);
        capture_samples(0U, sample_index);
        while (current_step < end_step) {
            std::uint64_t next_step = end_step;
            if (event_index < interventions.size())
                next_step = std::min(next_step, time_to_step(interventions[event_index].time));
            if (sample_index < sample_times.size())
                next_step = std::min(next_step, time_to_step(sample_times[sample_index]));
            if (next_step < current_step)
                throw std::logic_error("protocol scheduler moved backwards in time");
            if (next_step > current_step) {
                baseline_ = state_;
                ++baseline_state_id_;
                const double start_time = static_cast<double>(current_step) * dt;
                const double finish_time = static_cast<double>(next_step) * dt;
                ExecutionStats stats = strategy == ExecutionStrategy::Full
                    ? run_full(static_cast<std::uint32_t>(next_step - current_step), 0.0F)
                    : run_frontier(protocol_seeds_,
                                   static_cast<std::uint32_t>(next_step - current_step));
                result.numerical_error = result.numerical_error || stats.numerical_error;
                result.segments.push_back(ProtocolSegment{start_time, finish_time, stats});
                if (strategy == ExecutionStrategy::Frontier) {
                    for (const SccId active : current_frontier_) {
                        add_protocol_seed(active);
                    }
                    protocol_causal_queue_.clear();
                    std::sort(protocol_seeds_.begin(), protocol_seeds_.end(), [&](SccId left, SccId right) {
                        return protocol_topological_rank_[left] < protocol_topological_rank_[right];
                    });
                }
                current_step = next_step;
            }
            apply_events(current_step, event_index);
            capture_samples(current_step, sample_index);
            if (result.numerical_error) break;
        }
        // Events exactly at end_time are meaningful, even though no integration follows.
        apply_events(end_step, event_index);
        capture_samples(end_step, sample_index);
        if (sample_index != sample_times.size())
            throw std::logic_error("protocol scheduler did not produce every requested sample");
        std::copy(state_.begin(), state_.end(), baseline_.begin());
        ++baseline_state_id_;
        clear_cache();
    } catch (...) {
        std::copy(protocol_saved_state_.begin(), protocol_saved_state_.end(), state_.begin());
        std::copy(protocol_saved_baseline_.begin(), protocol_saved_baseline_.end(), baseline_.begin());
        std::copy(protocol_saved_reaction_multiplier_.begin(), protocol_saved_reaction_multiplier_.end(),
                  reaction_multiplier_.begin());
        std::copy(protocol_saved_extra_node_decay_.begin(), protocol_saved_extra_node_decay_.end(),
                  extra_node_decay_.begin());
        std::copy(protocol_saved_clamp_value_.begin(), protocol_saved_clamp_value_.end(), clamp_value_.begin());
        std::copy(protocol_saved_clamp_active_.begin(), protocol_saved_clamp_active_.end(), clamp_active_.begin());
        journal_.clear();
        journaling_ = false;
        active_perturbed_reaction_ = kInvalidReaction;
        active_perturbed_node_ = kInvalidNode;
        ++baseline_state_id_;
        clear_cache();
        throw;
    }
}

std::size_t Engine::CacheKeyHash::operator()(const CacheKey& key) const noexcept {
    std::size_t hash = 0xcbf29ce484222325ULL;
    auto mix = [&hash](std::uint64_t value) {
        hash ^= static_cast<std::size_t>(value + 0x9e3779b97f4a7c15ULL +
                                         (static_cast<std::uint64_t>(hash) << 6U) +
                                         (static_cast<std::uint64_t>(hash) >> 2U));
    };
    mix(key.network_version);
    mix(key.context_id);
    mix(key.baseline_state_id);
    mix(key.perturbation_type);
    mix(key.node);
    mix(key.reaction);
    mix(key.strength_bits);
    mix(key.strategy);
    return hash;
}

Engine::CacheKey Engine::make_cache_key(
    const Perturbation& perturbation,
    ExecutionStrategy strategy) const {
    CacheKey key;
    key.network_version = network_.network_version;
    key.context_id = context_id_;
    key.baseline_state_id = baseline_state_id_;
    key.perturbation_type = static_cast<std::uint32_t>(perturbation.type);
    key.node = perturbation.node;
    key.reaction = perturbation.reaction;
    key.strength_bits = std::bit_cast<std::uint32_t>(perturbation.strength);
    key.strategy = static_cast<std::uint32_t>(strategy);
    return key;
}

}  // namespace cellnet
