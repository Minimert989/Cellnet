#pragma once

#include "cellnet/context.hpp"
#include "cellnet/network.hpp"

#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

namespace cellnet {

struct SolverConfig {
    SimulationMode mode{SimulationMode::Continuous};
    float dt{0.05F};
    float epsilon{1.0e-6F};
    float result_epsilon{1.0e-5F};
    std::uint32_t baseline_max_steps{2'000};
    std::uint32_t perturbation_steps{100};
    bool clamp_normalized_state{true};
    bool apply_confidence{true};
    bool enable_hot_cache{false};
    bool enable_boundary_compression{true};
    std::size_t hot_cache_max_entries{1'024U};
    bool enable_dense_temporal{true};
};

struct Perturbation {
    PerturbationType type{PerturbationType::Clamp};
    NodeId node{kInvalidNode};
    ReactionId reaction{kInvalidReaction};
    float strength{1.0F};
};

struct NodeChange {
    NodeId id{kInvalidNode};
    float baseline{0.0F};
    float value{0.0F};
    float delta{0.0F};
};

struct PhenotypeChange {
    NodeId id{kInvalidNode};
    float delta{0.0F};
};

struct ExecutionStats {
    std::uint64_t visited_unique_scc{0};
    std::uint64_t scc_evaluations{0};
    std::uint64_t changed_nodes{0};
    std::uint64_t reaction_evaluations{0};
    std::uint64_t queue_pushes{0};
    // Frontier solver scheduler steps that performed at least one numerical SCC
    // evaluation; analytic linear power actions do not count as integration steps.
    std::uint64_t integration_steps{0};
    std::uint64_t analytic_scc_count{0};
    std::uint64_t numerical_scc_count{0};
    std::uint64_t analytic_updates{0};
    std::uint64_t numerical_steps{0};
    std::uint64_t avoided_integration_steps{0};
    std::uint64_t nonlinear_reaction_evaluations{0};
    std::uint64_t linear_fallback_scc_count{0};
    std::uint64_t static_scc_executions{0}, scalar_closed_form_executions{0};
    std::uint64_t precomputed_propagator_executions{0}, runtime_propagator_builds{0};
    std::uint64_t runtime_propagator_cache_hits{0}, boundary_driven_linear_executions{0};
    std::uint64_t boundary_input_steps{0}, boundary_segments{0};
    std::uint64_t boundary_steps{0}, response_kernel_applications{0};
    std::uint64_t fused_linear_regions{0}, avoided_linear_integration_steps{0};
    std::uint64_t sparse_linear_executions{0}, linear_operator_applications{0};
    std::uint64_t sparse_power_action_executions{0}, power_action_terms{0};
    std::uint64_t fused_region_count{0};
    std::uint64_t region_nodes_total{0}, region_nodes_touched{0};
    std::uint64_t region_edges_total{0}, region_edges_touched{0};
    std::uint64_t boundary_channels{0}, power_action_calls{0};
    std::uint64_t sparse_operator_passes{0}, spmv_equivalent_operations{0};
    std::uint64_t causal_cone_nodes{0}, causal_cone_edges{0};
    std::uint64_t boundary_forcing_aggregations{0}, region_passes_per_perturbation{0};
    std::uint64_t causal_cone_cache_hits{0}, causal_cone_cache_misses{0};
    std::uint64_t causal_cone_full_region_guard{0};
    std::uint64_t selected_sparse_recurrence_runs{0}, selected_sparse_power_action_runs{0};
    double time_boundary_recording_us{0.0}, time_boundary_compression_us{0.0};
    double time_region_execution_us{0.0}, time_power_action_us{0.0};
    double time_result_materialization_us{0.0};
    bool backend_selected_by_cost_model{false};
    std::uint64_t numerical_nonlinear_scc_count{0}, numerical_nonlinear_steps{0};
    std::uint64_t propagator_memory_bytes{0}, temporal_plan_memory_bytes{0};
    double trajectory_compression_ratio{1};
    float maximum_step_delta{0.0F};
    bool converged{false};
    bool numerical_error{false};
};

struct SimulationResult {
    double execution_time_us{0.0};
    bool cache_hit{false};
    ExecutionStats stats;
    std::vector<NodeChange> molecular_changes;
    std::vector<PhenotypeChange> phenotype_changes;
};

struct TimedIntervention {
    double time{0.0};
    Perturbation action;
};

struct ProtocolSample {
    double time{0.0};
    std::vector<float> state;
};

struct ProtocolSegment {
    double start_time{0.0};
    double end_time{0.0};
    ExecutionStats stats;
};

struct ProtocolResult {
    std::vector<ProtocolSample> samples;
    std::vector<ProtocolSegment> segments;
    bool numerical_error{false};
};

// Opaque value snapshot sufficient to restore dynamic state and persistent
// perturbation modifiers. Context is owned by the Engine and is not copied.
struct EngineStateSnapshot {
    std::vector<float> state;
    std::vector<float> reaction_multipliers;
    std::vector<float> extra_node_decay;
    std::vector<float> clamp_values;
    std::vector<std::uint8_t> clamp_active;
};

struct EngineMemoryBreakdown {
    std::size_t base_graph_bytes{0};
    std::size_t scc_metadata_bytes{0};
    std::size_t temporal_plan_bytes{0};
    std::size_t fused_region_bytes{0};
    std::size_t response_workspace_bytes{0};
    std::size_t scratch_workspace_bytes{0};
    std::size_t cache_bytes{0};

    [[nodiscard]] std::size_t total_bytes() const noexcept {
        return base_graph_bytes + scc_metadata_bytes + temporal_plan_bytes + fused_region_bytes +
               response_workspace_bytes + scratch_workspace_bytes + cache_bytes;
    }
};

class Engine {
public:
    using BaselineStepObserver =
        std::function<void(std::uint32_t, float, std::span<const float>)>;

    explicit Engine(CompiledNetwork network, SolverConfig config = {});

    void initialize_baseline(const BaselineStepObserver& observer = {});
    void set_context(const Context& context);
    void clear_cache();

    [[nodiscard]] SimulationResult run(
        const Perturbation& perturbation,
        ExecutionStrategy strategy = ExecutionStrategy::Frontier);

    // Reserve result storage once before a repeated workload. With hot-cache disabled and
    // adequately reserved output vectors, run_into performs no heap allocation during a
    // perturbation execution.
    void reserve_result_buffers(SimulationResult& result) const;

    // Execute into caller-owned storage. This is the allocation-free perturbation API when
    // reserve_result_buffers() has been called before the timed loop.
    void run_into(const Perturbation& perturbation,
                  SimulationResult& result,
                  ExecutionStrategy strategy = ExecutionStrategy::Frontier);

    [[nodiscard]] std::vector<SimulationResult> run_batch(
        const std::vector<Perturbation>& perturbations,
        ExecutionStrategy strategy = ExecutionStrategy::Frontier);

    [[nodiscard]] EngineStateSnapshot snapshot_state() const;
    void restore_state(const EngineStateSnapshot& snapshot);
    void reserve_protocol_result(ProtocolResult& result,
                                 std::size_t sample_count,
                                 std::size_t intervention_count = 0U) const;
    void run_protocol_into(std::span<const TimedIntervention> interventions,
                           double end_time,
                           std::span<const double> sample_times,
                           ProtocolResult& result,
                           ExecutionStrategy strategy = ExecutionStrategy::Frontier);

    [[nodiscard]] const CompiledNetwork& network() const noexcept { return network_; }
    [[nodiscard]] const std::vector<float>& baseline_state() const noexcept { return baseline_; }
    [[nodiscard]] const SolverConfig& config() const noexcept { return config_; }
    [[nodiscard]] EngineMemoryBreakdown memory_breakdown() const noexcept;
    [[nodiscard]] double region_construction_time_us() const noexcept {
        return region_construction_time_us_;
    }

private:
    struct PendingValue {
        NodeId node{kInvalidNode};
        float value{0.0F};
        SccId scc{kInvalidScc};
        float delta{0.0F};
    };

    struct JournalEntry {
        NodeId node{kInvalidNode};
        float old_value{0.0F};
    };

    CompiledNetwork network_;
    SolverConfig config_;
    std::vector<float> state_;
    std::vector<float> baseline_;
    std::vector<float> reaction_multiplier_;
    std::vector<float> context_reaction_multiplier_;
    std::vector<float> context_node_multiplier_;
    std::vector<float> extra_node_decay_;
    std::vector<float> clamp_value_;
    std::vector<std::uint8_t> clamp_active_;
    std::vector<std::uint32_t> journal_generation_;
    std::vector<std::uint32_t> queue_generation_;
    std::uint32_t current_journal_generation_{1};
    std::uint32_t current_queue_generation_{1};
    std::vector<JournalEntry> journal_;
    std::vector<PendingValue> pending_;
    std::vector<SccId> current_frontier_;
    std::vector<SccId> next_frontier_;
    std::vector<float> scratch_values_;
    std::vector<std::uint8_t> scratch_flags_;
    std::vector<SccId> seed_sccs_;
    // All protocol scratch/checkpoint buffers are sized at Engine construction so the
    // caller can repeatedly execute run_protocol_into() without allocating.
    std::vector<float> protocol_saved_state_, protocol_saved_baseline_;
    std::vector<float> protocol_saved_reaction_multiplier_, protocol_saved_extra_node_decay_;
    std::vector<float> protocol_saved_clamp_value_;
    std::vector<std::uint8_t> protocol_saved_clamp_active_;
    std::vector<SccId> protocol_seeds_, protocol_causal_queue_;
    std::vector<std::uint32_t> protocol_topological_rank_, protocol_seed_generation_;
    std::uint32_t current_protocol_seed_generation_{1U};
    std::vector<float> active_deltas_;
    std::vector<std::uint32_t> dirty_counts_;
    std::vector<std::uint32_t> dirty_count_generation_;
    std::vector<std::uint32_t> dirty_reaction_generation_;
    std::vector<std::uint32_t> scc_visit_generation_;
    std::vector<std::uint32_t> analytic_visit_generation_;
    std::vector<std::uint32_t> numerical_visit_generation_;
    std::vector<std::uint32_t> temporal_completed_generation_;
    std::uint32_t current_temporal_generation_{0};
    std::vector<NodeId> changed_nodes_this_step_;
    std::vector<double> linear_coefficients_, linear_scratch_;
    std::vector<double> power_action_a_, power_action_b_, power_action_result_;
    std::vector<double> power_action_diagonal_, power_action_offdiagonal_, power_action_rowsum_;
    std::vector<std::uint32_t> linear_generation_, boundary_generation_;
    std::vector<float> boundary_previous_;
    std::vector<std::uint32_t> boundary_previous_generation_;
    std::uint32_t current_boundary_generation_{0};
    std::vector<std::uint8_t> boundary_deferred_scc_;
    std::vector<NodeId> boundary_trace_nodes_;
    std::vector<std::uint32_t> boundary_trace_index_;
    std::vector<float> boundary_trajectory_, boundary_saved_state_;
    std::uint32_t boundary_trace_horizon_{0};
    struct LinearTemporalRegion {
        std::uint32_t scc_offset{0}, scc_count{0};
        std::uint32_t node_offset{0}, node_count{0};
        std::uint32_t entry_offset{0}, entry_count{0};
        std::uint32_t boundary_offset{0}, boundary_count{0};
        bool deferred{false};
    };
    std::vector<LinearTemporalRegion> linear_regions_;
    std::vector<SccId> linear_region_sccs_;
    std::vector<NodeId> linear_region_nodes_, linear_region_boundaries_;
    std::vector<LinearEntry> linear_region_entries_;
    std::vector<std::uint32_t> linear_region_of_scc_;
    std::vector<std::uint32_t> linear_region_node_slot_;
    std::vector<std::uint64_t> linear_region_dependency_offsets_;
    std::vector<std::uint32_t> linear_region_dependency_targets_;
    std::vector<std::uint64_t> linear_region_entry_row_offsets_;
    std::vector<std::uint32_t> linear_region_entry_row_indices_;
    std::vector<std::uint8_t> linear_region_baseline_stationary_;
    std::vector<std::uint32_t> causal_cone_generation_, causal_cone_queue_;
    std::vector<std::uint32_t> causal_cone_depth_;
    std::vector<std::uint32_t> active_region_node_slots_, active_region_entry_slots_;
    std::vector<NodeId> active_boundary_nodes_;
    std::vector<std::uint32_t> causal_cone_scc_generation_, active_region_sccs_;
    std::uint32_t current_causal_cone_generation_{1};
    std::uint32_t current_causal_scc_generation_{1};
    ReactionId active_perturbed_reaction_{kInvalidReaction};
    NodeId active_perturbed_node_{kInvalidNode};
    std::uint64_t causal_cone_cache_signature_{0};
    std::uint64_t causal_cone_cache_baseline_id_{0};
    std::uint64_t causal_cone_cache_context_id_{0};
    std::uint32_t causal_cone_cache_region_{std::numeric_limits<std::uint32_t>::max()};
    std::uint32_t causal_cone_cache_horizon_{0};
    bool causal_cone_cache_valid_{false};
    double region_construction_time_us_{0.0};
    std::vector<std::uint32_t> cached_causal_cone_nodes_, cached_causal_cone_entries_;
    std::vector<NodeId> cached_causal_boundary_nodes_;
    std::vector<double> region_accumulator_, region_initial_state_, region_segment_response_;
    std::vector<std::uint32_t> boundary_segment_start_;
    std::vector<double> dense_base_, dense_power_, dense_product_;
    std::vector<std::uint32_t> dense_free_index_;
    struct PropagatorCacheSlot {
        bool valid{false};
        std::uint64_t signature{0};
        SccId scc{kInvalidScc};
        std::uint32_t dimension{0};
        std::vector<double> matrix;
    };
    std::vector<PropagatorCacheSlot> propagator_cache_;
    struct PowerActionCacheSlot {
        bool valid{false};
        std::uint64_t signature{0};
        SccId scc{kInvalidScc};
        std::uint32_t horizon{0};
        std::uint32_t terms{0};
        std::vector<double> values;
    };
    std::vector<PowerActionCacheSlot> power_action_cache_;
    std::size_t next_propagator_slot_{0};
    std::size_t next_power_action_slot_{0};
    std::uint32_t current_dirty_generation_{1};
    std::uint32_t current_scc_visit_generation_{1};
    bool journaling_{false};
    std::uint64_t context_id_{0};
    std::uint64_t baseline_state_id_{0};

    struct CacheKey {
        std::uint64_t network_version{0};
        std::uint64_t context_id{0};
        std::uint64_t baseline_state_id{0};
        std::uint32_t perturbation_type{0};
        NodeId node{kInvalidNode};
        ReactionId reaction{kInvalidReaction};
        std::uint32_t strength_bits{0};
        std::uint32_t strategy{0};

        bool operator==(const CacheKey&) const = default;
    };

    struct CacheKeyHash {
        [[nodiscard]] std::size_t operator()(const CacheKey& key) const noexcept;
    };

    std::unordered_map<CacheKey, SimulationResult, CacheKeyHash> hot_cache_;

    [[nodiscard]] float combine_values(
        const NodeId* nodes,
        const float* weights,
        const float* K,
        const float* n,
        std::uint32_t count,
        RegulationOp op,
        float default_K,
        float default_n) const;
    [[nodiscard]] float reaction_signal(const Reaction& reaction) const;
    [[nodiscard]] float reaction_flux(const Reaction& reaction) const;
    [[nodiscard]] float evaluate_scc(SccId scc, ExecutionStats& stats);
    [[nodiscard]] float evaluate_scc_dirty(SccId scc, ExecutionStats& stats);
    [[nodiscard]] float evaluate_scc_temporal(SccId scc,
                                               std::uint32_t remaining_steps,
                                               std::uint32_t max_steps,
                                               ExecutionStats& stats);
    [[nodiscard]] bool can_evaluate_scc_temporal(SccId scc,
                                                  std::uint32_t remaining_steps,
                                                  std::uint32_t max_steps) const;
    [[nodiscard]] bool scalar_temporal_coefficients(SccId scc, double& p, double& q) const;
    void initialize_linear_workspace();
    void prepare_linear_coefficients(SccId scc);
    bool evaluate_dense_linear(SccId scc, std::uint32_t horizon, ExecutionStats& stats);
    bool evaluate_sparse_power_action(SccId scc, std::uint32_t horizon, ExecutionStats& stats);
    void evaluate_sparse_linear(SccId scc, ExecutionStats& stats);
    void execute_deferred_boundary_sccs(std::uint32_t horizon, ExecutionStats& stats);
    bool evaluate_sparse_power_action_region(std::uint32_t region, std::uint32_t horizon,
                                             ExecutionStats& stats,
                                             bool boundary_only = false,
                                             bool homogeneous_only = false,
                                             bool suppress_boundary = false,
                                             std::span<const std::uint32_t> active_nodes = {},
                                             std::span<const std::uint32_t> active_entries = {});
    void record_linear_execution(SccId scc, ExecutionStats& stats);
    [[nodiscard]] float commit_pending(ExecutionStats& stats, std::vector<float>* scc_delta = nullptr);
    void remember_before_change(NodeId node);
    void set_state(NodeId node, float value);
    void apply_perturbation(const Perturbation& perturbation);
    void rollback(const Perturbation& perturbation,
                  const std::optional<std::pair<ReactionId, float>>& reaction_change);
    [[nodiscard]] ExecutionStats run_full(
        std::uint32_t max_steps, float convergence_threshold,
        const BaselineStepObserver& observer = {});
    [[nodiscard]] ExecutionStats run_frontier(const std::vector<SccId>& seeds, std::uint32_t max_steps);
    void enqueue_unique(std::vector<SccId>& queue, SccId scc);
    void mark_scc_all_dirty(SccId scc);
    void mark_reaction_slot_dirty(SccId scc, std::uint32_t slot);
    void mark_node_dependents_dirty(NodeId node);
    void begin_dirty_generation();
    void record_scc_visit(SccId scc, ExecutionStats& stats);
    void record_scc_class(SccId scc, bool analytic, ExecutionStats& stats);
    void record_reaction_evaluation(SccId scc, ExecutionStats& stats);
    void collect_result_into(const ExecutionStats& stats,
                             double elapsed_us,
                             SimulationResult& result) const;
    [[nodiscard]] CacheKey make_cache_key(const Perturbation& perturbation,
                                          ExecutionStrategy strategy) const;
};

}  // namespace cellnet
