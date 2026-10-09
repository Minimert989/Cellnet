#pragma once

#include "cellnet/types.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace cellnet {

struct Node {
    NodeId id{kInvalidNode};
    NodeType type{NodeType::Protein};
    Compartment compartment{Compartment::Unknown};
    float initial_value{0.0F};
    std::uint32_t entity_id{0};
    std::uint32_t state_id{0};
};

struct Reaction {
    ReactionId id{kInvalidReaction};
    Primitive primitive{Primitive::Modulate};

    std::uint32_t input_offset{0};
    std::uint32_t input_count{0};
    std::uint32_t output_offset{0};
    std::uint32_t output_count{0};
    std::uint32_t catalyst_offset{0};
    std::uint32_t catalyst_count{0};
    std::uint32_t pos_reg_offset{0};
    std::uint32_t pos_reg_count{0};
    std::uint32_t neg_reg_offset{0};
    std::uint32_t neg_reg_count{0};

    RegulationOp input_op{RegulationOp::Product};
    RegulationOp pos_reg_op{RegulationOp::Product};
    RegulationOp neg_reg_op{RegulationOp::Product};

    float k{1.0F};
    float Km{0.5F};
    float hill_n{1.0F};
    float confidence{1.0F};
    float context_weight{1.0F};
};

struct SccBlock {
    SccId id{kInvalidScc};
    std::uint32_t node_offset{0};
    std::uint32_t node_count{0};
    std::uint32_t reaction_offset{0};
    std::uint32_t reaction_count{0};
    std::uint32_t downstream_offset{0};
    std::uint32_t downstream_count{0};
    bool cyclic{false};
};

// Ranges into compiler-generated SCC-local execution blocks. The optimized frontier executor
// iterates these ranges directly; it does not dispatch on Primitive for common reactions.
struct SccExecutionPlan {
    std::uint32_t modulate_offset{0};
    std::uint32_t modulate_count{0};
    std::uint32_t transfer_offset{0};
    std::uint32_t transfer_count{0};
    std::uint32_t produce_offset{0};
    std::uint32_t produce_count{0};
    std::uint32_t destroy_offset{0};
    std::uint32_t destroy_count{0};
    std::uint32_t generic_offset{0};
    std::uint32_t generic_count{0};
};

enum class SccTemporalClass : std::uint8_t {
    Static,
    LinearHomogeneous,
    LinearAffine,
    SimpleNonlinear,
    GenericNonlinear,
};

enum class TemporalBackend : std::uint8_t {
    StaticDirect, ScalarClosedForm, PrecomputedDensePropagator,
    RuntimeDensePropagator, SparseLinear, BoundaryDrivenLinear, NumericalNonlinear
};
enum class LinearCoefficientKind : std::uint8_t { Fixed, Signal, Capped };
struct LinearEntry {
    NodeId row{kInvalidNode}, column{kInvalidNode}; // invalid column denotes constant
    ReactionId reaction{kInvalidReaction};
    float gain{0}, stoichiometry{1}, input_stoichiometry{1};
    LinearCoefficientKind kind{LinearCoefficientKind::Fixed};
};
struct LinearExecutionPlan {
    TemporalBackend backend{TemporalBackend::NumericalNonlinear};
    std::uint32_t entry_offset{0}, entry_count{0};
    std::uint32_t boundary_offset{0}, boundary_count{0};
    double density{0};
};

// Compiler-generated finite-horizon temporal plan. Matrix and integral ranges are flattened
// row-major arrays; for the default benchmark horizon they encode the exact explicit-Euler
// propagator x(T) = M x(0) + S b used by the reference equations.
struct SccTemporalPlan {
    SccTemporalClass classification{SccTemporalClass::GenericNonlinear};
    std::uint8_t bounded{0};
    std::uint8_t propagator{0};
    std::uint8_t external_inputs{0};
    std::uint8_t reserved{0};
    std::uint32_t node_count{0};
    std::uint32_t matrix_offset{0};
    std::uint32_t matrix_count{0};
    std::uint32_t integral_offset{0};
    std::uint32_t integral_count{0};
    std::uint32_t affine_offset{0};
    std::uint32_t affine_count{0};
};

struct NodeSpec {
    std::uint64_t external_id{0};
    std::string name;
    NodeType type{NodeType::Protein};
    Compartment compartment{Compartment::Unknown};
    float initial_value{0.0F};
    std::string entity;
    std::string state;
};

struct WeightedNodeSpec {
    std::uint64_t node_external_id{0};
    float coefficient{1.0F};
};

struct RegulatorSpec {
    std::uint64_t node_external_id{0};
    float weight{1.0F};
    float K{0.5F};
    float n{1.0F};
};

struct ReactionSpec {
    std::uint64_t external_id{0};
    Primitive primitive{Primitive::Modulate};
    RegulationOp input_op{RegulationOp::Product};
    RegulationOp pos_reg_op{RegulationOp::Product};
    RegulationOp neg_reg_op{RegulationOp::Product};
    float k{1.0F};
    float Km{0.5F};
    float hill_n{1.0F};
    float confidence{1.0F};
    float context_weight{1.0F};
    std::vector<WeightedNodeSpec> inputs;
    std::vector<WeightedNodeSpec> outputs;
    std::vector<std::uint64_t> catalysts;
    std::vector<RegulatorSpec> positive_regulators;
    std::vector<RegulatorSpec> negative_regulators;
};

struct NetworkSpec {
    std::vector<NodeSpec> nodes;
    std::vector<ReactionSpec> reactions;
};

struct CompiledNetwork {
    std::vector<Node> nodes;
    std::vector<Reaction> reactions;

    std::vector<NodeId> reaction_inputs;
    std::vector<float> input_stoichiometry;
    std::vector<NodeId> reaction_outputs;
    std::vector<float> output_stoichiometry;
    std::vector<NodeId> catalysts;
    std::vector<NodeId> positive_regulators;
    std::vector<float> positive_weights;
    std::vector<float> positive_K;
    std::vector<float> positive_n;
    std::vector<NodeId> negative_regulators;
    std::vector<float> negative_weights;
    std::vector<float> negative_K;
    std::vector<float> negative_n;

    // Node dependency graph in CSR form.
    std::vector<std::uint64_t> dependency_offsets;
    std::vector<NodeId> dependency_targets;

    std::vector<SccId> node_to_scc;
    std::vector<std::uint32_t> node_position_in_scc;
    std::vector<SccBlock> sccs;
    std::vector<NodeId> scc_nodes;
    std::vector<ReactionId> scc_reactions;
    std::vector<SccId> downstream_sccs;
    std::vector<SccId> topological_sccs;
    std::vector<std::uint64_t> reaction_scc_offsets;
    std::vector<SccId> reaction_sccs;

    // SCC-local execution plan blocks. Each slot is an index into the owning SCC's
    // scc_reactions range and is also used by the dirty-reaction generation array.
    std::vector<SccExecutionPlan> execution_plans;
    std::vector<ReactionId> plan_modulate_reactions;
    std::vector<std::uint32_t> plan_modulate_slots;
    std::vector<std::uint32_t> plan_modulate_output_locals;
    std::vector<float> plan_modulate_output_stoich;
    std::vector<ReactionId> plan_transfer_reactions;
    std::vector<std::uint32_t> plan_transfer_slots;
    std::vector<ReactionId> plan_produce_reactions;
    std::vector<std::uint32_t> plan_produce_slots;
    std::vector<ReactionId> plan_destroy_reactions;
    std::vector<std::uint32_t> plan_destroy_slots;
    std::vector<ReactionId> plan_generic_reactions;
    std::vector<std::uint32_t> plan_generic_slots;

    // For each flattened node in scc_nodes, the SCC-local reaction slots that depend on it.
    std::vector<std::uint64_t> scc_node_reaction_offsets;
    std::vector<std::uint32_t> scc_node_reaction_slots;

    // Temporal solver compilation metadata and finite-horizon propagators.
    std::vector<SccTemporalPlan> temporal_plans;
    std::vector<float> temporal_propagators;
    std::vector<float> temporal_integrals;
    std::vector<float> temporal_affine;
    float temporal_compile_dt{0.05F};
    std::uint32_t temporal_compile_steps{100};
    std::vector<LinearExecutionPlan> linear_plans;
    std::vector<LinearEntry> linear_entries;
    std::vector<NodeId> linear_boundaries;

    std::vector<std::string> node_names;
    std::vector<std::uint64_t> node_external_ids;
    std::vector<std::uint64_t> reaction_external_ids;
    std::vector<NodeId> phenotype_nodes;

    std::uint64_t network_version{1};

    [[nodiscard]] std::size_t estimated_bytes() const noexcept;
};

void compile_linear_execution(CompiledNetwork& network);

}  // namespace cellnet
