#include "cellnet/binary.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <span>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cellnet {
namespace {

constexpr std::array<std::uint8_t, 8> kMagic = {'C', 'E', 'L', 'L', 'N', 'E', 'T', '1'};
constexpr std::uint32_t kFileVersion = 3U;
constexpr std::uint32_t kEndianMarker = 0x01020304U;
constexpr std::size_t kHeaderBytes = 32U;

class Encoder {
public:
    void u8(std::uint8_t value) { data_.push_back(value); }

    void u32(std::uint32_t value) {
        for (unsigned shift = 0U; shift < 32U; shift += 8U) {
            data_.push_back(static_cast<std::uint8_t>((value >> shift) & 0xFFU));
        }
    }

    void u64(std::uint64_t value) {
        for (unsigned shift = 0U; shift < 64U; shift += 8U) {
            data_.push_back(static_cast<std::uint8_t>((value >> shift) & 0xFFU));
        }
    }

    void f32(float value) { u32(std::bit_cast<std::uint32_t>(value)); }

    void text(const std::string& value) {
        u64(value.size());
        data_.insert(data_.end(), value.begin(), value.end());
    }

    void bytes(std::span<const std::uint8_t> values) {
        data_.insert(data_.end(), values.begin(), values.end());
    }

    [[nodiscard]] const std::vector<std::uint8_t>& data() const noexcept { return data_; }

private:
    std::vector<std::uint8_t> data_;
};

class Decoder {
public:
    explicit Decoder(std::span<const std::uint8_t> data) : data_(data) {}

    [[nodiscard]] std::uint8_t u8() {
        require(1U);
        return data_[position_++];
    }

    [[nodiscard]] std::uint32_t u32() {
        require(4U);
        std::uint32_t value = 0U;
        for (unsigned shift = 0U; shift < 32U; shift += 8U) {
            value |= static_cast<std::uint32_t>(data_[position_++]) << shift;
        }
        return value;
    }

    [[nodiscard]] std::uint64_t u64() {
        require(8U);
        std::uint64_t value = 0U;
        for (unsigned shift = 0U; shift < 64U; shift += 8U) {
            value |= static_cast<std::uint64_t>(data_[position_++]) << shift;
        }
        return value;
    }

    [[nodiscard]] float f32() { return std::bit_cast<float>(u32()); }

    [[nodiscard]] std::string text() {
        const std::uint64_t count = u64();
        if (count > remaining()) {
            throw BinaryFormatError("string extends past the binary payload");
        }
        const std::size_t size = static_cast<std::size_t>(count);
        const char* begin = reinterpret_cast<const char*>(data_.data() + position_);
        std::string value(begin, size);
        position_ += size;
        return value;
    }

    [[nodiscard]] std::span<const std::uint8_t> bytes(std::size_t count) {
        require(count);
        const auto result = data_.subspan(position_, count);
        position_ += count;
        return result;
    }

    [[nodiscard]] std::size_t remaining() const noexcept { return data_.size() - position_; }
    [[nodiscard]] bool empty() const noexcept { return remaining() == 0U; }

    [[nodiscard]] std::size_t count(std::size_t minimum_element_bytes) {
        const std::uint64_t raw_count = u64();
        if (raw_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw BinaryFormatError("vector count exceeds platform size capacity");
        }
        if (minimum_element_bytes > 0U &&
            raw_count > static_cast<std::uint64_t>(remaining() / minimum_element_bytes)) {
            throw BinaryFormatError("vector count extends past the binary payload");
        }
        return static_cast<std::size_t>(raw_count);
    }

private:
    void require(std::size_t count) const {
        if (count > remaining()) {
            throw BinaryFormatError("unexpected end of compiled network file");
        }
    }

    std::span<const std::uint8_t> data_;
    std::size_t position_{0U};
};

std::uint64_t checksum(std::span<const std::uint8_t> bytes) {
    std::uint64_t hash = 14695981039346656037ULL;
    for (const std::uint8_t byte : bytes) {
        hash ^= byte;
        hash *= 1099511628211ULL;
    }
    return hash;
}

template <typename T, typename Writer>
void encode_vector(Encoder& encoder, const std::vector<T>& values, Writer writer) {
    encoder.u64(values.size());
    for (const auto& value : values) writer(encoder, value);
}

template <typename T, typename Reader>
std::vector<T> decode_vector(Decoder& decoder, std::size_t minimum_bytes, Reader reader) {
    const std::size_t count = decoder.count(minimum_bytes);
    std::vector<T> values;
    values.reserve(count);
    for (std::size_t index = 0; index < count; ++index) values.push_back(reader(decoder));
    return values;
}

void encode_node(Encoder& encoder, const Node& node) {
    encoder.u32(node.id);
    encoder.u8(static_cast<std::uint8_t>(node.type));
    encoder.u8(static_cast<std::uint8_t>(node.compartment));
    encoder.f32(node.initial_value);
    encoder.u32(node.entity_id);
    encoder.u32(node.state_id);
}

Node decode_node(Decoder& decoder) {
    Node node;
    node.id = decoder.u32();
    const auto type = decoder.u8();
    const auto compartment = decoder.u8();
    if (type > static_cast<std::uint8_t>(NodeType::ExternalSignal) ||
        compartment > static_cast<std::uint8_t>(Compartment::Golgi)) {
        throw BinaryFormatError("invalid node enum in compiled network");
    }
    node.type = static_cast<NodeType>(type);
    node.compartment = static_cast<Compartment>(compartment);
    node.initial_value = decoder.f32();
    node.entity_id = decoder.u32();
    node.state_id = decoder.u32();
    return node;
}

void encode_reaction(Encoder& encoder, const Reaction& reaction) {
    encoder.u32(reaction.id);
    encoder.u8(static_cast<std::uint8_t>(reaction.primitive));
    encoder.u32(reaction.input_offset);
    encoder.u32(reaction.input_count);
    encoder.u32(reaction.output_offset);
    encoder.u32(reaction.output_count);
    encoder.u32(reaction.catalyst_offset);
    encoder.u32(reaction.catalyst_count);
    encoder.u32(reaction.pos_reg_offset);
    encoder.u32(reaction.pos_reg_count);
    encoder.u32(reaction.neg_reg_offset);
    encoder.u32(reaction.neg_reg_count);
    encoder.u8(static_cast<std::uint8_t>(reaction.input_op));
    encoder.u8(static_cast<std::uint8_t>(reaction.pos_reg_op));
    encoder.u8(static_cast<std::uint8_t>(reaction.neg_reg_op));
    encoder.f32(reaction.k);
    encoder.f32(reaction.Km);
    encoder.f32(reaction.hill_n);
    encoder.f32(reaction.confidence);
    encoder.f32(reaction.context_weight);
}

Reaction decode_reaction(Decoder& decoder) {
    Reaction reaction;
    reaction.id = decoder.u32();
    const auto primitive = decoder.u8();
    if (primitive > static_cast<std::uint8_t>(Primitive::Delay)) {
        throw BinaryFormatError("invalid reaction primitive in compiled network");
    }
    reaction.primitive = static_cast<Primitive>(primitive);
    reaction.input_offset = decoder.u32();
    reaction.input_count = decoder.u32();
    reaction.output_offset = decoder.u32();
    reaction.output_count = decoder.u32();
    reaction.catalyst_offset = decoder.u32();
    reaction.catalyst_count = decoder.u32();
    reaction.pos_reg_offset = decoder.u32();
    reaction.pos_reg_count = decoder.u32();
    reaction.neg_reg_offset = decoder.u32();
    reaction.neg_reg_count = decoder.u32();
    const auto input_op = decoder.u8();
    const auto pos_op = decoder.u8();
    const auto neg_op = decoder.u8();
    const auto maximum_op = static_cast<std::uint8_t>(RegulationOp::Xor);
    if (input_op > maximum_op || pos_op > maximum_op || neg_op > maximum_op) {
        throw BinaryFormatError("invalid regulation operator in compiled network");
    }
    reaction.input_op = static_cast<RegulationOp>(input_op);
    reaction.pos_reg_op = static_cast<RegulationOp>(pos_op);
    reaction.neg_reg_op = static_cast<RegulationOp>(neg_op);
    reaction.k = decoder.f32();
    reaction.Km = decoder.f32();
    reaction.hill_n = decoder.f32();
    reaction.confidence = decoder.f32();
    reaction.context_weight = decoder.f32();
    return reaction;
}

void encode_scc(Encoder& encoder, const SccBlock& block) {
    encoder.u32(block.id);
    encoder.u32(block.node_offset);
    encoder.u32(block.node_count);
    encoder.u32(block.reaction_offset);
    encoder.u32(block.reaction_count);
    encoder.u32(block.downstream_offset);
    encoder.u32(block.downstream_count);
    encoder.u8(block.cyclic ? 1U : 0U);
}

void encode_plan(Encoder& encoder, const SccExecutionPlan& plan) {
    encoder.u32(plan.modulate_offset);
    encoder.u32(plan.modulate_count);
    encoder.u32(plan.transfer_offset);
    encoder.u32(plan.transfer_count);
    encoder.u32(plan.produce_offset);
    encoder.u32(plan.produce_count);
    encoder.u32(plan.destroy_offset);
    encoder.u32(plan.destroy_count);
    encoder.u32(plan.generic_offset);
    encoder.u32(plan.generic_count);
}

SccExecutionPlan decode_plan(Decoder& decoder) {
    SccExecutionPlan plan;
    plan.modulate_offset = decoder.u32();
    plan.modulate_count = decoder.u32();
    plan.transfer_offset = decoder.u32();
    plan.transfer_count = decoder.u32();
    plan.produce_offset = decoder.u32();
    plan.produce_count = decoder.u32();
    plan.destroy_offset = decoder.u32();
    plan.destroy_count = decoder.u32();
    plan.generic_offset = decoder.u32();
    plan.generic_count = decoder.u32();
    return plan;
}

void encode_temporal_plan(Encoder& encoder, const SccTemporalPlan& plan) {
    encoder.u32(static_cast<std::uint32_t>(plan.classification));
    encoder.u32(plan.bounded);
    encoder.u32(plan.propagator);
    encoder.u32(plan.external_inputs);
    encoder.u32(plan.node_count);
    encoder.u32(plan.matrix_offset);
    encoder.u32(plan.matrix_count);
    encoder.u32(plan.integral_offset);
    encoder.u32(plan.integral_count);
    encoder.u32(plan.affine_offset);
    encoder.u32(plan.affine_count);
}

SccTemporalPlan decode_temporal_plan(Decoder& decoder) {
    SccTemporalPlan plan;
    const auto classification = decoder.u32();
    const auto bounded = decoder.u32();
    const auto propagator = decoder.u32();
    const auto external_inputs = decoder.u32();
    if (classification > static_cast<std::uint32_t>(SccTemporalClass::GenericNonlinear) ||
        bounded > 1U || propagator > 1U || external_inputs > 1U) {
        throw BinaryFormatError("invalid SCC temporal plan enum");
    }
    plan.classification = static_cast<SccTemporalClass>(classification);
    plan.bounded = static_cast<std::uint8_t>(bounded);
    plan.propagator = static_cast<std::uint8_t>(propagator);
    plan.external_inputs = static_cast<std::uint8_t>(external_inputs);
    plan.node_count = decoder.u32();
    plan.matrix_offset = decoder.u32();
    plan.matrix_count = decoder.u32();
    plan.integral_offset = decoder.u32();
    plan.integral_count = decoder.u32();
    plan.affine_offset = decoder.u32();
    plan.affine_count = decoder.u32();
    return plan;
}

SccBlock decode_scc(Decoder& decoder) {
    SccBlock block;
    block.id = decoder.u32();
    block.node_offset = decoder.u32();
    block.node_count = decoder.u32();
    block.reaction_offset = decoder.u32();
    block.reaction_count = decoder.u32();
    block.downstream_offset = decoder.u32();
    block.downstream_count = decoder.u32();
    const std::uint8_t cyclic = decoder.u8();
    if (cyclic > 1U) throw BinaryFormatError("invalid SCC cyclic flag");
    block.cyclic = cyclic != 0U;
    return block;
}

template <typename T>
void encode_unsigned_vector(Encoder& encoder, const std::vector<T>& values) {
    static_assert(std::is_unsigned_v<T>);
    encode_vector(encoder, values, [](Encoder& output, T value) {
        if constexpr (sizeof(T) == 4U) output.u32(value);
        else output.u64(value);
    });
}

template <typename T>
std::vector<T> decode_unsigned_vector(Decoder& decoder) {
    static_assert(std::is_unsigned_v<T>);
    return decode_vector<T>(decoder, sizeof(T), [](Decoder& input) {
        if constexpr (sizeof(T) == 4U) return static_cast<T>(input.u32());
        else return static_cast<T>(input.u64());
    });
}

void encode_float_vector(Encoder& encoder, const std::vector<float>& values) {
    encode_vector(encoder, values, [](Encoder& output, float value) { output.f32(value); });
}

std::vector<float> decode_float_vector(Decoder& decoder) {
    return decode_vector<float>(decoder, 4U, [](Decoder& input) { return input.f32(); });
}

void validate_offsets(const std::vector<std::uint64_t>& offsets, std::size_t owner_count,
                      std::size_t target_count, const char* label) {
    if (offsets.size() != owner_count + 1U || offsets.empty() || offsets.front() != 0U ||
        offsets.back() != target_count ||
        !std::is_sorted(offsets.begin(), offsets.end())) {
        throw BinaryFormatError(std::string("invalid ") + label + " offsets");
    }
}

void validate_network(const CompiledNetwork& network) {
    const std::size_t node_count = network.nodes.size();
    const std::size_t reaction_count = network.reactions.size();
    const std::size_t scc_count = network.sccs.size();
    if (network.node_names.size() != node_count || network.node_external_ids.size() != node_count ||
        network.node_to_scc.size() != node_count ||
        network.node_position_in_scc.size() != node_count ||
        network.reaction_external_ids.size() != reaction_count) {
        throw BinaryFormatError("compiled network metadata counts do not match");
    }
    if (network.reaction_inputs.size() != network.input_stoichiometry.size() ||
        network.reaction_outputs.size() != network.output_stoichiometry.size() ||
        network.positive_regulators.size() != network.positive_weights.size() ||
        network.positive_regulators.size() != network.positive_K.size() ||
        network.positive_regulators.size() != network.positive_n.size() ||
        network.negative_regulators.size() != network.negative_weights.size() ||
        network.negative_regulators.size() != network.negative_K.size() ||
        network.negative_regulators.size() != network.negative_n.size()) {
        throw BinaryFormatError("compiled network parallel arrays do not match");
    }
    if (network.execution_plans.size() != scc_count ||
        network.plan_modulate_reactions.size() != network.plan_modulate_slots.size() ||
        network.plan_modulate_reactions.size() != network.plan_modulate_output_locals.size() ||
        network.plan_modulate_reactions.size() != network.plan_modulate_output_stoich.size() ||
        network.plan_transfer_reactions.size() != network.plan_transfer_slots.size() ||
        network.plan_produce_reactions.size() != network.plan_produce_slots.size() ||
        network.plan_destroy_reactions.size() != network.plan_destroy_slots.size() ||
        network.plan_generic_reactions.size() != network.plan_generic_slots.size()) {
        throw BinaryFormatError("compiled execution plan arrays do not match");
    }
    if (network.temporal_plans.size() != scc_count ||
        !std::isfinite(network.temporal_compile_dt) || network.temporal_compile_dt <= 0.0F ||
        network.temporal_compile_steps == 0U) {
        throw BinaryFormatError("compiled temporal plan metadata does not match");
    }
    validate_offsets(network.scc_node_reaction_offsets, network.scc_nodes.size(),
                     network.scc_node_reaction_slots.size(), "SCC node reaction");
    validate_offsets(network.dependency_offsets, node_count, network.dependency_targets.size(),
                     "dependency");
    validate_offsets(network.reaction_scc_offsets, reaction_count, network.reaction_sccs.size(),
                     "reaction SCC");

    auto validate_temporal_range = [](std::uint32_t offset,
                                      std::uint32_t count,
                                      std::size_t size,
                                      const char* label) {
        if (static_cast<std::uint64_t>(offset) + count > size) {
            throw BinaryFormatError(std::string("invalid ") + label + " range");
        }
    };
    for (std::size_t index = 0; index < node_count; ++index) {
        if (network.nodes[index].id != index || network.node_to_scc[index] >= scc_count ||
            !std::isfinite(network.nodes[index].initial_value) ||
            network.nodes[index].initial_value < 0.0F ||
            network.nodes[index].initial_value > 1.0F) {
            throw BinaryFormatError("invalid node id or SCC mapping");
        }
    }
    auto check_nodes = [node_count](const std::vector<NodeId>& nodes, const char* label) {
        if (std::any_of(nodes.begin(), nodes.end(), [node_count](NodeId node) {
                return node >= node_count;
            })) {
            throw BinaryFormatError(std::string("out-of-range node in ") + label);
        }
    };
    check_nodes(network.reaction_inputs, "reaction inputs");
    check_nodes(network.reaction_outputs, "reaction outputs");
    check_nodes(network.catalysts, "catalysts");
    check_nodes(network.positive_regulators, "positive regulators");
    check_nodes(network.negative_regulators, "negative regulators");
    check_nodes(network.dependency_targets, "dependency targets");
    check_nodes(network.scc_nodes, "SCC nodes");
    check_nodes(network.phenotype_nodes, "phenotype nodes");

    for (std::size_t index = 0; index < reaction_count; ++index) {
        const auto& reaction = network.reactions[index];
        if (reaction.id != index ||
            static_cast<std::uint64_t>(reaction.input_offset) + reaction.input_count >
                network.reaction_inputs.size() ||
            static_cast<std::uint64_t>(reaction.output_offset) + reaction.output_count >
                network.reaction_outputs.size() ||
            static_cast<std::uint64_t>(reaction.catalyst_offset) + reaction.catalyst_count >
                network.catalysts.size() ||
            static_cast<std::uint64_t>(reaction.pos_reg_offset) + reaction.pos_reg_count >
                network.positive_regulators.size() ||
            static_cast<std::uint64_t>(reaction.neg_reg_offset) + reaction.neg_reg_count >
                network.negative_regulators.size() ||
            !std::isfinite(reaction.k) || !std::isfinite(reaction.Km) ||
            !std::isfinite(reaction.hill_n) || !std::isfinite(reaction.confidence) ||
            !std::isfinite(reaction.context_weight) || reaction.k < 0.0F ||
            reaction.Km < 0.0F || reaction.hill_n <= 0.0F || reaction.confidence < 0.0F ||
            reaction.confidence > 1.0F || reaction.context_weight < 0.0F ||
            reaction.context_weight > 1.0F) {
            throw BinaryFormatError("invalid reaction offset or id");
        }
    }

    auto valid_positive_floats = [](const std::vector<float>& values) {
        return std::all_of(values.begin(), values.end(), [](float value) {
            return std::isfinite(value) && value > 0.0F;
        });
    };
    auto valid_finite_floats = [](const std::vector<float>& values) {
        return std::all_of(values.begin(), values.end(), [](float value) {
            return std::isfinite(value);
        });
    };
    if (!valid_positive_floats(network.input_stoichiometry) ||
        !valid_positive_floats(network.output_stoichiometry) ||
        !valid_finite_floats(network.positive_weights) ||
        !valid_finite_floats(network.negative_weights) ||
        !valid_positive_floats(network.positive_n) ||
        !valid_positive_floats(network.negative_n) ||
        !std::all_of(network.positive_K.begin(), network.positive_K.end(), [](float value) {
            return std::isfinite(value) && value >= 0.0F;
        }) ||
        !std::all_of(network.negative_K.begin(), network.negative_K.end(), [](float value) {
            return std::isfinite(value) && value >= 0.0F;
        })) {
        throw BinaryFormatError("compiled network contains invalid numeric arrays");
    }

    if (network.topological_sccs.size() != scc_count) {
        throw BinaryFormatError("topological SCC count does not match");
    }
    std::vector<std::uint8_t> seen_scc(scc_count, 0U);
    for (const SccId scc : network.topological_sccs) {
        if (scc >= scc_count || seen_scc[scc] != 0U) {
            throw BinaryFormatError("invalid topological SCC order");
        }
        seen_scc[scc] = 1U;
    }
    std::vector<std::uint8_t> seen_nodes(node_count, 0U);
    for (std::size_t index = 0; index < scc_count; ++index) {
        const auto& block = network.sccs[index];
        if (block.id != index ||
            static_cast<std::uint64_t>(block.node_offset) + block.node_count >
                network.scc_nodes.size() ||
            static_cast<std::uint64_t>(block.reaction_offset) + block.reaction_count >
                network.scc_reactions.size() ||
            static_cast<std::uint64_t>(block.downstream_offset) + block.downstream_count >
                network.downstream_sccs.size()) {
            throw BinaryFormatError("invalid SCC block offset or id");
        }
        for (std::uint32_t local = 0; local < block.node_count; ++local) {
            const NodeId node = network.scc_nodes[block.node_offset + local];
            if (network.node_to_scc[node] != block.id ||
                network.node_position_in_scc[node] != local || seen_nodes[node] != 0U) {
                throw BinaryFormatError("SCC node mapping is inconsistent");
            }
            seen_nodes[node] = 1U;
        }
        const auto& plan = network.execution_plans[index];
        if (static_cast<std::uint64_t>(plan.modulate_offset) + plan.modulate_count >
                network.plan_modulate_reactions.size() ||
            static_cast<std::uint64_t>(plan.transfer_offset) + plan.transfer_count >
                network.plan_transfer_reactions.size() ||
            static_cast<std::uint64_t>(plan.produce_offset) + plan.produce_count >
                network.plan_produce_reactions.size() ||
            static_cast<std::uint64_t>(plan.destroy_offset) + plan.destroy_count >
                network.plan_destroy_reactions.size() ||
            static_cast<std::uint64_t>(plan.generic_offset) + plan.generic_count >
                network.plan_generic_reactions.size()) {
            throw BinaryFormatError("invalid SCC execution plan range");
        }
        const auto validate_plan_slots = [&](std::uint32_t offset,
                                              std::uint32_t count,
                                              const std::vector<ReactionId>& reactions,
                                              const std::vector<std::uint32_t>& slots,
                                              const char* label) {
            for (std::uint32_t local = 0; local < count; ++local) {
                const std::size_t plan_index = static_cast<std::size_t>(offset) + local;
                const std::uint32_t slot = slots[plan_index];
                if (slot < block.reaction_offset ||
                    slot >= block.reaction_offset + block.reaction_count ||
                    reactions[plan_index] != network.scc_reactions[slot]) {
                    throw BinaryFormatError(std::string("invalid ") + label + " SCC slot");
                }
            }
        };
        validate_plan_slots(plan.modulate_offset, plan.modulate_count,
                            network.plan_modulate_reactions, network.plan_modulate_slots,
                            "modulate plan");
        validate_plan_slots(plan.transfer_offset, plan.transfer_count,
                            network.plan_transfer_reactions, network.plan_transfer_slots,
                            "transfer plan");
        validate_plan_slots(plan.produce_offset, plan.produce_count,
                            network.plan_produce_reactions, network.plan_produce_slots,
                            "produce plan");
        validate_plan_slots(plan.destroy_offset, plan.destroy_count,
                            network.plan_destroy_reactions, network.plan_destroy_slots,
                            "destroy plan");
        validate_plan_slots(plan.generic_offset, plan.generic_count,
                            network.plan_generic_reactions, network.plan_generic_slots,
                            "generic plan");
        for (std::uint32_t local = 0; local < plan.modulate_count; ++local) {
            const std::size_t plan_index = static_cast<std::size_t>(plan.modulate_offset) + local;
            if (network.plan_modulate_output_locals[plan_index] >= block.node_count ||
                !std::isfinite(network.plan_modulate_output_stoich[plan_index]) ||
                network.plan_modulate_output_stoich[plan_index] <= 0.0F) {
                throw BinaryFormatError("invalid modulate plan output");
            }
        }
        const auto& temporal = network.temporal_plans[index];
        if (static_cast<std::uint8_t>(temporal.classification) >
                static_cast<std::uint8_t>(SccTemporalClass::GenericNonlinear) ||
            temporal.bounded > 1U || temporal.propagator > 1U || temporal.external_inputs > 1U ||
            temporal.node_count != block.node_count ||
            (temporal.propagator == 0U &&
             (temporal.matrix_count != 0U || temporal.integral_count != 0U ||
              temporal.affine_count != 0U)) ||
            (temporal.propagator != 0U &&
             (temporal.matrix_count != static_cast<std::uint64_t>(block.node_count) * block.node_count ||
              temporal.integral_count != static_cast<std::uint64_t>(block.node_count) * block.node_count ||
              temporal.affine_count != block.node_count))) {
            throw BinaryFormatError("invalid SCC temporal plan dimensions");
        }
        validate_temporal_range(temporal.matrix_offset, temporal.matrix_count,
                                network.temporal_propagators.size(), "temporal propagator");
        validate_temporal_range(temporal.integral_offset, temporal.integral_count,
                                network.temporal_integrals.size(), "temporal integral");
        validate_temporal_range(temporal.affine_offset, temporal.affine_count,
                                network.temporal_affine.size(), "temporal affine");
    }
    if (!std::all_of(network.temporal_propagators.begin(), network.temporal_propagators.end(),
                     [](float value) { return std::isfinite(value); }) ||
        !std::all_of(network.temporal_integrals.begin(), network.temporal_integrals.end(),
                     [](float value) { return std::isfinite(value); }) ||
        !std::all_of(network.temporal_affine.begin(), network.temporal_affine.end(),
                     [](float value) { return std::isfinite(value); })) {
        throw BinaryFormatError("compiled temporal plan contains a non-finite value");
    }
    if (network.scc_nodes.size() != node_count) {
        throw BinaryFormatError("SCC node list does not cover all nodes");
    }
    if (std::any_of(network.scc_reactions.begin(), network.scc_reactions.end(),
                    [reaction_count](ReactionId reaction) { return reaction >= reaction_count; }) ||
        std::any_of(network.downstream_sccs.begin(), network.downstream_sccs.end(),
                    [scc_count](SccId scc) { return scc >= scc_count; }) ||
        std::any_of(network.reaction_sccs.begin(), network.reaction_sccs.end(),
                    [scc_count](SccId scc) { return scc >= scc_count; })) {
        throw BinaryFormatError("compiled network contains out-of-range SCC/reaction ids");
    }
    auto check_plan_slots = [total_slots = network.scc_reactions.size()](
                                const std::vector<std::uint32_t>& slots) {
        return std::all_of(slots.begin(), slots.end(), [total_slots](std::uint32_t slot) {
            return slot < total_slots;
        });
    };
    if (!check_plan_slots(network.plan_modulate_slots) ||
        !check_plan_slots(network.plan_transfer_slots) ||
        !check_plan_slots(network.plan_produce_slots) ||
        !check_plan_slots(network.plan_destroy_slots) ||
        !check_plan_slots(network.plan_generic_slots) ||
        std::any_of(network.scc_node_reaction_slots.begin(), network.scc_node_reaction_slots.end(),
                    [total_slots = network.scc_reactions.size()](std::uint32_t slot) {
                        return slot >= total_slots;
                    })) {
        throw BinaryFormatError("compiled execution plan contains an invalid slot");
    }
}

}  // namespace

void save_compiled_network(const CompiledNetwork& network, const std::filesystem::path& path) {
    validate_network(network);
    Encoder payload;
    payload.u64(network.network_version);
    encode_vector(payload, network.nodes, encode_node);
    encode_vector(payload, network.reactions, encode_reaction);
    encode_unsigned_vector(payload, network.reaction_inputs);
    encode_float_vector(payload, network.input_stoichiometry);
    encode_unsigned_vector(payload, network.reaction_outputs);
    encode_float_vector(payload, network.output_stoichiometry);
    encode_unsigned_vector(payload, network.catalysts);
    encode_unsigned_vector(payload, network.positive_regulators);
    encode_float_vector(payload, network.positive_weights);
    encode_float_vector(payload, network.positive_K);
    encode_float_vector(payload, network.positive_n);
    encode_unsigned_vector(payload, network.negative_regulators);
    encode_float_vector(payload, network.negative_weights);
    encode_float_vector(payload, network.negative_K);
    encode_float_vector(payload, network.negative_n);
    encode_unsigned_vector(payload, network.dependency_offsets);
    encode_unsigned_vector(payload, network.dependency_targets);
    encode_unsigned_vector(payload, network.node_to_scc);
    encode_unsigned_vector(payload, network.node_position_in_scc);
    encode_vector(payload, network.sccs, encode_scc);
    encode_unsigned_vector(payload, network.scc_nodes);
    encode_unsigned_vector(payload, network.scc_reactions);
    encode_unsigned_vector(payload, network.downstream_sccs);
    encode_unsigned_vector(payload, network.topological_sccs);
    encode_unsigned_vector(payload, network.reaction_scc_offsets);
    encode_unsigned_vector(payload, network.reaction_sccs);
    encode_vector(payload, network.execution_plans, encode_plan);
    encode_unsigned_vector(payload, network.plan_modulate_reactions);
    encode_unsigned_vector(payload, network.plan_modulate_slots);
    encode_unsigned_vector(payload, network.plan_modulate_output_locals);
    encode_float_vector(payload, network.plan_modulate_output_stoich);
    encode_unsigned_vector(payload, network.plan_transfer_reactions);
    encode_unsigned_vector(payload, network.plan_transfer_slots);
    encode_unsigned_vector(payload, network.plan_produce_reactions);
    encode_unsigned_vector(payload, network.plan_produce_slots);
    encode_unsigned_vector(payload, network.plan_destroy_reactions);
    encode_unsigned_vector(payload, network.plan_destroy_slots);
    encode_unsigned_vector(payload, network.plan_generic_reactions);
    encode_unsigned_vector(payload, network.plan_generic_slots);
    encode_unsigned_vector(payload, network.scc_node_reaction_offsets);
    encode_unsigned_vector(payload, network.scc_node_reaction_slots);
    encode_vector(payload, network.temporal_plans, encode_temporal_plan);
    encode_float_vector(payload, network.temporal_propagators);
    encode_float_vector(payload, network.temporal_integrals);
    encode_float_vector(payload, network.temporal_affine);
    payload.f32(network.temporal_compile_dt);
    payload.u32(network.temporal_compile_steps);
    encode_vector(payload, network.node_names,
                  [](Encoder& output, const std::string& value) { output.text(value); });
    encode_unsigned_vector(payload, network.node_external_ids);
    encode_unsigned_vector(payload, network.reaction_external_ids);
    encode_unsigned_vector(payload, network.phenotype_nodes);

    Encoder header;
    header.bytes(kMagic);
    header.u32(kFileVersion);
    header.u32(kEndianMarker);
    header.u64(payload.data().size());
    header.u64(checksum(payload.data()));

    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) throw BinaryFormatError("cannot open binary output: " + path.string());
    output.write(reinterpret_cast<const char*>(header.data().data()),
                 static_cast<std::streamsize>(header.data().size()));
    output.write(reinterpret_cast<const char*>(payload.data().data()),
                 static_cast<std::streamsize>(payload.data().size()));
    if (!output) throw BinaryFormatError("failed while writing binary output: " + path.string());
}

CompiledNetwork load_compiled_network(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary | std::ios::ate);
    if (!input) throw BinaryFormatError("cannot open compiled network: " + path.string());
    const std::streamoff end = input.tellg();
    if (end < 0 || static_cast<std::uint64_t>(end) >
                       static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw BinaryFormatError("compiled network file is too large for this platform");
    }
    const std::size_t file_size = static_cast<std::size_t>(end);
    if (file_size < kHeaderBytes) throw BinaryFormatError("compiled network header is truncated");
    std::vector<std::uint8_t> bytes(file_size);
    input.seekg(0, std::ios::beg);
    input.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!input) throw BinaryFormatError("failed while reading compiled network: " + path.string());

    Decoder file(bytes);
    const auto magic = file.bytes(kMagic.size());
    if (!std::equal(magic.begin(), magic.end(), kMagic.begin())) {
        throw BinaryFormatError("compiled network magic does not match");
    }
    if (file.u32() != kFileVersion) throw BinaryFormatError("unsupported compiled network version");
    if (file.u32() != kEndianMarker) throw BinaryFormatError("compiled network endian marker is invalid");
    const std::uint64_t payload_size = file.u64();
    const std::uint64_t expected_checksum = file.u64();
    if (payload_size != file.remaining()) {
        throw BinaryFormatError("compiled network payload size does not match file size");
    }
    const auto payload_bytes = file.bytes(static_cast<std::size_t>(payload_size));
    if (checksum(payload_bytes) != expected_checksum) {
        throw BinaryFormatError("compiled network checksum mismatch");
    }

    Decoder payload(payload_bytes);
    CompiledNetwork network;
    network.network_version = payload.u64();
    network.nodes = decode_vector<Node>(payload, 18U, decode_node);
    network.reactions = decode_vector<Reaction>(payload, 64U, decode_reaction);
    network.reaction_inputs = decode_unsigned_vector<NodeId>(payload);
    network.input_stoichiometry = decode_float_vector(payload);
    network.reaction_outputs = decode_unsigned_vector<NodeId>(payload);
    network.output_stoichiometry = decode_float_vector(payload);
    network.catalysts = decode_unsigned_vector<NodeId>(payload);
    network.positive_regulators = decode_unsigned_vector<NodeId>(payload);
    network.positive_weights = decode_float_vector(payload);
    network.positive_K = decode_float_vector(payload);
    network.positive_n = decode_float_vector(payload);
    network.negative_regulators = decode_unsigned_vector<NodeId>(payload);
    network.negative_weights = decode_float_vector(payload);
    network.negative_K = decode_float_vector(payload);
    network.negative_n = decode_float_vector(payload);
    network.dependency_offsets = decode_unsigned_vector<std::uint64_t>(payload);
    network.dependency_targets = decode_unsigned_vector<NodeId>(payload);
    network.node_to_scc = decode_unsigned_vector<SccId>(payload);
    network.node_position_in_scc = decode_unsigned_vector<std::uint32_t>(payload);
    network.sccs = decode_vector<SccBlock>(payload, 29U, decode_scc);
    network.scc_nodes = decode_unsigned_vector<NodeId>(payload);
    network.scc_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.downstream_sccs = decode_unsigned_vector<SccId>(payload);
    network.topological_sccs = decode_unsigned_vector<SccId>(payload);
    network.reaction_scc_offsets = decode_unsigned_vector<std::uint64_t>(payload);
    network.reaction_sccs = decode_unsigned_vector<SccId>(payload);
    network.execution_plans = decode_vector<SccExecutionPlan>(payload, 40U, decode_plan);
    network.plan_modulate_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.plan_modulate_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.plan_modulate_output_locals = decode_unsigned_vector<std::uint32_t>(payload);
    network.plan_modulate_output_stoich = decode_float_vector(payload);
    network.plan_transfer_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.plan_transfer_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.plan_produce_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.plan_produce_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.plan_destroy_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.plan_destroy_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.plan_generic_reactions = decode_unsigned_vector<ReactionId>(payload);
    network.plan_generic_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.scc_node_reaction_offsets = decode_unsigned_vector<std::uint64_t>(payload);
    network.scc_node_reaction_slots = decode_unsigned_vector<std::uint32_t>(payload);
    network.temporal_plans = decode_vector<SccTemporalPlan>(payload, 44U, decode_temporal_plan);
    network.temporal_propagators = decode_float_vector(payload);
    network.temporal_integrals = decode_float_vector(payload);
    network.temporal_affine = decode_float_vector(payload);
    network.temporal_compile_dt = payload.f32();
    network.temporal_compile_steps = payload.u32();
    network.node_names = decode_vector<std::string>(
        payload, 8U, [](Decoder& input_decoder) { return input_decoder.text(); });
    network.node_external_ids = decode_unsigned_vector<std::uint64_t>(payload);
    network.reaction_external_ids = decode_unsigned_vector<std::uint64_t>(payload);
    network.phenotype_nodes = decode_unsigned_vector<NodeId>(payload);
    if (!payload.empty()) throw BinaryFormatError("compiled network has trailing payload bytes");
    validate_network(network);
    compile_linear_execution(network);
    return network;
}

}  // namespace cellnet
