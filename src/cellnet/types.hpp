#pragma once

#include <cstdint>
#include <limits>

namespace cellnet {

using NodeId = std::uint32_t;
using ReactionId = std::uint32_t;
using SccId = std::uint32_t;

inline constexpr NodeId kInvalidNode = std::numeric_limits<NodeId>::max();
inline constexpr ReactionId kInvalidReaction = std::numeric_limits<ReactionId>::max();
inline constexpr SccId kInvalidScc = std::numeric_limits<SccId>::max();

enum class NodeType : std::uint8_t {
    Gene,
    RNA,
    Protein,
    ProteinState,
    Complex,
    Metabolite,
    Phenotype,
    ExternalSignal,
};

enum class Compartment : std::uint8_t {
    Unknown,
    Extracellular,
    Membrane,
    Cytosol,
    Nucleus,
    Mitochondria,
    ER,
    Golgi,
};

enum class Primitive : std::uint8_t {
    Transfer,
    Produce,
    Destroy,
    Bind,
    Unbind,
    Modulate,
    Move,
    Delay,
};

enum class RegulationOp : std::uint8_t {
    None,
    Sum,
    Product,
    Min,
    Max,
    HillPositive,
    HillNegative,
    ContinuousAnd,
    ContinuousOr,
    Threshold,
    Xor,
};

enum class SimulationMode : std::uint8_t {
    Continuous,
    Boolean,
};

enum class ExecutionStrategy : std::uint8_t {
    Full,
    Frontier,
};

enum class PerturbationType : std::uint8_t {
    GeneKnockout,
    GeneKnockdown,
    ProteinDegradation,
    ActivityInhibition,
    ActivityActivation,
    ReactionInhibition,
    ReactionActivation,
    EdgeBlock,
    Clamp,
};

}  // namespace cellnet
