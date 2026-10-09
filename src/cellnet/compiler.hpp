#pragma once

#include "cellnet/network.hpp"

#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>

namespace cellnet {

class CompileError : public std::runtime_error {
public:
    explicit CompileError(const std::string& message) : std::runtime_error(message) {}
};

struct CompileMetrics {
    double total_ms{0.0};
    double reaction_lowering_ms{0.0};
    double dependency_csr_ms{0.0};
    double scc_condensation_ms{0.0};
    double execution_plan_ms{0.0};
    double temporal_plan_ms{0.0};
};

struct TemporalCompileConfig {
    float dt{0.05F};
    std::uint32_t horizon_steps{100U};
};

class Compiler {
public:
    [[nodiscard]] CompiledNetwork compile(const NetworkSpec& spec,
                                          std::uint64_t network_version = 1U,
                                          CompileMetrics* metrics = nullptr,
                                          TemporalCompileConfig temporal = {}) const;
};

}  // namespace cellnet
