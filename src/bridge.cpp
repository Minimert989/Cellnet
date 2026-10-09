#include <Rcpp.h>

#include "cellnet/compiler.hpp"
#include "cellnet/engine.hpp"

#include <algorithm>
#include <atomic>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct CellNetHandle {
    cellnet::Engine engine;
    std::vector<std::string> names;
    std::unordered_map<std::string, cellnet::NodeId> ids;

    CellNetHandle(cellnet::CompiledNetwork network, cellnet::SolverConfig config,
                  std::vector<std::string> node_names)
        : engine(std::move(network), config), names(std::move(node_names)) {
        for (std::size_t i = 0; i < names.size(); ++i) {
            ids.emplace(names[i], static_cast<cellnet::NodeId>(i));
        }
        engine.initialize_baseline();
    }
};

CellNetHandle& checked_handle(SEXP pointer) {
    Rcpp::XPtr<CellNetHandle> handle(pointer);
    if (handle.get() == nullptr) Rcpp::stop("CellNet pointer is null");
    return *handle;
}

Rcpp::DataFrame materialize(const CellNetHandle& handle,
                            const cellnet::SimulationResult& result) {
    const R_xlen_t n = static_cast<R_xlen_t>(result.molecular_changes.size());
    Rcpp::CharacterVector node(n);
    Rcpp::NumericVector baseline(n), perturbed(n), delta(n);
    Rcpp::LogicalVector affected(n, true);
    Rcpp::IntegerVector distance(n, NA_INTEGER);
    for (R_xlen_t i = 0; i < n; ++i) {
        const auto& change = result.molecular_changes[static_cast<std::size_t>(i)];
        node[i] = handle.names.at(change.id);
        baseline[i] = change.baseline;
        perturbed[i] = change.value;
        delta[i] = change.delta;
    }
    Rcpp::DataFrame output = Rcpp::DataFrame::create(
        Rcpp::_["node"] = node,
        Rcpp::_["baseline"] = baseline,
        Rcpp::_["perturbed"] = perturbed,
        Rcpp::_["delta"] = delta,
        Rcpp::_["affected"] = affected,
        Rcpp::_["causal_distance"] = distance);
    output.attr("execution_time_us") = result.execution_time_us;
    output.attr("visited_unique_scc") = static_cast<double>(result.stats.visited_unique_scc);
    output.attr("scc_evaluations") = static_cast<double>(result.stats.scc_evaluations);
    output.attr("reaction_evaluations") = static_cast<double>(result.stats.reaction_evaluations);
    output.attr("integration_steps") = static_cast<double>(result.stats.integration_steps);
    return output;
}

cellnet::ExecutionStrategy strategy_from_string(const std::string& strategy) {
    if (strategy == "frontier") return cellnet::ExecutionStrategy::Frontier;
    if (strategy == "full") return cellnet::ExecutionStrategy::Full;
    Rcpp::stop("strategy must be 'frontier' or 'full'");
}

cellnet::Perturbation make_perturbation(CellNetHandle& handle,
                                         const std::string& target,
                                         double effect) {
    const auto found = handle.ids.find(target);
    if (found == handle.ids.end()) Rcpp::stop("unknown target node: %s", target);
    cellnet::Perturbation p;
    p.node = found->second;
    p.strength = static_cast<float>(std::min(std::abs(effect), 1.0));
    p.type = effect < 0.0 ? cellnet::PerturbationType::ActivityInhibition
                          : cellnet::PerturbationType::ActivityActivation;
    return p;
}

}  // namespace

// [[Rcpp::export]]
SEXP cn_compile(Rcpp::DataFrame edges, Rcpp::CharacterVector node_names,
                Rcpp::NumericVector initial_values, double dt, int horizon,
                double epsilon, int baseline_max_steps) {
    if (!R_finite(dt) || dt <= 0.0 || horizon <= 0)
        Rcpp::stop("dt and horizon must be positive");
    if (!R_finite(epsilon) || epsilon <= 0.0)
        Rcpp::stop("epsilon must be finite and positive");
    if (baseline_max_steps <= 0)
        Rcpp::stop("baseline_max_steps must be positive");
    const R_xlen_t n_nodes = node_names.size();
    if (initial_values.size() != n_nodes) Rcpp::stop("initial values must match node_names");
    const Rcpp::CharacterVector sources = edges["source"];
    const Rcpp::CharacterVector targets = edges["target"];
    const Rcpp::NumericVector signs = edges["sign"];
    const Rcpp::NumericVector weights = edges["weight"];
    const Rcpp::NumericVector confidences = edges["confidence"];
    const R_xlen_t n_edges = edges.nrows();
    if (targets.size() != n_edges || signs.size() != n_edges || weights.size() != n_edges ||
        confidences.size() != n_edges) Rcpp::stop("edge columns have inconsistent lengths");

    cellnet::NetworkSpec spec;
    spec.nodes.reserve(static_cast<std::size_t>(n_nodes));
    std::unordered_map<std::string, std::uint64_t> ids;
    for (R_xlen_t i = 0; i < n_nodes; ++i) {
        const std::string name = Rcpp::as<std::string>(node_names[i]);
        ids.emplace(name, static_cast<std::uint64_t>(i));
        cellnet::NodeSpec node;
        node.external_id = static_cast<std::uint64_t>(i);
        node.name = name;
        node.initial_value = static_cast<float>(initial_values[i]);
        spec.nodes.push_back(std::move(node));
    }
    spec.reactions.reserve(static_cast<std::size_t>(n_edges));
    for (R_xlen_t i = 0; i < n_edges; ++i) {
        const std::string source = Rcpp::as<std::string>(sources[i]);
        const std::string target = Rcpp::as<std::string>(targets[i]);
        if (!ids.contains(source) || !ids.contains(target)) Rcpp::stop("edge contains unknown node");
        if (!R_finite(signs[i]) || signs[i] == 0.0) Rcpp::stop("edge sign must be nonzero and finite");
        if (!R_finite(weights[i]) || weights[i] < 0.0) Rcpp::stop("edge weight must be finite and nonnegative");
        if (!R_finite(confidences[i]) || confidences[i] < 0.0 || confidences[i] > 1.0)
            Rcpp::stop("confidence must be between 0 and 1");
        cellnet::ReactionSpec reaction;
        reaction.external_id = static_cast<std::uint64_t>(i);
        reaction.primitive = cellnet::Primitive::Modulate;
        reaction.k = static_cast<float>(weights[i]);
        reaction.confidence = static_cast<float>(confidences[i]);
        const auto source_id = static_cast<std::uint64_t>(ids.at(source));
        if (signs[i] > 0.0) reaction.positive_regulators.push_back({source_id, 1.0F, 0.5F, 1.0F});
        else reaction.negative_regulators.push_back({source_id, 1.0F, 0.5F, 1.0F});
        reaction.outputs.push_back({static_cast<std::uint64_t>(ids.at(target)), 1.0F});
        spec.reactions.push_back(std::move(reaction));
    }
    cellnet::SolverConfig config;
    config.dt = static_cast<float>(dt);
    config.epsilon = static_cast<float>(epsilon);
    config.baseline_max_steps = static_cast<std::uint32_t>(baseline_max_steps);
    config.perturbation_steps = static_cast<std::uint32_t>(horizon);
    auto compiled = cellnet::Compiler{}.compile(
        spec, 1U, nullptr, cellnet::TemporalCompileConfig{config.dt, config.perturbation_steps});
    std::vector<std::string> names;
    names.reserve(static_cast<std::size_t>(n_nodes));
    for (const auto& item : node_names) names.push_back(Rcpp::as<std::string>(item));
    Rcpp::XPtr<CellNetHandle> pointer(
        new CellNetHandle(std::move(compiled), config, std::move(names)), true);
    return pointer;
}

// [[Rcpp::export]]
Rcpp::DataFrame cn_perturb(SEXP pointer, std::string target, double effect,
                           std::string strategy) {
    auto& handle = checked_handle(pointer);
    const auto p = make_perturbation(handle, target, effect);
    return materialize(handle, handle.engine.run(p, strategy_from_string(strategy)));
}

// [[Rcpp::export]]
Rcpp::DataFrame cn_perturb_reaction(SEXP pointer, int reaction_id, double effect,
                                   std::string strategy) {
    auto& handle = checked_handle(pointer);
    if (reaction_id < 1 || static_cast<std::size_t>(reaction_id) > handle.engine.network().reactions.size())
        Rcpp::stop("reaction_id is a 1-based edge row index and is out of range");
    if (!R_finite(effect) || effect == 0.0 || std::abs(effect) > 1.0)
        Rcpp::stop("effect must be nonzero and within [-1,1]");
    cellnet::Perturbation p;
    p.reaction = static_cast<cellnet::ReactionId>(reaction_id - 1);
    p.strength = static_cast<float>(std::abs(effect));
    p.type = effect < 0.0 ? cellnet::PerturbationType::ReactionInhibition
                          : cellnet::PerturbationType::ReactionActivation;
    return materialize(handle, handle.engine.run(p, strategy_from_string(strategy)));
}

// [[Rcpp::export]]
Rcpp::List cn_perturb_batch(SEXP pointer, Rcpp::CharacterVector targets,
                            double effect, std::string strategy) {
    auto& handle = checked_handle(pointer);
    std::vector<cellnet::Perturbation> perturbations;
    perturbations.reserve(static_cast<std::size_t>(targets.size()));
    for (const auto& target : targets)
        perturbations.push_back(make_perturbation(handle, Rcpp::as<std::string>(target), effect));
    const auto results = handle.engine.run_batch(perturbations, strategy_from_string(strategy));
    Rcpp::List output(results.size());
    for (std::size_t i = 0; i < results.size(); ++i) output[i] = materialize(handle, results[i]);
    output.attr("names") = targets;
    return output;
}

// [[Rcpp::export]]
void cn_set_context(SEXP pointer, Rcpp::CharacterVector nodes, Rcpp::NumericVector multipliers) {
    if (nodes.size() != multipliers.size()) Rcpp::stop("nodes and multipliers must have same length");
    auto& handle = checked_handle(pointer);
    cellnet::Context context;
    static std::atomic<std::uint64_t> next_context_id{1U};
    context.context_id = next_context_id.fetch_add(1U);
    for (R_xlen_t i = 0; i < nodes.size(); ++i) {
        const auto found = handle.ids.find(Rcpp::as<std::string>(nodes[i]));
        if (found == handle.ids.end()) Rcpp::stop("unknown context node");
        if (!R_finite(multipliers[i]) || multipliers[i] < 0.0) Rcpp::stop("multipliers must be finite and nonnegative");
        context.node_multipliers.push_back({found->second, static_cast<float>(multipliers[i])});
    }
    handle.engine.set_context(context);
    // A context changes the baseline dynamics. Rebuild the baseline before
    // exposing the compiled engine to perturbation/protocol calls.
    handle.engine.initialize_baseline();
}

// [[Rcpp::export]]
Rcpp::List cn_run_protocol(SEXP pointer, Rcpp::DataFrame events, double end_time,
                           Rcpp::NumericVector sample_times, std::string strategy) {
    auto& handle = checked_handle(pointer);
    Rcpp::NumericVector times = events["time"];
    Rcpp::CharacterVector targets = events["target"];
    Rcpp::NumericVector effects = events["effect"];
    if (times.size() != targets.size() || times.size() != effects.size())
        Rcpp::stop("protocol columns have inconsistent lengths");
    std::vector<cellnet::TimedIntervention> schedule;
    schedule.reserve(static_cast<std::size_t>(times.size()));
    for (R_xlen_t i = 0; i < times.size(); ++i) {
        const std::string target = Rcpp::as<std::string>(targets[i]);
        const auto found = handle.ids.find(target);
        if (found == handle.ids.end()) Rcpp::stop("unknown protocol node: %s", target);
        cellnet::TimedIntervention item;
        item.time = times[i];
        item.action.node = found->second;
        item.action.strength = static_cast<float>(std::min(std::abs(effects[i]), 1.0));
        item.action.type = effects[i] < 0.0 ? cellnet::PerturbationType::ActivityInhibition
                                             : cellnet::PerturbationType::ActivityActivation;
        schedule.push_back(item);
    }
    std::vector<double> sample_vector(sample_times.begin(), sample_times.end());
    cellnet::ProtocolResult result;
    handle.engine.reserve_protocol_result(result, sample_vector.size(), schedule.size());
    handle.engine.run_protocol_into(schedule, end_time, sample_vector, result,
                                    strategy_from_string(strategy));
    const R_xlen_t n_time = static_cast<R_xlen_t>(result.samples.size());
    const R_xlen_t n_node = static_cast<R_xlen_t>(handle.names.size());
    Rcpp::NumericMatrix states(n_time, n_node);
    Rcpp::NumericVector output_times(n_time);
    Rcpp::CharacterVector node_names(n_node);
    for (R_xlen_t j = 0; j < n_node; ++j) node_names[j] = handle.names[static_cast<std::size_t>(j)];
    for (R_xlen_t i = 0; i < n_time; ++i) {
        const auto& sample = result.samples[static_cast<std::size_t>(i)];
        output_times[i] = sample.time;
        for (R_xlen_t j = 0; j < n_node; ++j)
            states(i, j) = sample.state[static_cast<std::size_t>(j)];
    }
    states.attr("dimnames") = Rcpp::List::create(R_NilValue, node_names);
    return Rcpp::List::create(Rcpp::_["times"] = output_times,
                              Rcpp::_["states"] = states,
                              Rcpp::_["numerical_error"] = result.numerical_error);
}
