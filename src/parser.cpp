#include "cellnet/parser.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cellnet {
namespace {

std::string trim(std::string value) {
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char character) {
        return std::isspace(character) != 0;
    });
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char character) {
        return std::isspace(character) != 0;
    }).base();
    if (first >= last) {
        return {};
    }
    return std::string(first, last);
}

std::string normalized_token(const std::string& text) {
    std::string result;
    result.reserve(text.size());
    for (const char raw_character : text) {
        const auto character = static_cast<unsigned char>(raw_character);
        if (std::isalnum(character) != 0) {
            result.push_back(static_cast<char>(std::toupper(character)));
        }
    }
    return result;
}

std::vector<std::string> parse_csv_line(const std::string& line, const std::filesystem::path& path,
                                        std::size_t line_number) {
    std::vector<std::string> fields;
    std::string field;
    bool quoted = false;
    for (std::size_t index = 0; index < line.size(); ++index) {
        const char character = line[index];
        if (quoted) {
            if (character == '"') {
                if (index + 1U < line.size() && line[index + 1U] == '"') {
                    field.push_back('"');
                    ++index;
                } else {
                    quoted = false;
                }
            } else {
                field.push_back(character);
            }
            continue;
        }
        if (character == '"') {
            if (!field.empty()) {
                throw ParseError(path.string() + ":" + std::to_string(line_number) +
                                 ": quote inside unquoted field");
            }
            quoted = true;
        } else if (character == ',') {
            fields.push_back(trim(std::move(field)));
            field.clear();
        } else {
            field.push_back(character);
        }
    }
    if (quoted) {
        throw ParseError(path.string() + ":" + std::to_string(line_number) +
                         ": unterminated quoted field");
    }
    fields.push_back(trim(std::move(field)));
    return fields;
}

using HeaderMap = std::unordered_map<std::string, std::size_t>;
using RowCallback = std::function<void(const HeaderMap&, const std::vector<std::string>&, std::size_t)>;

void read_csv(const std::filesystem::path& path, bool required, const RowCallback& callback) {
    std::ifstream input(path);
    if (!input) {
        if (required) {
            throw ParseError("cannot open required CSV: " + path.string());
        }
        return;
    }
    std::string line;
    std::size_t line_number = 0U;
    HeaderMap headers;
    std::size_t column_count = 0U;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line_number == 1U && line.size() >= 3U &&
            static_cast<unsigned char>(line[0]) == 0xEFU &&
            static_cast<unsigned char>(line[1]) == 0xBBU &&
            static_cast<unsigned char>(line[2]) == 0xBFU) {
            line.erase(0U, 3U);
        }
        if (trim(line).empty()) {
            continue;
        }
        const auto fields = parse_csv_line(line, path, line_number);
        if (headers.empty()) {
            column_count = fields.size();
            for (std::size_t index = 0; index < fields.size(); ++index) {
                const std::string key = normalized_token(fields[index]);
                if (key.empty() || !headers.emplace(key, index).second) {
                    throw ParseError(path.string() + ": invalid or duplicate header");
                }
            }
            continue;
        }
        if (fields.size() != column_count) {
            throw ParseError(path.string() + ":" + std::to_string(line_number) +
                             ": column count does not match header");
        }
        callback(headers, fields, line_number);
    }
    if (headers.empty()) {
        throw ParseError("CSV has no header: " + path.string());
    }
}

const std::string& field(const HeaderMap& headers, const std::vector<std::string>& row,
                         const std::string& name, const std::filesystem::path& path,
                         std::size_t line_number, bool required = true) {
    const auto found = headers.find(normalized_token(name));
    if (found == headers.end()) {
        if (!required) {
            static const std::string empty;
            return empty;
        }
        throw ParseError(path.string() + ": missing required column " + name);
    }
    const std::string& value = row[found->second];
    if (required && value.empty()) {
        throw ParseError(path.string() + ":" + std::to_string(line_number) +
                         ": empty required field " + name);
    }
    return value;
}

std::uint64_t parse_u64(const std::string& text, const std::filesystem::path& path,
                        std::size_t line_number, const char* label) {
    try {
        std::size_t consumed = 0U;
        const unsigned long long value = std::stoull(text, &consumed, 10);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return static_cast<std::uint64_t>(value);
    } catch (const std::exception&) {
        throw ParseError(path.string() + ":" + std::to_string(line_number) + ": invalid " + label +
                         ": " + text);
    }
}

float parse_float(const std::string& text, const std::filesystem::path& path,
                  std::size_t line_number, const char* label, float default_value,
                  bool required = false) {
    if (text.empty() && !required) {
        return default_value;
    }
    try {
        std::size_t consumed = 0U;
        const float value = std::stof(text, &consumed);
        if (consumed != text.size()) {
            throw std::invalid_argument("trailing characters");
        }
        return value;
    } catch (const std::exception&) {
        throw ParseError(path.string() + ":" + std::to_string(line_number) + ": invalid " + label +
                         ": " + text);
    }
}

ReactionSpec& find_reaction(NetworkSpec& spec,
                            const std::unordered_map<std::uint64_t, std::size_t>& reaction_index,
                            std::uint64_t external_id, const std::filesystem::path& path,
                            std::size_t line_number) {
    const auto found = reaction_index.find(external_id);
    if (found == reaction_index.end()) {
        throw ParseError(path.string() + ":" + std::to_string(line_number) +
                         ": unknown reaction_id " + std::to_string(external_id));
    }
    return spec.reactions[found->second];
}

}  // namespace

NodeType parse_node_type(const std::string& text) {
    const std::string token = normalized_token(text);
    if (token == "GENE") return NodeType::Gene;
    if (token == "RNA" || token == "MRNA") return NodeType::RNA;
    if (token == "PROTEIN") return NodeType::Protein;
    if (token == "PROTEINSTATE") return NodeType::ProteinState;
    if (token == "COMPLEX") return NodeType::Complex;
    if (token == "METABOLITE") return NodeType::Metabolite;
    if (token == "PHENOTYPE") return NodeType::Phenotype;
    if (token == "EXTERNALSIGNAL") return NodeType::ExternalSignal;
    throw ParseError("unknown node type: " + text);
}

Compartment parse_compartment(const std::string& text) {
    const std::string token = normalized_token(text);
    if (token.empty() || token == "UNKNOWN") return Compartment::Unknown;
    if (token == "EXTRACELLULAR") return Compartment::Extracellular;
    if (token == "MEMBRANE") return Compartment::Membrane;
    if (token == "CYTOSOL" || token == "CYTOPLASM") return Compartment::Cytosol;
    if (token == "NUCLEUS" || token == "NUCLEAR") return Compartment::Nucleus;
    if (token == "MITOCHONDRIA" || token == "MITOCHONDRION") return Compartment::Mitochondria;
    if (token == "ER" || token == "ENDOPLASMICRETICULUM") return Compartment::ER;
    if (token == "GOLGI") return Compartment::Golgi;
    throw ParseError("unknown compartment: " + text);
}

Primitive parse_primitive(const std::string& text) {
    const std::string token = normalized_token(text);
    if (token == "TRANSFER" || token == "PHOSPHORYLATION" ||
        token == "DEPHOSPHORYLATION" || token == "STATECHANGE") return Primitive::Transfer;
    if (token == "PRODUCE" || token == "SYNTHESIS" || token == "TRANSCRIPTION" ||
        token == "TRANSLATION") return Primitive::Produce;
    if (token == "DESTROY" || token == "DEGRADATION") return Primitive::Destroy;
    if (token == "BIND" || token == "BINDING") return Primitive::Bind;
    if (token == "UNBIND" || token == "DISSOCIATION") return Primitive::Unbind;
    if (token == "MODULATE" || token == "REGULATION") return Primitive::Modulate;
    if (token == "MOVE" || token == "TRANSLOCATION" || token == "TRANSPORT") return Primitive::Move;
    if (token == "DELAY") return Primitive::Delay;
    throw ParseError("unknown primitive: " + text);
}

RegulationOp parse_regulation_op(const std::string& text) {
    const std::string token = normalized_token(text);
    if (token.empty() || token == "NONE") return RegulationOp::None;
    if (token == "SUM" || token == "WEIGHTEDSUM") return RegulationOp::Sum;
    if (token == "PRODUCT") return RegulationOp::Product;
    if (token == "MIN") return RegulationOp::Min;
    if (token == "MAX") return RegulationOp::Max;
    if (token == "HILLPOS" || token == "HILLPOSITIVE") return RegulationOp::HillPositive;
    if (token == "HILLNEG" || token == "HILLNEGATIVE") return RegulationOp::HillNegative;
    if (token == "AND" || token == "CONTINUOUSAND") return RegulationOp::ContinuousAnd;
    if (token == "OR" || token == "CONTINUOUSOR") return RegulationOp::ContinuousOr;
    if (token == "THRESHOLD") return RegulationOp::Threshold;
    if (token == "XOR") return RegulationOp::Xor;
    throw ParseError("unknown regulation operator: " + text);
}

NetworkSpec parse_csv_directory(const std::filesystem::path& directory) {
    NetworkSpec spec;
    std::unordered_map<std::uint64_t, std::size_t> reaction_index;

    const auto nodes_path = directory / "nodes.csv";
    read_csv(nodes_path, true, [&](const HeaderMap& headers, const std::vector<std::string>& row,
                                   std::size_t line_number) {
        NodeSpec node;
        node.external_id = parse_u64(field(headers, row, "id", nodes_path, line_number), nodes_path,
                                     line_number, "node id");
        node.name = field(headers, row, "name", nodes_path, line_number);
        node.type = parse_node_type(field(headers, row, "type", nodes_path, line_number));
        node.entity = field(headers, row, "entity", nodes_path, line_number, false);
        node.state = field(headers, row, "state", nodes_path, line_number, false);
        node.compartment = parse_compartment(
            field(headers, row, "compartment", nodes_path, line_number, false));
        node.initial_value = parse_float(
            field(headers, row, "initial", nodes_path, line_number), nodes_path, line_number,
            "initial value", 0.0F, true);
        spec.nodes.push_back(std::move(node));
    });

    const auto reactions_path = directory / "reactions.csv";
    read_csv(reactions_path, true,
             [&](const HeaderMap& headers, const std::vector<std::string>& row,
                 std::size_t line_number) {
        ReactionSpec reaction;
        reaction.external_id = parse_u64(
            field(headers, row, "id", reactions_path, line_number), reactions_path, line_number,
            "reaction id");
        if (reaction_index.contains(reaction.external_id)) {
            throw ParseError(reactions_path.string() + ":" + std::to_string(line_number) +
                             ": duplicate reaction id");
        }
        reaction.primitive =
            parse_primitive(field(headers, row, "primitive", reactions_path, line_number));
        reaction.k = parse_float(field(headers, row, "k", reactions_path, line_number, false),
                                 reactions_path, line_number, "k", 1.0F);
        reaction.Km = parse_float(field(headers, row, "Km", reactions_path, line_number, false),
                                  reactions_path, line_number, "Km", 0.5F);
        reaction.hill_n = parse_float(
            field(headers, row, "hill_n", reactions_path, line_number, false), reactions_path,
            line_number, "hill_n", 1.0F);
        reaction.confidence = parse_float(
            field(headers, row, "confidence", reactions_path, line_number, false), reactions_path,
            line_number, "confidence", 1.0F);
        reaction.context_weight = parse_float(
            field(headers, row, "context_weight", reactions_path, line_number, false),
            reactions_path, line_number, "context_weight", 1.0F);
        const std::string& input_op =
            field(headers, row, "input_op", reactions_path, line_number, false);
        if (!input_op.empty()) {
            reaction.input_op = parse_regulation_op(input_op);
        }
        reaction_index.emplace(reaction.external_id, spec.reactions.size());
        spec.reactions.push_back(std::move(reaction));
    });

    const auto inputs_path = directory / "reaction_inputs.csv";
    read_csv(inputs_path, false, [&](const HeaderMap& headers, const std::vector<std::string>& row,
                                    std::size_t line_number) {
        const std::uint64_t reaction_id = parse_u64(
            field(headers, row, "reaction_id", inputs_path, line_number), inputs_path, line_number,
            "reaction_id");
        const std::uint64_t node_id = parse_u64(
            field(headers, row, "node_id", inputs_path, line_number), inputs_path, line_number,
            "node_id");
        const float stoich = parse_float(
            field(headers, row, "stoich", inputs_path, line_number, false), inputs_path, line_number,
            "stoich", 1.0F);
        find_reaction(spec, reaction_index, reaction_id, inputs_path, line_number)
            .inputs.push_back(WeightedNodeSpec{node_id, stoich});
    });

    const auto outputs_path = directory / "reaction_outputs.csv";
    read_csv(outputs_path, false, [&](const HeaderMap& headers, const std::vector<std::string>& row,
                                     std::size_t line_number) {
        const std::uint64_t reaction_id = parse_u64(
            field(headers, row, "reaction_id", outputs_path, line_number), outputs_path, line_number,
            "reaction_id");
        const std::uint64_t node_id = parse_u64(
            field(headers, row, "node_id", outputs_path, line_number), outputs_path, line_number,
            "node_id");
        const float stoich = parse_float(
            field(headers, row, "stoich", outputs_path, line_number, false), outputs_path,
            line_number, "stoich", 1.0F);
        find_reaction(spec, reaction_index, reaction_id, outputs_path, line_number)
            .outputs.push_back(WeightedNodeSpec{node_id, stoich});
    });

    const auto catalysts_path = directory / "catalysts.csv";
    read_csv(catalysts_path, false,
             [&](const HeaderMap& headers, const std::vector<std::string>& row,
                 std::size_t line_number) {
        const std::uint64_t reaction_id = parse_u64(
            field(headers, row, "reaction_id", catalysts_path, line_number), catalysts_path,
            line_number, "reaction_id");
        const std::uint64_t node_id = parse_u64(
            field(headers, row, "node_id", catalysts_path, line_number), catalysts_path, line_number,
            "node_id");
        find_reaction(spec, reaction_index, reaction_id, catalysts_path, line_number)
            .catalysts.push_back(node_id);
    });

    const auto regulation_path = directory / "regulation.csv";
    read_csv(regulation_path, false,
             [&](const HeaderMap& headers, const std::vector<std::string>& row,
                 std::size_t line_number) {
        const std::uint64_t reaction_id = parse_u64(
            field(headers, row, "reaction_id", regulation_path, line_number), regulation_path,
            line_number, "reaction_id");
        const std::uint64_t node_id = parse_u64(
            field(headers, row, "node_id", regulation_path, line_number), regulation_path,
            line_number, "node_id");
        const std::string type = normalized_token(
            field(headers, row, "type", regulation_path, line_number));
        const RegulationOp op = parse_regulation_op(
            field(headers, row, "operator", regulation_path, line_number));
        RegulatorSpec regulator;
        regulator.node_external_id = node_id;
        regulator.K = parse_float(field(headers, row, "K", regulation_path, line_number, false),
                                  regulation_path, line_number, "K", 0.5F);
        regulator.n = parse_float(field(headers, row, "n", regulation_path, line_number, false),
                                  regulation_path, line_number, "n", 1.0F);
        regulator.weight = parse_float(
            field(headers, row, "weight", regulation_path, line_number, false), regulation_path,
            line_number, "weight", 1.0F);
        ReactionSpec& reaction =
            find_reaction(spec, reaction_index, reaction_id, regulation_path, line_number);
        if (type == "POSITIVE" || type == "POS") {
            if (!reaction.positive_regulators.empty() && reaction.pos_reg_op != op) {
                throw ParseError(regulation_path.string() + ":" + std::to_string(line_number) +
                                 ": mixed positive operators in one reaction");
            }
            reaction.pos_reg_op = op;
            reaction.positive_regulators.push_back(regulator);
        } else if (type == "NEGATIVE" || type == "NEG") {
            if (!reaction.negative_regulators.empty() && reaction.neg_reg_op != op) {
                throw ParseError(regulation_path.string() + ":" + std::to_string(line_number) +
                                 ": mixed negative operators in one reaction");
            }
            reaction.neg_reg_op = op;
            reaction.negative_regulators.push_back(regulator);
        } else {
            throw ParseError(regulation_path.string() + ":" + std::to_string(line_number) +
                             ": regulation type must be positive or negative");
        }
    });

    return spec;
}

}  // namespace cellnet
