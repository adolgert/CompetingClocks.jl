# config.jl - Phase 2 of the gauntlet: the configuration system.
#
# Defines `TestConfiguration` (one CONDITION cell of the matrix), the
# combinatorial design generator `generate_full_design`, and the sampler
# capability filters `filter_for_sampler` / `exponential_only_configs`.
#
# A CONDITION is the cross product of three axes:
#   - distribution family: :exponential | :weibull | :shifted
#   - memory:              :forget | :remember   (TravelMemory)
#   - graph:               :cycle  | :complete   (TravelGraph)
#
# Each `TestConfiguration` maps to a concrete `TravelConfig` plus a model builder
# via `travel_config_and_builder`. The three families cover the plan's
# "distribution family" discussion:
#   - :exponential -> exponential per-destination rates, no enabling delay.
#   - :weibull     -> Weibull per-destination rates (weibull_model_builder), the
#                     non-exponential family where the two-sample test has power.
#   - :shifted     -> exponential rates with a RIGHT enabling-time shift
#                     (TravelRateDelay.right): a clock cannot fire until a later
#                     time. This is the plan's "shift in enabling times" and is a
#                     NON-exponential condition (excluded for exponential-only
#                     samplers, per gauntlet.md).
#
# This file assumes travel.jl (TravelModel) and runner.jl (default_model_builder,
# weibull_model_builder, and the MultiSampler specs) are already in scope.

using CompetingClocks
using CompetingClocks: NextReactionMethod, FirstReactionMethod, FirstToFireMethod,
    DirectMethod, RejectionMethod, PartialPropensityMethod


"""
    TestConfiguration(family, memory, graph; state_cnt, history_steps, doob_steps)

One condition cell of the gauntlet matrix. `family ∈ (:exponential,:weibull,
:shifted)`, `memory ∈ (:forget,:remember)`, `graph ∈ (:cycle,:complete)`.
`state_cnt` is the number of sites; `history_steps` the length of the fixed
history for the two-sample test; `doob_steps` the Doob-Meyer trajectory length.
"""
struct TestConfiguration
    family::Symbol
    memory::Symbol
    graph::Symbol
    state_cnt::Int
    history_steps::Int
    doob_steps::Int
end

function TestConfiguration(family::Symbol, memory::Symbol, graph::Symbol;
        state_cnt::Int=5, history_steps::Int=5, doob_steps::Int=1000)
    return TestConfiguration(family, memory, graph, state_cnt, history_steps, doob_steps)
end


"A short, stable key for a condition, e.g. \"weibull/forget/cycle\"."
config_key(tc::TestConfiguration) = string(tc.family, "/", tc.memory, "/", tc.graph)

"A human-readable label for the verdict table / report."
config_label(tc::TestConfiguration) =
    string(tc.family, ", memory=", tc.memory, ", graph=", tc.graph)


"""
    parse_config_key(key; state_cnt, history_steps, doob_steps) -> TestConfiguration

Inverse of `config_key`. Parses "family/memory/graph" back into a
`TestConfiguration` so a single flagged cell can be reproduced from a command.
"""
function parse_config_key(key::AbstractString;
        state_cnt::Int=5, history_steps::Int=5, doob_steps::Int=1000)
    parts = split(key, "/")
    length(parts) == 3 || error("bad config key: $key (want family/memory/graph)")
    return TestConfiguration(Symbol(parts[1]), Symbol(parts[2]), Symbol(parts[3]),
        state_cnt, history_steps, doob_steps)
end


"""
    travel_config_and_builder(tc) -> (TravelConfig, model_builder)

Map a `TestConfiguration` to the concrete `TravelConfig` and the model-builder
function the runner uses to construct the `Travel` model.
"""
function travel_config_and_builder(tc::TestConfiguration)
    memory = tc.memory == :forget ? TravelMemory.forget : TravelMemory.remember
    graph = tc.graph == :cycle ? TravelGraph.cycle : TravelGraph.complete
    if tc.family == :exponential
        cfg = TravelConfig(memory, graph, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none)
        return cfg, default_model_builder
    elseif tc.family == :shifted
        cfg = TravelConfig(memory, graph, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.right)
        return cfg, default_model_builder
    elseif tc.family == :weibull
        # weibull_model_builder ignores config.dist/count/delay (it builds pure
        # Weibull per destination with no delay) but honors graph and memory.
        cfg = TravelConfig(memory, graph, TravelRateDist.general,
            TravelRateCount.destination, TravelRateDelay.none)
        return cfg, weibull_model_builder
    else
        error("unknown distribution family: $(tc.family)")
    end
end


# --- Combinatorial design ---------------------------------------------------

const DEFAULT_FAMILIES = (:exponential, :weibull, :shifted)
const DEFAULT_MEMORIES = (:forget, :remember)
const DEFAULT_GRAPHS = (:cycle, :complete)


"""
    generate_full_design(; families, memories, graphs, state_cnt, history_steps, doob_steps)

Full combinatorial design over the three condition axes. Returns a
`Vector{TestConfiguration}` in a fixed, deterministic order (families outermost,
then memories, then graphs) so cell seeds are reproducible.
"""
function generate_full_design(;
        families=DEFAULT_FAMILIES, memories=DEFAULT_MEMORIES, graphs=DEFAULT_GRAPHS,
        state_cnt::Int=5, history_steps::Int=5, doob_steps::Int=1000)
    configs = TestConfiguration[]
    for family in families, memory in memories, graph in graphs
        push!(configs, TestConfiguration(family, memory, graph,
            state_cnt, history_steps, doob_steps))
    end
    return configs
end


"""
    exponential_only_configs(configs) -> Vector{TestConfiguration}

Keep only the conditions an exponential-only sampler (Direct*, RSSA, PSSACR) can
run: the `:exponential` distribution family. The `:weibull` and `:shifted`
families are non-exponential and are excluded.
"""
exponential_only_configs(configs) =
    filter(c -> c.family == :exponential, configs)


"""
    sampler_capability(sampler_key) -> Symbol

`:exponential` for samplers restricted to exponential distributions,
`:general` otherwise. Driven by the sampler's registry entry (see
`matrix_samplers`).
"""
function sampler_capability(sampler_key::AbstractString)
    for (key, _spec, cap) in matrix_samplers()
        key == sampler_key && return cap
    end
    error("unknown sampler key: $sampler_key")
end


"""
    filter_for_sampler(configs, capability) -> Vector{TestConfiguration}

Restrict a design to the conditions a sampler of the given `capability`
(`:exponential` or `:general`) can run. Exponential-only samplers get only the
exponential family; general samplers get everything.
"""
function filter_for_sampler(configs, capability::Symbol)
    capability == :exponential && return exponential_only_configs(configs)
    return collect(configs)
end


# --- Sampler registry -------------------------------------------------------

"""
    matrix_samplers() -> Vector{Tuple{String,spec,Symbol}}

The samplers under test, as `(display_key, sampler_spec, capability)` triples.
`capability` is `:exponential` (exponential-only) or `:general`.
FirstReactionMethod is both the trusted reference AND its own SUT row: as an SUT
its two-sample comparison is FirstReaction-vs-FirstReaction, a built-in null
control.
"""
function matrix_samplers()
    return Tuple{String,Any,Symbol}[
        ("FirstReactionMethod", FirstReactionMethod(), :general),   # null-control SUT
        ("NextReactionMethod", NextReactionMethod(), :general),
        ("FirstToFireMethod", FirstToFireMethod(), :general),
        ("DirectMethod(:keep,:array)", DirectMethod(:keep, :array), :exponential),
        ("DirectMethod(:keep,:tree)", DirectMethod(:keep, :tree), :exponential),
        ("DirectMethod(:remove,:array)", DirectMethod(:remove, :array), :exponential),
        ("DirectMethod(:remove,:tree)", DirectMethod(:remove, :tree), :exponential),
        ("RejectionMethod(RSSA)", RejectionMethod(), :exponential),
        ("PartialPropensityMethod(PSSACR)", PartialPropensityMethod(), :exponential),
        ("MultiSampler(all=>CNR)", multisampler_control_spec(), :general),
        ("MultiSampler(even=>CNR,odd=>CNR)", multisampler_split_spec(), :general),
        ("MultiSampler(even=>CNR,odd=>FirstToFire)", multisampler_mixed_spec(), :general),
    ]
end


"Look up a sampler spec by its display key (for single-cell reproduction)."
function sampler_spec_for(sampler_key::AbstractString)
    for (key, spec, _cap) in matrix_samplers()
        key == sampler_key && return spec
    end
    error("unknown sampler key: $sampler_key")
end
