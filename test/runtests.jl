using CompetingClocks
using SafeTestsets
using Test

# Not tests. These are helper functions for tests.
include("vas.jl")
include("test_utility.jl")
include("erlang_loss.jl")

# All available test files
all_tests = [
    "test_package.jl",
    "test_combinednr.jl",
    "test_context.jl",
    "test_direct.jl",
    "test_firstreaction.jl",
    "test_firsttofire.jl",
    "test_interface.jl",
    "test_keyedprefixsearch.jl",
    "test_likelihood.jl",
    "test_multiple_direct.jl",
    "test_nrtransition.jl",
    "test_petri.jl",
    "test_prefixsearch.jl",
    "test_pssa_cr.jl",
    "test_rssa.jl",
    "test_sampler.jl",
    "test_samplerspec.jl",
    "test_sampler_builder.jl",
    "test_setofsets.jl",
    "test_track.jl",
    "test_vas_integrate.jl",
    "test_vas.jl",
    "test_with_common_random.jl",
    "gauntlet/test_travel.jl",
]

# Filter tests based on command-line arguments (ARGS)
# If no arguments provided, run all tests
# Otherwise, run only tests whose names contain any of the arguments
selected_tests = if isempty(ARGS)
    all_tests
else
    filtered = filter(test_file -> any(arg -> occursin(arg, test_file), ARGS), all_tests)
    if isempty(filtered)
        @warn "No tests matched arguments: $(ARGS). Running all tests."
        all_tests
    else
        filtered
    end
end

@info "Running $(length(selected_tests)) of $(length(all_tests)) test file(s)"
if length(selected_tests) < length(all_tests)
    @info "Selected tests: $(selected_tests)"
end

# Run selected tests
for test_file in selected_tests
    @testset "$test_file" begin
        include(test_file)
    end
end
