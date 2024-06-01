using CompetingClocks
using SafeTestsets
using Test

# Not tests. These are helper functions for tests.
include("vas.jl")
include("test_utility.jl")

@testset "test_commonrandom.jl" begin
    include("test_commonrandom.jl")
end


@testset "test_prefixsearch.jl" begin
    include("test_prefixsearch.jl")
end


@testset "test_keyedprefixsearch.jl" begin
    include("test_keyedprefixsearch.jl")
end

@testset "test_track.jl" begin
    include("test_track.jl")
end

@testset "test_combinednr.jl" begin
    include("test_combinednr.jl")
end


@testset "test_firstreaction.jl" begin
    include("test_firstreaction.jl")
end


@testset "test_nrtransition.jl" begin
    include("test_nrtransition.jl")
end


@testset "test_direct.jl" begin
    include("test_direct.jl")
end


@testset "test_firsttofire.jl" begin
    include("test_firsttofire.jl")
end


@testset "test_sampler.jl" begin
    include("test_sampler.jl")
end


@testset "test_multiple_direct.jl" begin
    include("test_multiple_direct.jl")
end


@testset "test_vas.jl" begin
    include("test_vas.jl")
end


@testset "test_vas_integrate.jl" begin
    include("test_vas_integrate.jl")
end
