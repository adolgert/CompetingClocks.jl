using Fleck
using SafeTestsets
using Test

# Not tests. These are helper functions for tests.
include("test_utility.jl")

@testset "test_prefixsearch.jl" begin
    include("test_prefixsearch.jl")
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


@testset "test_nextreaction.jl" begin
    include("test_nextreaction.jl")
end


@testset "test_firsttofire.jl" begin
    include("test_firsttofire.jl")
end


@testset "test_vas.jl" begin
    include("test_vas.jl")
end

@testset "test_vas_integrate.jl" begin
    include("test_vas_integrate.jl")
end
