using Fleck
using SafeTestsets
using Test

# Not tests. These are helper functions for tests.
include("test_utility.jl")

@testset "test_prefixsearch.jl" begin
    include("test_prefixsearch.jl")
end

@testset "test_direct.jl" begin
    include("test_direct.jl")
end

@testset "test_vas.jl" begin
    include("test_vas.jl")
end

@testset "test_vas_integrate.jl" begin
    include("test_vas_integrate.jl")
end
