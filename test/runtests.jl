using Fleck
using SafeTestsets
using Test

@testset "test_prefixsearch.jl" begin
    include("test_prefixsearch.jl")
end

@testset "test_direct.jl" begin
    include("test_direct.jl")
end

@testset "test_vas.jl" begin
    include("test_vas.jl")
end

#include("test_vas_integrate.jl")
