using Fleck
using SafeTestsets
using Test

@testset "Fleck.jl" begin
    include("test_prefixsearch.jl")
    include("test_direct.jl")
    include("test_vas.jl")
    #include("test_vas_integrate.jl")
end
