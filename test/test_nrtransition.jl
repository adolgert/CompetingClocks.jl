using SafeTestsets


@safetestset nrtransition_smoke = "nr transition can be created and compared" begin
    using Fleck: NRTransition

    a = NRTransition{Int}(3, 2.2)
    b = NRTransition{Int}(1, 2.5)
    c = NRTransition{Int}(2, 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    a = NRTransition{String}("S", 2.2)
    b = NRTransition{String}("I", 2.5)
    c = NRTransition{String}("R", 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    # You can index it with a tuple if you want.
    a = NRTransition{Tuple{String,Int}}(("S", 3), 2.2)
    b = NRTransition{Tuple{String,Int}}(("I", 7), 2.5)
    c = NRTransition{Tuple{String,Int}}(("R", 4), 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c
end


@safetestset nrtransition_heapable = "nr transition can be used in a heap" begin
    using DataStructures
    using Fleck: NRTransition

    # This shows that the custom isless() operator is what we need in order to
    # use this data structure for sampling.
    heap = MutableBinaryMinHeap{NRTransition{Int}}()
    push!(heap, NRTransition{Int}(3, 2.2))
    push!(heap, NRTransition{Int}(1, 3.9))
    push!(heap, NRTransition{Int}(7, 0.05))
    push!(heap, NRTransition{Int}(5, 1.7))
    @test pop!(heap).key == 7
    @test pop!(heap).key == 5
    @test pop!(heap).key == 3
    @test pop!(heap).key == 1
end
