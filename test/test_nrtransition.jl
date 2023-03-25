using SafeTestsets


@safetestset OrderedSample_smoke = "nr transition can be created and compared" begin
    using Fleck: OrderedSample

    a = OrderedSample{Int}(3, 2.2)
    b = OrderedSample{Int}(1, 2.5)
    c = OrderedSample{Int}(2, 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    a = OrderedSample{String}("S", 2.2)
    b = OrderedSample{String}("I", 2.5)
    c = OrderedSample{String}("R", 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    # You can index it with a tuple if you want.
    a = OrderedSample{Tuple{String,Int}}(("S", 3), 2.2)
    b = OrderedSample{Tuple{String,Int}}(("I", 7), 2.5)
    c = OrderedSample{Tuple{String,Int}}(("R", 4), 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c
end


@safetestset OrderedSample_heapable = "nr transition can be used in a heap" begin
    using DataStructures
    using Fleck: OrderedSample

    # This shows that the custom isless() operator is what we need in order to
    # use this data structure for sampling.
    heap = MutableBinaryMinHeap{OrderedSample{Int}}()
    push!(heap, OrderedSample{Int}(3, 2.2))
    push!(heap, OrderedSample{Int}(1, 3.9))
    push!(heap, OrderedSample{Int}(7, 0.05))
    push!(heap, OrderedSample{Int}(5, 1.7))
    @test pop!(heap).key == 7
    @test pop!(heap).key == 5
    @test pop!(heap).key == 3
    @test pop!(heap).key == 1
end
