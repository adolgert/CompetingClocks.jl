using SafeTestsets


@safetestset OrderedSample_smoke = "nr transition can be created and compared" begin
    using CompetingClocks: OrderedSample

    a = OrderedSample{Int,Float64}(3, 2.2)
    b = OrderedSample{Int,Float64}(1, 2.5)
    c = OrderedSample{Int,Float64}(2, 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    a = OrderedSample{String,Float64}("S", 2.2)
    b = OrderedSample{String,Float64}("I", 2.5)
    c = OrderedSample{String,Float64}("R", 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c

    # You can index it with a tuple if you want.
    a = OrderedSample{Tuple{String,Int},Float64}(("S", 3), 2.2)
    b = OrderedSample{Tuple{String,Int},Float64}(("I", 7), 2.5)
    c = OrderedSample{Tuple{String,Int},Float64}(("R", 4), 2.5)
    @test a < b
    @test isless(a, b)
    @test b > a
    @test b == c
end


@safetestset OrderedSample_heapable = "nr transition can be used in a heap" begin
    using DataStructures
    using CompetingClocks: OrderedSample

    # This shows that the custom isless() operator is what we need in order to
    # use this data structure for sampling.
    heap = MutableBinaryMinHeap{OrderedSample{Int,Float64}}()
    push!(heap, OrderedSample{Int,Float64}(3, 2.2))
    push!(heap, OrderedSample{Int,Float64}(1, 3.9))
    push!(heap, OrderedSample{Int,Float64}(7, 0.05))
    push!(heap, OrderedSample{Int,Float64}(5, 1.7))
    @test pop!(heap).key == 7
    @test pop!(heap).key == 5
    @test pop!(heap).key == 3
    @test pop!(heap).key == 1
end
