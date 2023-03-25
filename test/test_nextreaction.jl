using SafeTestsets


@safetestset NextReaction_handle_always_positive = "heap handle is always positive" begin
    using DataStructures
    using Fleck: OrderedSample

    # The NextReaction algorithm relies on the heap handle always being positive
    # so this test checks that is the case.
    heap = MutableBinaryMinHeap{OrderedSample{Int}}()
    enabled = Set{Int}()
    for i in 1:10000
        if rand() < 0.2 && length(heap) > 0
            v, handle = top_with_handle(heap)
            @test(handle > 0)
            delete!(enabled, handle)
        elseif rand() < 0.4 && length(heap) > 0 && length(enabled) > 0
            modify = rand(enabled)
            update!(heap, modify, OrderedSample{Int}(i, rand()))
        else
            handle = push!(heap, OrderedSample{Int}(i, rand()))
            push!(enabled, handle)
            @test(handle > 0)
        end
    end
end


@safetestset NextReactionSmoke = "next reaction does basic things" begin
    using Distributions
    using Random
    using Fleck: NextReaction, next, enable!, disable!

    rng = MersenneTwister(349827)
    for i in 1:100
        sampler = NextReaction{String}()
        @test next(sampler, 3.0, rng)[2] === nothing
        enable!(sampler, "walk home", Exponential(1.5), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] == "walk home"
        enable!(sampler, "run", Gamma(1, 3), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        enable!(sampler, "walk to sandwich shop", Weibull(2, 1), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run", "walk to sandwich shop"]
        disable!(sampler, "walk to sandwich shop", 1.7)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
    end
end
