
@safetestset track_Petri_smoke = "Petri smoke" begin
    using Distributions
    using CompetingClocks
    using Random
    using Base
    rng = Xoshiro(3242234)
    tw = Petri{Int,Float64}()
    # Show that distributions with very different rates
    # are all sampled equally by Petri.
    enable!(tw, 3, Exponential(100.0), 0.0, 0.0, rng)
    @test length(tw) == 1 && 3 ∈ keys(tw)
    enable!(tw, 4, Exponential(0.001), 0.0, 3.0, rng)
    @test length(tw) == 2 && 4 ∈ keys(tw)
    enable!(tw, 7, Exponential(0.001), 5.0, 5.0, rng)
    @test length(tw) == 3 && 7 ∈ keys(tw)

    counts = Dict{Int,Int}(3 => 0, 4 => 0, 7 => 0)
    for i in 1:1000
        when = 100 * rand(rng)
        time_out, which = next(tw, when, rng)
        counts[which] += 1
        @assert abs(time_out - when - 1) < 1e-9
    end
    hi, lo = (maximum(values(counts)), minimum(values(counts)))
    @test (hi - lo) / lo < 0.3

    disable!(tw, 4, 9.0)
    @test length(tw) == 2 && 4 ∉ keys(tw)

    dst = Petri{Int,Float64}()
    enable!(dst, 11, Exponential(), 5.0, 5.0, rng)
    copy!(dst, tw)
    @test length(tw) == 2 && 11 ∉ keys(tw)
end


@safetestset track_Petri_more = "Petri more functions" begin
    using Distributions
    using CompetingClocks
    using Random
    using Base
    rng = Xoshiro(3242234)
    tw = Petri{Int,Float64}()
    # Show that distributions with very different rates
    # are all sampled equally by Petri.
    enable!(tw, 3, Exponential(100.0), 0.0, 0.0, rng)
    @test tw[3].clock == 3
    @test tw[3].distribution == Exponential(100.0)
    @test haskey(tw, 3)
    @test !haskey(tw, 4)
end
