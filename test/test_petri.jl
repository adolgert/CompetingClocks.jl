
@safetestset track_Petri_smoke = "Petri smoke" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: Petri
    using Random
    using Base
    rng = Xoshiro(3242234)
    tw = Petri{Int,Float64}()
    # Show that distributions with very different rates
    # are all sampled equally by Petri.
    enable!(tw, 3, Exponential(100.0), 0.0, 0.0)
    @test length(tw) == 1 && 3 ∈ keys(tw)
    enable!(tw, 4, Exponential(0.001), 0.0, 3.0)
    @test length(tw) == 2 && 4 ∈ keys(tw)
    enable!(tw, 7, Exponential(0.001), 5.0, 5.0)
    @test length(tw) == 3 && 7 ∈ keys(tw)

    counts = Dict{Int,Int}(3 => 0, 4 => 0, 7 => 0)
    for i in 1:1000
        when = 100 * rand(rng)
        time_out, which = next(tw, when)
        counts[which] += 1
        @assert abs(time_out - when - 1) < 1e-9
    end
    hi, lo = (maximum(values(counts)), minimum(values(counts)))
    @test (hi - lo) / lo < 0.3

    disable!(tw, 4, 9.0)
    @test length(tw) == 2 && 4 ∉ keys(tw)

    dst = Petri{Int,Float64}()
    enable!(dst, 11, Exponential(), 5.0, 5.0)
    copy_clocks!(dst, tw)
    @test length(tw) == 2 && 11 ∉ keys(tw)
end


@safetestset track_Petri_more = "Petri more functions" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: Petri
    using Random
    using Base
    rng = Xoshiro(3242234)
    tw = Petri{Int,Float64}()
    # Show that distributions with very different rates
    # are all sampled equally by Petri.
    enable!(tw, 3, Exponential(100.0), 0.0, 0.0)
    @test tw[3].clock == 3
    @test tw[3].distribution == Exponential(100.0)
    @test haskey(tw, 3)
    @test !haskey(tw, 4)
end


@safetestset track_Petri_clone = "Petri clone" begin
    using Distributions: Exponential
    using CompetingClocks: Petri, enable!, similar_sampler
    using Random: Xoshiro

    rng = Xoshiro(234567)
    sampler = Petri{Int,Float64}(2.5)  # custom time_duration
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0)

    cloned = similar_sampler(sampler)
    @test length(cloned) == 0  # cloned is empty
    @test cloned.time_duration == 2.5  # time_duration is preserved
    @test length(sampler) == 2  # original unchanged
end


@safetestset track_Petri_fire = "Petri fire! disables the fired clock" begin
    using CompetingClocks: Petri, enable!, fire!, next, isenabled
    using Random: Xoshiro
    using Distributions: Exponential

    # Petri retains no residual draw randomness, so fire! is the interface
    # fallback and acts as disable!. This pins that the fallback reaches it.
    rng = Xoshiro(979697)
    petri = Petri{Int,Float64}(1.0)
    enable!(petri, 1, Exponential(1.0), 0.0, 0.0)
    enable!(petri, 2, Exponential(1.0), 0.0, 0.0)
    when, which = next(petri, 0.0)
    fire!(petri, which, when)
    @test !isenabled(petri, which)
    @test isenabled(petri, which == 1 ? 2 : 1)
    @test length(petri) == 1
end
