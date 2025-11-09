@safetestset available_samplers_test = "SamplerSpec list all" begin
    using CompetingClocks
    samplist = available_samplers()
    tofind = ["NextReaction", "FirstReaction", "Direct", "Petri", "FirstToFire"]
    for sstring in tofind
        @test any(occursin(sstring, docstring) for docstring in samplist)
    end
end


@safetestset samplerspec_nextreaction = "SamplerSpec next reaction" begin
    using CompetingClocks

    spec = NextReactionMethod()
    sampler = spec(Int, Float64)
    @test isa(sampler, SSA)
end


@safetestset samplerspec_direct = "SamplerSpec direct method" begin
    using CompetingClocks

    spec = DirectMethod()
    sampler = spec(Int, Float64)
    @test isa(sampler, SSA)
    
    spec2 = DirectMethod(:keep, :array)
    sampler2 = spec2(Int, Float64)
    @test isa(sampler2, SSA)

    spec3 = DirectMethod(:remove)
    sampler3 = spec3(Int, Float64)
    @test isa(sampler3, SSA)

    spec4 = DirectMethod(:tree)
    sampler4 = spec4(Int, Float64)
    @test isa(sampler4, SSA)

    @test_throws "Specify one of" DirectMethod(:funnybusiness)
    @test_throws "Specify one of" DirectMethod(:keep, :keep)
    @test_throws "Specify one of" DirectMethod(:tree, :array)
end


@safetestset samplerspec_firstreaction = "SamplerSpec first reaction" begin
    using CompetingClocks

    spec = FirstReactionMethod()
    sampler = spec(Int, Float64)
    @test isa(sampler, SSA)
end


@safetestset samplerspec_firsttofire = "SamplerSpec first to fire" begin
    using CompetingClocks

    spec = FirstToFireMethod()
    sampler = spec(Int, Float64)
    @test isa(sampler, SSA)
end


@safetestset samplerspec_petri = "SamplerSpec petri method" begin
    using CompetingClocks

    spec = PetriMethod()
    sampler = spec(Int, Float64)
    @test isa(sampler, SSA)
    
    # Test with custom dt
    spec2 = PetriMethod(0.5)
    sampler2 = spec2(Int, Float64)
    @test isa(sampler2, SSA)
end
