using CompetingClocks

@safetestset never_fair = "Never fair weather" begin
    using CompetingClocks: Never

    n = Never()
    p = params(n)
    @test p == ()
    @test isa(n, Never)
end
