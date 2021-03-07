using Fleck

@safetestset never_fair = "Never fair weather" begin
    using Fleck: Never

    n = Never()
    p = params(n)
    @test p == ()
    @test isa(n, Never)
end
