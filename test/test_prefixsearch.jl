using SafeTestsets


@safetestset cumsumprefix_smoke = "cumsum prefix smoke" begin
    using Fleck: CumSumPrefixSearch, choose
    ps = CumSumPrefixSearch(Float64)
    for idx in 1:20
        push!(ps, 0.5)
    end
    @test length(ps) == 20
    maximum = sum!(ps)
    @test abs(maximum - 20 * 0.5) < 1e-13
    for check in 1:200
        variate = (check - 1) * maximum / 200
        idx = choose(ps, variate)
        @test idx >= 1
        @test idx <= 20
    end
end


@safetestset cumsumprefix_rand = "cumsum prefix has rand" begin
    using Fleck: CumSumPrefixSearch, choose
    using Random
    rng = Xoshiro(23423)
    ps = CumSumPrefixSearch(Float64)
    for idx in 1:20
        push!(ps, 0.5)
    end
    for i in 1:100
        idx = rand(rng, ps)
        @test idx >= 1
        @test idx <= 20
    end
end


@safetestset prefixsearch_single = "test single value" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
	t = BinaryTreePrefixSearch([3])
	@test sum(t) == 3
	@test choose(t, 2)[1] == 1
	push!(t, [(1, 4)])
	@test sum(t) == 4
	@test choose(t, 3.3)[1] == 1
end


@safetestset prefixsearch_double = "test double value" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
    t = BinaryTreePrefixSearch([3, 1])
	@test sum(t) == 4
	@test choose(t, 2)[1] == 1
	@test choose(t, 3)[1] == 2
	@test choose(t, 3.7)[1] == 2
	push!(t, [(1, 4)])
	v = [(2, 1), (3.3, 1), (4.1, 2)]
	for (guess, result) in v
		@test choose(t, guess)[1] == result
	end
	push!(t, [(2, 2)])
	v = [(2, 1), (3.3, 1), (5.1, 2)]
	for (guess, result) in v
		@test choose(t, guess)[1] == result
	end
end


@safetestset prefixsearch_three = "test three values" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
    t = BinaryTreePrefixSearch([3, 1, 2])
    @test sum(t) == 6
    @test choose(t, 2.1)[1] == 1
    @test choose(t, 3)[1] == 2
    @test choose(t, 3.1)[1] == 2
    @test choose(t, 5)[1] == 3
    push!(t, ((2, 2), (3, 3)))
    v = [(2, 1), (3, 2), (4.8, 2), (5.1, 3), (7.9, 3)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end

@safetestset prefixsearch_three_floats = "three floats" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
    t = BinaryTreePrefixSearch([3.5, 1.5, 2.5])
    @test abs(sum(t) - 7.5) < 1e-9
    @test choose(t, 2.1)[1] == 1
    @test choose(t, 3.5)[1] == 2
    @test choose(t, 3.7)[1] == 2
    @test choose(t, 5.5)[1] == 3
    push!(t, ((2, 2.5), (3, 3.5)))
    v = [(2, 1),(3.5, 2),(4.8, 2),(6.1, 3),(7.9, 3)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_four_floats = "four floats" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
    vals = [3, 1, 2, 4]
    t = BinaryTreePrefixSearch(vals)
    @test sum(t) == 10
    v = [(2, 1), (3, 2),(4.1, 3), (9.9, 4)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_five = "five values" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum
    vals = [3, 0, 2, 4, 1]
    t = BinaryTreePrefixSearch(vals)
    @test sum(t) == 10
    v = [(2, 1), (3, 3), (4.1, 3), (9.9, 5)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
    push!(t, [(2, 7), (3, 3)])
    @test sum(t) == 18
    v = [(17.7, 5)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end
