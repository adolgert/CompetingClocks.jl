using SafeTestsets


@safetestset prefixsearch_single = "test single value" begin
    using Fleck: PrefixSearchTree, choose, sum
	t = PrefixSearchTree([3])
	@test sum(t) == 3
	@test choose(t, 2)[1] == 1
	push!(t, [(1, 4)])
	@test sum(t) == 4
	@test choose(t, 3.3)[1] == 1
end


@safetestset prefixsearch_double = "test double value" begin
    using Fleck: PrefixSearchTree, choose, sum
    t = PrefixSearchTree([3, 1])
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
    using Fleck: PrefixSearchTree, choose, sum
    t = PrefixSearchTree([3, 1, 2])
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
    using Fleck: PrefixSearchTree, choose, sum
    t = PrefixSearchTree([3.5, 1.5, 2.5])
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
    using Fleck: PrefixSearchTree, choose, sum
    vals = [3, 1, 2, 4]
    t = PrefixSearchTree(vals)
    @test sum(t) == 10
    v = [(2, 1), (3, 2),(4.1, 3), (9.9, 4)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_five = "five values" begin
    using Fleck: PrefixSearchTree, choose, sum
    vals = [3, 0, 2, 4, 1]
    t = PrefixSearchTree(vals)
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
