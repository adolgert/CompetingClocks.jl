using SafeTestsets

@safetestset prefixsearch_log2 = "log2 internal" begin
    using Fleck: ceil_log2

    # I'm replacing one clear expression with one that's more obscure but fast.
    # This tests that relationship.
    for i in 1:500
        @test ceil_log2(i) == Int(ceil(log2(i)))
    end
end


@safetestset cumsumprefix_smoke = "cumsum prefix smoke" begin
    using Fleck: CumSumPrefixSearch, choose
    ps = CumSumPrefixSearch{Float64}()
    for idx in 1:20
        push!(ps, 0.5)
    end
    @test length(ps) == 20
    maximum = sum!(ps)
    @test abs(maximum - 20 * 0.5) < 1e-13
    for check in 1:200
        variate = (check - 1) * maximum / 200
        idx, value = choose(ps, variate)
        @test idx >= 1
        @test idx <= 20
    end
end


@safetestset cumsumprefix_rand = "cumsum prefix has rand" begin
    using Fleck: CumSumPrefixSearch, choose
    using Random
    rng = Xoshiro(23423)
    ps = CumSumPrefixSearch{Float64}()
    for idx in 1:20
        push!(ps, 0.5)
    end
    for i in 1:100
        idx, value = rand(rng, ps)
        @test idx >= 1
        @test idx <= 20
    end
end


@safetestset prefixsearch_structure = "test structure" begin
    using Fleck: _btps_sizes
    # depth offset array_cnt
    @test _btps_sizes(1) == (1, 1, 1)
    @test _btps_sizes(2) == (2, 2, 3)
    @test _btps_sizes(3) == (3, 4, 7) # values in positions [4, 5, 6, 7]
    @test _btps_sizes(4) == _btps_sizes(3)
    @test _btps_sizes(5) == (4, 8, 15)
end


@safetestset prefixsearch_single = "test single value" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
	t = BinaryTreePrefixSearch{Int64}()
    push!(t, 3)
	@test sum!(t) == 3
	@test choose(t, 2)[1] == 1
    t[1] = 4
	@test sum!(t) == 4
	@test choose(t, 3.3)[1] == 1
end


@safetestset prefixsearch_add_values = "binarytreeprefix add values" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!, allocated
    initial_allocation = 1
	t = BinaryTreePrefixSearch{Int64}(initial_allocation)
    push!(t, 3)
    @test t.cnt == 1
    @test length(t) == 1
    @test allocated(t) == 1
    push!(t, 4)
    @test length(t) == 2
    @test allocated(t) == 2
    push!(t, 2)
    @test length(t) == 3
    @test allocated(t) == 4
    push!(t, 3)
    @test length(t) == 4
    @test allocated(t) == 4
	@test sum!(t) == 3+4+2+3
    @test t.array[1] == 3+4+2+3
    @test t.array[2] == 3+4
    @test t.array[3] == 2+3
    @test t.array[4] == 3
	@test choose(t, 2)[1] == 1
	@test choose(t, 3.3)[1] == 2
end


@safetestset prefixsearch_double = "test double value" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
    t = BinaryTreePrefixSearch{Int64}(2)
    push!(t, 3)
    push!(t, 1)
	@test sum!(t) == 4
	@test choose(t, 2)[1] == 1
	@test choose(t, 3)[1] == 2
	@test choose(t, 3.7)[1] == 2
    t[1] = 4
	v = [(2, 1), (3.3, 1), (4.1, 2)]
	for (guess, result) in v
		@test choose(t, guess)[1] == result
	end
    t[2] = 2
	v = [(2, 1), (3.3, 1), (5.1, 2)]
	for (guess, result) in v
		@test choose(t, guess)[1] == result
	end
end


@safetestset prefixsearch_three = "test three values" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
    t = BinaryTreePrefixSearch{Int64}(3)
    for v in [3, 1, 2]
        push!(t, v)
    end
    @test sum!(t) == 6
    @test choose(t, 2.1)[1] == 1
    @test choose(t, 3)[1] == 2
    @test choose(t, 3.1)[1] == 2
    @test choose(t, 5)[1] == 3
    t[2] = 2
    t[3] = 3
    v = [(2, 1), (3, 2), (4.8, 2), (5.1, 3), (7.9, 3)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end

@safetestset prefixsearch_three_floats = "three floats" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
    t = BinaryTreePrefixSearch{Float64}(3)
    for v in [3.5, 1.5, 2.5]
        push!(t, v)
    end
    @test abs(sum!(t) - 7.5) < 1e-9
    @test choose(t, 2.1)[1] == 1
    @test choose(t, 3.5)[1] == 2
    @test choose(t, 3.7)[1] == 2
    @test choose(t, 5.5)[1] == 3
    t[2] = 2.5
    t[3] = 3.5
    v = [(2, 1),(3.5, 2),(4.8, 2),(6.1, 3),(7.9, 3)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_four_floats = "four floats" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
    vals = [3, 1, 2, 4]
    t = BinaryTreePrefixSearch{Int64}(4)
    for v in vals
        push!(t, v)
    end
    @test sum!(t) == 10
    v = [(2, 1), (3, 2),(4.1, 3), (9.9, 4)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_five = "five values" begin
    using Fleck: BinaryTreePrefixSearch, choose, sum!, update!
    vals = [3, 0, 2, 4, 1]
    t = BinaryTreePrefixSearch{Int64}(5)
    for v in vals
        push!(t, v)
    end
    @test sum!(t) == 10
    v = [(2, 1), (3, 3), (4.1, 3), (9.9, 5)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
    t[2] = 7
    t[3] = 3
    @test sum!(t) == 18
    v = [(17.7, 5)]
    for (guess, result) in v
        @test choose(t, guess)[1] == result
    end
end


@safetestset prefixsearch_getindex = "setindex getindex" begin
    using Random
    using Fleck: BinaryTreePrefixSearch, set_multiple!
    @enum Modify Append Zero Reset

    rng = Xoshiro(9347243)
    for trialidx in 1:20
        cnt = rand(rng, 1:20)
        vals = rand(rng, 0:500, cnt)
        btps = BinaryTreePrefixSearch{Int64}(cnt)
        set_multiple!(btps, collect(enumerate(vals)))
        for mod_and_check_idx in 1:500
            modify = rand(rng, instances(Modify))
            if modify == Append
                val = rand(rng, 1:500)
                push!(vals, val)
                push!(btps, val)
            elseif modify == Zero
                idx = rand(rng, 1:length(vals))
                vals[idx] = 0
                btps[idx] = 0
            elseif modify == Reset
                idx = rand(rng, 1:length(vals))
                vals[idx] = rand(rng, 1:500)
                btps[idx] = vals[idx]
            end
            same = all(btps[i] == vals[i] for i in eachindex(vals))
            @test same
            if !same
                println("modify $modify not same $(length(vals))")
                return
            end
        end
    end
end
