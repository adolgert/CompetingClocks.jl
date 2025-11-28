using SafeTestsets

@safetestset prefixsearch_log2 = "log2 internal" begin
    using CompetingClocks: ceil_log2

    # I'm replacing one clear expression with one that's more obscure but fast.
    # This tests that relationship.
    for i in 1:500
        @test ceil_log2(i) == Int(ceil(log2(i)))
    end
end


@safetestset cumsumprefix_smoke = "cumsum prefix smoke" begin
    using CompetingClocks: CumSumPrefixSearch, choose
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
    using CompetingClocks: CumSumPrefixSearch, choose
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

@safetestset cumsumprefix_copy = "cumsum prefix copy!" begin
    using CompetingClocks: CumSumPrefixSearch

    src = CumSumPrefixSearch{Float64}()
    for v in [1.0, 2.0, 3.0]
        push!(src, v)
    end
    sum!(src)  # populate cumulant

    dst = CumSumPrefixSearch{Float64}()
    for v in [0.0, 0.0, 0.0]
        push!(dst, v)
    end

    copy!(dst, src)

    @test dst.array == src.array
    @test dst.cumulant == src.cumulant
    @test dst.dirty == src.dirty
    @test length(dst) == length(src)
end

@safetestset cumsumprefix_time_type = "cumsum prefix time_type instance" begin
    using CompetingClocks: CumSumPrefixSearch, time_type

    ps_float = CumSumPrefixSearch{Float64}()
    @test time_type(ps_float) == Float64

    ps_float32 = CumSumPrefixSearch{Float32}()
    @test time_type(ps_float32) == Float32
end

@safetestset cumsumprefix_isenabled = "cumsum prefix isenabled" begin
    using CompetingClocks: CumSumPrefixSearch, isenabled

    ps = CumSumPrefixSearch{Float64}()
    push!(ps, 1.0)
    push!(ps, 0.0)  # zero value
    push!(ps, 3.0)

    # Enabled clocks (positive value)
    @test isenabled(ps, 1) == true
    @test isenabled(ps, 3) == true

    # Disabled clock (zero value)
    @test isenabled(ps, 2) == false

    # Out of bounds
    @test isenabled(ps, 0) == false
    @test isenabled(ps, 4) == false
    @test isenabled(ps, -1) == false
end

@safetestset prefixsearch_structure = "test structure" begin
    using CompetingClocks: _btps_sizes
    # depth offset array_cnt
    @test _btps_sizes(1) == (1, 1, 1)
    @test _btps_sizes(2) == (2, 2, 3)
    @test _btps_sizes(3) == (3, 4, 7) # values in positions [4, 5, 6, 7]
    @test _btps_sizes(4) == _btps_sizes(3)
    @test _btps_sizes(5) == (4, 8, 15)
end


@safetestset prefixsearch_single = "test single value" begin
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
	t = BinaryTreePrefixSearch{Int64}()
    push!(t, 3)
	@test sum!(t) == 3
	@test choose(t, 2)[1] == 1
    t[1] = 4
	@test sum!(t) == 4
	@test choose(t, 3.3)[1] == 1
end


@safetestset prefixsearch_add_values = "binarytreeprefix add values" begin
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!, allocated
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
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
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
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
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
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
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
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
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
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!
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
    using CompetingClocks: BinaryTreePrefixSearch, set_multiple!
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


@safetestset prefixsearch_empty = "testing empty!" begin
    using Random
    using CompetingClocks: BinaryTreePrefixSearch, set_multiple!, allocated

    rng = Xoshiro(3297920)
    cnt = 20
    btps1 = BinaryTreePrefixSearch{Int64}(cnt)
    # The initial allocation is 32.
    @test allocated(btps1) == Int(2^ceil(log2(cnt-1)))
    vals = rand(rng, 0:500, 2 * cnt)
    set_multiple!(btps1, collect(enumerate(vals)))
    @test length(btps1) == 2 * cnt
    # That allocation has to grow.
    @test allocated(btps1) == Int(2^ceil(log2(2 * cnt-1)))
    empty!(btps1)
    # After emptying, the allocaiton returns to 32.
    @test length(btps1) == 0
    @test allocated(btps1) == Int(2^ceil(log2(cnt-1)))
end

@safetestset prefixsearch_time_type = "time_type instance method" begin
    using CompetingClocks: BinaryTreePrefixSearch, time_type

    btps_float = BinaryTreePrefixSearch{Float64}(4)
    @test time_type(btps_float) == Float64

    btps_int = BinaryTreePrefixSearch{Int64}(4)
    @test time_type(btps_int) == Int64
end

@safetestset prefixsearch_choose_error = "choose errors on invalid value" begin
    using CompetingClocks: BinaryTreePrefixSearch, choose, sum!

    btps = BinaryTreePrefixSearch{Float64}(4)
    push!(btps, 1.0)
    push!(btps, 2.0)
    push!(btps, 3.0)

    total = sum!(btps)
    @test total == 6.0

    # Value equal to total should error
    @test_throws ErrorException choose(btps, 6.0)

    # Value greater than total should error
    @test_throws ErrorException choose(btps, 7.0)
end

@safetestset prefixsearch_haskey = "haskey methods" begin
    using CompetingClocks: BinaryTreePrefixSearch

    btps = BinaryTreePrefixSearch{Float64}(4)
    push!(btps, 1.0)
    push!(btps, 0.0)  # zero value
    push!(btps, 3.0)

    # haskey with Int clock - positive value exists
    @test haskey(btps, 1) == true
    @test haskey(btps, 3) == true

    # haskey with Int clock - zero value does not count as "has"
    @test haskey(btps, 2) == false

    # haskey with Int clock - out of bounds
    @test haskey(btps, 0) == false
    @test haskey(btps, 4) == false
    @test haskey(btps, -1) == false

    # haskey works with any Integer type
    @test haskey(btps, Int32(1)) == true
    @test haskey(btps, UInt(3)) == true
end

@safetestset prefixsearch_rand_sampler = "rand with SamplerTrivial" begin
    using Random
    using CompetingClocks: BinaryTreePrefixSearch

    btps = BinaryTreePrefixSearch{Float64}(4)
    push!(btps, 1.0)
    push!(btps, 2.0)
    push!(btps, 3.0)
    push!(btps, 4.0)

    rng = Xoshiro(12345)

    # Sample multiple times and verify results are valid
    for _ in 1:100
        idx, remainder = rand(rng, btps)
        @test 1 <= idx <= 4
        @test remainder >= 0.0
    end
end
