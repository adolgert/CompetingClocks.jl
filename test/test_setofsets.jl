using SafeTestsets

@safetestset setofsets_smoke = "SetOfSets smoke" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    c = SetOfSets([a, b])
    @test c isa AbstractSet{Int64}
    @test 4 in c
    @test 3 in c
    @test 24 ∉ c
end

@safetestset setofsets_empty = "SetOfSets empty" begin
    using CompetingClocks: SetOfSets

    # Empty vector of sets
    empty_sos = SetOfSets(Set{Int}[])
    @test isempty(collect(empty_sos))
    @test length(empty_sos) == 0
    @test !(1 in empty_sos)

    # Vector with empty sets
    a = Set{Int}()
    b = Set{Int}()
    sos = SetOfSets([a, b])
    @test isempty(collect(sos))
    @test length(sos) == 0
end

@safetestset setofsets_iteration = "SetOfSets iteration" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    c = Set([5, 8])
    sos = SetOfSets([a, b, c])

    # Collect all elements
    elements = collect(sos)
    @test length(elements) == 8  # Fixed: 8 elements total
    @test Set(elements) == Set([1, 2, 3, 4, 5, 6, 7, 8])

    # Test that iteration works multiple times
    elements2 = collect(sos)
    @test elements2 == elements || Set(elements2) == Set(elements)
end

@safetestset setofsets_length = "SetOfSets length" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    c = Set([5])
    sos = SetOfSets([a, b, c])

    @test length(sos) == 7
    @test length(a) + length(b) + length(c) == length(sos)
end

@safetestset setofsets_membership = "SetOfSets membership" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    c = Set([5, 8, 9])
    sos = SetOfSets([a, b, c])

    # Test all elements are found
    for x in 1:9
        @test x in sos
    end

    # Test non-members
    @test !(10 in sos)
    @test !(0 in sos)
    @test !(-1 in sos)
end

@safetestset setofsets_eltype = "SetOfSets eltype" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    sos = SetOfSets([a, b])

    @test eltype(sos) == Int
    @test eltype(typeof(sos)) == Int
end

@safetestset setofsets_union = "SetOfSets union" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    sos = SetOfSets([a, b])

    # Union with another set
    other = Set([8, 9, 10])
    result = union(sos, other)
    @test result == Set([1, 2, 3, 4, 6, 7, 8, 9, 10])

    # Union with multiple sets
    other2 = Set([11, 12])
    result2 = union(sos, other, other2)
    @test result2 == Set([1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12])
end

@safetestset setofsets_intersect = "SetOfSets intersect" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7, 10, 15])
    b = Set([2, 4, 6, 10, 15])
    sos = SetOfSets([a, b])

    # Intersect with another set
    other = Set([10, 15, 20])
    result = intersect(sos, other)
    @test result == Set([10, 15])
end

@safetestset setofsets_setdiff = "SetOfSets setdiff" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7])
    b = Set([2, 4, 6])
    sos = SetOfSets([a, b])

    # Set difference: setdiff splatted means setdiff(a, b, other)
    # which removes elements of b and other from a
    other = Set([7])
    result = setdiff(sos, other)
    # setdiff(a, b, other) = setdiff(Set([1,3,7]), Set([2,4,6]), Set([7]))
    # = Set([1,3])
    @test result == Set([1, 3])
end

@safetestset setofsets_issubset = "SetOfSets issubset" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3])
    b = Set([2, 4])
    sos = SetOfSets([a, b])

    # Test subset relationships
    superset = Set([1, 2, 3, 4, 5, 6])
    @test issubset(sos, superset)

    not_superset = Set([1, 2, 3])
    @test !issubset(sos, not_superset)

    exact_set = Set([1, 2, 3, 4])
    @test issubset(sos, exact_set)
end

@safetestset setofsets_skip_empty = "SetOfSets skip empty subsets" begin
    using CompetingClocks: SetOfSets

    # Mix of empty and non-empty sets
    a = Set([1, 3])
    b = Set{Int}()
    c = Set([5, 7])
    d = Set{Int}()
    e = Set([9])
    sos = SetOfSets([a, b, c, d, e])

    elements = collect(sos)
    @test Set(elements) == Set([1, 3, 5, 7, 9])
    @test length(sos) == 5
end

@safetestset setofsets_single_set = "SetOfSets single set" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 3, 7, 9, 11])
    sos = SetOfSets([a])

    @test length(sos) == 5
    @test 7 in sos
    @test !(8 in sos)
    @test collect(sos) |> Set == a
end

@safetestset setofsets_string_type = "SetOfSets with strings" begin
    using CompetingClocks: SetOfSets

    a = Set(["hello", "world"])
    b = Set(["foo", "bar"])
    sos = SetOfSets([a, b])

    @test eltype(sos) == String
    @test "hello" in sos
    @test "bar" in sos
    @test !("baz" in sos)
    @test length(sos) == 4
end

@safetestset setofsets_union_two_sos = "SetOfSets union of two SetOfSets" begin
    using CompetingClocks: SetOfSets

    a = Set([1, 2])
    b = Set([3, 4])
    sos1 = SetOfSets([a, b])

    c = Set([5, 6])
    d = Set([7, 8])
    sos2 = SetOfSets([c, d])

    # Union of two SetOfSets
    result = union(sos1, sos2)
    @test result isa SetOfSets
    @test Set(collect(result)) == Set(1:8)
    @test length(result) == 8
end

@safetestset setofsets_isempty_method = "SetOfSets isempty" begin
    using CompetingClocks: SetOfSets

    # Non-empty SetOfSets
    a = Set([1, 2])
    b = Set([3, 4])
    sos = SetOfSets([a, b])
    @test !isempty(sos)

    # Empty SetOfSets (no subsets)
    empty_sos = SetOfSets(Set{Int}[])
    @test isempty(empty_sos)

    # SetOfSets with only empty subsets
    empty_a = Set{Int}()
    empty_b = Set{Int}()
    sos_empty_subsets = SetOfSets([empty_a, empty_b])
    @test isempty(sos_empty_subsets)
end
