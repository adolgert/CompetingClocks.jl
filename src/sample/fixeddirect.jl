
"""
This Direct method assumes that there are a fixed total number
of transitions and that each clock is marked with :Index,
numbered from 1, as an extra argument to the AddTransition!.
"""
struct FixedDirect
    tree::PrefixSearchTree{Float64}
    N::Int
    clock_index::Dict{Int, Any}
    init::Bool
    FixedDirect(N::Int) = new(
        PrefixSearchTree(Float64, N), N,
        Dict{Int,Any}(), true
        )
end


function fd_indexof(kind::Array{Any,1})
    for (symbol, value) in kind
        if symbol==:index
            return value
        end
    end
    error("Need an index=<int> for AddTransition!")
end


function Next(propagator::FixedDirect, process, rng)
    if propagator.init
        hazards=Array{Tuple{Int, Float64}, 1}()
        Hazards(process, rng) do clock, now, enabled, rng2
            lambda=Parameters(Distribution(clock.intensity))[1]
            index=fd_indexof(clock.kind)
            propagator.clock_index[index]=clock
            push!(hazards, (index, lambda))
        end
        push!(propagator.tree, hazards)
        propagator.init=false
    end
    total=sum(propagator.tree)
    if total>eps(Float64)
        (index, value)=choose(propagator.tree, rand(rng)*total)
        clock=propagator.clock_index[index]
        return (Time(process)-log(rand(rng))/total, clock)
    else
        return (Inf, nothing)
    end
end


function Observer(propagator::FixedDirect)
    function fdobserve(clock, time, updated, rng)
        if updated != :Disabled && updated != :Fired
            index = fd_indexof(clock.kind)
            propagator.clock_index[index] = clock
            lambda = Parameters(Distribution(clock.intensity))[1]
            push!(propagator.tree, index, lambda)
        else
            index=fd_indexof(clock.kind)
            push!(propagator.tree, index, 0.0)
        end
    end
end
