using Random: AbstractRNG

"""
    hazards(::Function, process, ::AbstractRNG)

A process can pass transition distributions by creating a `hazards`
function that follows a visitor pattern. The first argument is a
callback function that takes three arguments: an identifier for a
clock, a distribution, and whether to enable or disable it.
"""
function hazards(::Function, process, ::AbstractRNG)
end

export hazards
