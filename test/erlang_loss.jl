module ErlangLoss
using Random
using Distributions: Exponential

"""
An Erlang loss system, also known as an M/M/c/c queue, is
a queuing model where there is no queue for arriving customers. If all c servers are busy, any new customer is "lost" or "blocked" and leaves the system immediately. The model's name comes from its use of the Erlang B formula to calculate the probability of loss. 
Key characteristics

 * Arrivals: Customers arrive according to a Poisson process.
 * Service: Service times are exponentially distributed.
 * Servers: There are c identical servers.
 * Queue: There is no waiting space, which is represented by the second c in the M/M/c/c notation.
 * Customer behavior: If a customer arrives and all c servers are busy, the customer is blocked and does not enter the system. 
 
"""
mutable struct ErlangLoss
    c::Int  # number of servers
    λ::Float64 # arrival rate
    μ::Float64 # service rate per servers
    busy::Set{Int} # busy servers
end


function BasicErlangLoss()
    ErlangLoss(10, 8.0, 1.0, Set{Int}())
end

function step!(model, sampler, when, which, rng)
    # We assume fire!() has been called on `which`.
    if which == (:arrival, 0)
        if length(model.busy) < model.c
            server_id = rand(rng, setdiff(1:model.c, model.busy))
            enable!(sampler, (:server, server_id), Exponential(1 / model.μ), when)
            push!(model.busy, sever_id)
            # else let the arrival go.
        end
        enable!(sampler, (:arrival, 0), Exponential(1 / model.λ, when))
    else
        delete!(model.busy, which[2])
    end
end


mutable struct ErlangLossObserver
    busy::Vector{Float64}
    time::Float64
    missed::Int
    arrivals::Int
    ErlangLossObserver(el::ErlangLoss) = new(zeros(Float64, el.c), zero(Float64), 0, 0)
end


function observe(elo::ErlangLossObserver, model, which, when)
    N = length(model.busy)
    if which == (:arrival, 0)
        elo.arrivals += 1
        if N == model.c
            elo.missed += 1
        end
    end
    duration = when - elo.time
    for i in 1:N
        elo.busy[i] += duration
    end
    elo.time = when
end

function run_erlangloss(sampler, rng)
    model = BasicErlangLoss()
    observer = ErlangLossObserver(model)
    for i in 1:100
        (when, which) = next(sampler)
        fire!(sampler, which, when)
        step!(model, sampler, when, which, rng)
        observe(observer, model, which, when)
    end
    return observer
end

end
