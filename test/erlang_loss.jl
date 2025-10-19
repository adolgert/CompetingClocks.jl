module ErlangLoss
using Random
using Distributions: Exponential
using CompetingClocks

export ErlangLoss, BasicErlangLoss, ErlangLossObserver, observe, step!, run_erlangloss


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

function step!(model::ErlangLoss, sampler, which, when, rng)
    # We assume fire!() has been called on `which`.
    if which == (:arrival, 0)
        if length(model.busy) < model.c
            server_id = rand(rng, setdiff(1:model.c, model.busy))
            enable!(sampler, (:server, server_id), Exponential(1 / model.μ), when, when, rng)
            push!(model.busy, server_id)
            # else let the arrival go.
        end
        enable!(sampler, (:arrival, 0), Exponential(1 / model.λ), when, when, rng)
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
    enable!(sampler, (:arrival, 0), Exponential(1 / model.λ), when, when, rng)
    for i in 1:100
        (when, which) = next(sampler)
        fire!(sampler, which, when)
        step!(model, sampler, when, which, rng)
        observe(observer, model, which, when)
    end
    return observer
end

end

module NonErlangLoss
using Random
using Distributions: Exponential
using CompetingClocks

export ErlangLoss, BasicErlangLoss, ErlangLossObserver, observe, step!, run_erlangloss


"""
This NonErlangLoss gives each server a different rate. We do this so that we
can compare across Direct algorithms whether they give the exact same answers.

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
    slope::Float64
    busy::Set{Int} # busy servers
end


function BasicErlangLoss()
    ErlangLoss(10, 8.0, 1.0, 0.01, Set{Int}())
end

function step!(model::ErlangLoss, samplers, which, when, rng, server_rng)
    # We assume fire!() has been called on `which`.
    if which == (:arrival, 0)
        if length(model.busy) < model.c
            # Using server_rng to make this deterministic.
            available_servers = collect(sort(setdiff(1:model.c, model.busy)))
            server_id = available_servers[1+mod(server_rng, length(available_servers))]
            rate = model.μ + model.slope * server_id
            for sampler in samplers
                enable!(sampler, (:server, server_id), Exponential(1 / rate), when, when, rng)
            end
            push!(model.busy, server_id)
            # else let the arrival go.
        end
        for sampler in samplers
            enable!(sampler, (:arrival, 0), Exponential(1 / model.λ), when, when, rng)
        end
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
    enable!(sampler, (:arrival, 0), Exponential(1 / model.λ), when, when, rng)
    for i in 1:100
        (when, which) = next(sampler)
        fire!(sampler, which, when)
        randint = rand(rng, Int)
        step!(model, sampler, when, which, rng, randint)
        observe(observer, model, which, when)
    end
    return observer
end

end
