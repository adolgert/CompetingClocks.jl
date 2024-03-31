using SimpleDiffEq, SimpleNonlinearSolve, Distributions

# gene regulation ex from https://github.com/CharlotteJana/pdmpsim?tab=readme-ov-file#a-simple-example
# use algo from https://arxiv.org/abs/1504.06873
# especially the change of time variable, eqn 3.1 to avoid costly inversion

const pars = (b = 0.5, a0 = 1, a1 = 3, k10 = 1, k01 = 0.5)

function simple_ode_analytic(a,b,x0,t)
    C = x0 - a/b
    xt = a/b + C*exp(-b*t)
    return xt
end

function simple_ode_numeric(a,b,k,x0,t)
    function dfdt!(du, u, p, t)
        f, Λ = u
        du[1] = a - b*f # the smooth dynamics
        du[2] = k*f # the rate function
        return nothing
    end
    odeprob = ODEProblem(dfdt!, [x0,0.0], (0.0,t), [])
    sol = solve(odeprob, SimpleRK4(), dt=1e-4)
    return sol
end

x0 = 5.0

sol1 = simple_ode_analytic(pars.a0, pars.b, x0, 20.0)
sol2 = simple_ode_numeric(pars.a0, pars.b, pars.k01, x0, 8.0)

sol1 ≈ only(sol2[end])

# the above just gives the flow, to sample the time of the next jump, we need to plug into
# the rate function and solve the cumulative hazard 
u = rand(Exponential(1))

function solve_jump_time_analytic(p)
    f(u,p) = p.u - simple_ode_analytic(p.a, p.b, p.x0, only(u))
    nlprob = NonlinearProblem(f, [1.0])
    sol = solve(nlprob, SimpleNewtonRaphson(), abstol = 1e-9)
    return sol
end



p = NamedTuple{(keys(pars)..., :x0)}((pars..., x0))
sol = solve_jump_time_analytic()

mutable struct GeneRegulation{T<:NamedTuple}
    parameters::T
    f::Float64
    d::Int
    when::Float64
end

# examplePDMP <- new("pdmpModel",
#                   descr = "Gene regulation with positive feedback",
#                   parms = list(b = 0.5, a0 = 1, a1 = 3, k10 = 1, k01 = 0.5), 
#                   init = c(f = 1, d = 1),
#                   discStates = list(d = 0:1),
#                   dynfunc = function(t, x, parms) {
#                     df <- with(as.list(c(x, parms)), {
#                       switch(d+1, a0 - b*f, a1 - b*f)
#                     })
#                     return(c(df, 0))
#                   }, 
#                   ratefunc = function(t, x, parms) {
#                     return(with(as.list(c(x, parms)), switch(d + 1, k01*f, k10)))
#                   }, 
#                   jumpfunc = function(t, x, parms, jtype) {
#                     c(x[1], 1 - x[2])
#                   }, 
#                   times = c(from = 0, to = 100, by = 0.1), 
#                   solver = "lsodar")