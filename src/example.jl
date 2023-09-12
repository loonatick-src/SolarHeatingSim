using UnPack, Plots
using OrdinaryDiffEq

function run_example(;tspan=(0.0, 1000.0), dt = 1.0e-3, time_values=LinRange(0, 50_000, 51))
    model = SWHSModel(;ρe=1.0e3)
    prob = SWHSProblem(;Tf0=280.0, Tp0=300.0)
    tl = @view time_values[begin:end-1]
    tr = @view time_values[begin+1:end]
    timeseries_data = Vector{Vector{value_type(prob)}}()
    for tspan in zip(tl,tr)
        push!(runge_kutta_4(f!, tspan, dt, prob))
    end
    (time_values, timeseries_data)
end  # module

ode_f!(du, u, cache, t) = f!(du, u, cache)

function using_DifferentialEquations_jl(;tspan=(0.0,10.0), dt=1.0e-2)
    model =SWHSModel()
    _prob = SWHSProblem(model)
    u0 = _prob.u
    cache = _prob.cache
    prob = ODEProblem(ode_f!, u0, tspan, cache)
    sol = solve(prob, Tsit5(); saveat=(last(tspan)-first(tspan))/10.0)
end
