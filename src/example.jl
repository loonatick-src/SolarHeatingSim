using UnPack, Plots, LaTeXStrings

function run_example(;tspan=(0.0, 100.0), dt = 1.0e-3, saveat=1.0)
    model = SWHSModel(;ρe=1.0e3)
    prob = SWHSProblem(model;Tf0=280.0, Tp0=300.0)
    timeseries_data = Vector{Vector{value_type(model)}}()
    t0 = first(tspan)
    tN = last(tspan)
    time_values = typeof(dt)[]
    push!(time_values, t0)
    push!(timeseries_data, copy(prob.u))
    for ti in (t0+saveat):saveat:tN
        push!(timeseries_data, copy(runge_kutta_4(f!, (t0,ti), dt, prob)))
        push!(time_values, ti)
    end
    (time_values, timeseries_data, prob)
end  # module

# function using_DifferentialEquations_jl(;tspan=(0.0,10.0), dt=1.0e-2)
#     model =SWHSModel()
#     _prob = SWHSProblem(model)
#     u0 = _prob.u
#     cache = _prob.cache
#     prob = ODEProblem(sciml_f!, u0, tspan, cache)
#     solve(prob, Tsit5(); saveat=(last(tspan)-first(tspan))/10.0)
# end
