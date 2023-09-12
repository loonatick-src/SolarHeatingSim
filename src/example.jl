using UnPack, Plots

function run_example(;tspan=(0.0, 1000.0), dt = 1.0e-3, time_values=LinRange(0, 50_000, 51))
    model = SWHSModel(;ρe=1.0e3)
    prob = SWHSProblem(;Tf0=280.0, Tp0=300.0)
    tl = @view time_values[begin:end-1]
    tr = @view time_values[begin+1:end]
    timeseries_data = Vector{Vector{value_type(prob)}}()
    for tspan in zip(tl,tr)
        push!(RK4(f!, tspan, dt, prob))
    end
    (time_values, timeseries_data)
end  # module

function plot_example(time_values, timeseries_data)
    
end
