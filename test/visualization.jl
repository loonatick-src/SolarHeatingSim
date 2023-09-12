using Plots, LaTeXStrings
using LaTeXStrings

plotsdir() = "plots"

flatten_time(time_values, data) = [t for t in time_values  for _ in eachindex(data[1]) ]
flatten_space(y_values, data) = [y  for _ in eachindex(data) for y in y_values]

@testset "Plotting" begin
time_values, timeseries_data, prob = Simulation.run_example(;tspan=(0.0, 100.0), dt = 1.0e-3, saveat=1.0)

function plot_surface(time_values, timeseries_data, prob::SWHSProblem)
    @unpack L, Nc, Ns = prob.cache.model
    y_values = LinRange(0, L, Nc)
    strats = 1:Ns
    plate_data = map(timeseries_data) do data
        plate_view(data, prob)
    end
    plate_data_flat = vcat(plate_data...)
    fluid_data = map(timeseries_data) do data
        fluid_view(data, prob)
    end
    fluid_data_flat = vcat(fluid_data...)
    tank_data = map(timeseries_data) do data
         tank_view(data, prob)
    end
    tank_data_flat = vcat(tank_data...)
    tp = flatten_time(time_values, plate_data)
    tf = flatten_time(time_values, fluid_data)
    @assert length(tp) == length(tf)
    y = flatten_space(y_values, plate_data)
    @assert length(tf) == length(y)
    tt = flatten_time(time_values, tank_data)
    nodes = flatten_space(strats, tank_data)
    @assert length(tt) == length(nodes)
    p_fluid = surface(tf, y, fluid_data_flat; xlabel="time (s)", ylabel="y (m)",
                        zlabel="Temperature (K)", title="Fluid in collector")
    p_plate = surface(tp, y, plate_data_flat; xlabel="time (s)", ylabel="y (m)",
                        zlabel="Temperature (K)", title="Collector Plate")
    p_tank = surface(tt, nodes, tank_data_flat, xlabel="time (s)", ylabel="Stratification layer",
                        zlabel="Temperature (K)", title="Storage Tank")
    p_fluid,p_plate,p_tank
end

@info "Generating plots"
p_fluid, p_plate, p_tank = plot_surface(time_values, timeseries_data, prob)

@info "saving plots"
savefig(p_fluid, joinpath(plotsdir(), "collector_fluid.png"))
savefig(p_plate, joinpath(plotsdir(), "plate.png"))
savefig(p_tank, joinpath(plotsdir(), "tank.png"))
@info "Done. Plots saved in \"/plots\""
end  # testset