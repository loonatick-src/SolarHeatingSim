module SimluationViz

using Simulation, UnPack, Plots

time_values, timeseries_data = Simulation.run_example(;tspan=(0.0, 100.0), dt = 1.0e-3, saveat=1.0)

end # module