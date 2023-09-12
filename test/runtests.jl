module SimulationTests

using Simulation
using Test
using UnPack
using OrdinaryDiffEq, SteadyStateDiffEq

model = SWHSModel()
_prob = SWHSProblem(model)
u0 = copy(_prob.u)
tspan = (0.0, 1.0e3)
cache = _prob.cache
prob = ODEProblem(sciml_f!, u0, tspan, cache)
dt = 0.1
nsave = 10

@testset "SciML ODE solvers on discretized system" begin
    sol = solve(prob, Tsit5())
    @test sol.retcode == SciMLBase.ReturnCode.Success
end

tspan = (0.0, 100.0)
prob = ODEProblem(sciml_f!, u0, tspan, cache)
@testset "SciML Steady state solvers on discretized system" begin
    prob_ss = SteadyStateProblem(prob)
    sol = solve(prob_ss, DynamicSS(Rodas5(autodiff=false)))
    @test sol.retcode == SciMLBase.ReturnCode.Success
end

# @testset "Fully automatic solver using ModelingToolkit.jl" begin
    # TODO: these solvers take the PDE system with boundary conditions and perform the discretization for you
    # leaving this one out for now because ModelingToolkit.jl is a big bear of a dependency    
# end

end # module