module Simulation

using UnPack

export SWHSModel, SWHSProblem, solve

@kwdef struct SWHSModel{T}
    ρp::T  = 8.0e3    # density of plate
    δ::T   = 0.1      # plate thickness
    W::T   = 1.0      # plate width
    L::T   = 2.0      # plate length
    Cpp::T = 450.0    # specific heat capacity of plate
    S::T   = 800.0    # radiation flux
    hpf::T = 1000.0   # heat transfer coefficient (plate and working fluid)
    hpa::T = 100.0    # heat transfer coefficient (plate and atmosphere)
    Ta::T  = 300.0    # ambient temperature
    T∞::T  = 295.0    # sky temperature
    α::T   = 5.5e-8   # radiation coefficient
    ρf::T  = 1.0e3    # working fluid density
    Af::T  = 0.6      # transverse cross section area of risers
    Cpf::T = 4.2e3    # specific heat capacity of working fluid
    mdot::T = 0.5 * Af * ρf  # mass flow rate (assuming 0.5m per s in risers
    ρe::T  = 1.5e3    # effective density (tank + wokring fluid)
    Cpe::T = 5000.0   # effective specific heat of tank + working fluid
    V::T   = 10.0     # volume of each stratified tank layer
    A::T   = 5.0      # contact area between tank and atmosphere per stratification
    hta::T = 500.0    # loss coefficient (storage tank)
    Ns::Int= 10       # stratifications count for tank
end

value_type(::SWHSModel{T}) where T = T

struct SWHSProblem{T}
    model::SWHSModel{T}
    Δy::T
    Δt::T
    Nc::Int
    c1::T
    c2::T
    c3::T
    c4::T
    c5::T
    Tf0::T
    Tp0::T
    function SWHSProblem(model::SWHSModel; Nc = 100, Δt = 0.1, Tf0 = 290.0, Tp0 = 300.0)
        T = value_type(model)
        @unpack Cpp, ρp, δ, L, Ns, mdot, V, Cpe = model
        @unpack W, hpf, ρf, Af, Cpf, ρe, hta, A = model
        Δy = L / Nc
        c1 = 1/(Cpp * ρp * δ)
        c2 = 1/(ρf * Af)
        c3 = W*hpf/Cpf
        c4 = 1/(ρe*V)
        c5 = hta * A / Cpe
        new{Float64}(model, Δy, Δt, Nc, c1, c2, c3, c4, c5, Tf0, Tp0)
    end
end

function step!(timeseries_data, prob::SWHSProblem, i)
    @unpack Δt, Δy, Nc, c1, c2, c3, c4, c5 = prob
    @unpack Ns, S, mdot, hpf, hpa, hta, Ta, T∞, α = prob.model
    Tcurr = @view timeseries_data[:,i]
    Tprev = @view timeseries_data[:,i-1]
    
    plate_begin = firstindex(Tcurr)
    plate_end = Nc
    fluid_begin = plate_end+1
    fluid_end = 2Nc
    tank_begin = fluid_end + 1
    tank_end = lastindex(Tcurr)
    
    Tp_curr = @view Tcurr[plate_begin:plate_end]
    Tp_prev = @view Tprev[plate_begin:plate_end]
    Tf_curr = @view Tcurr[fluid_begin:fluid_end]
    Tf_prev = @view Tprev[fluid_begin:fluid_end]
    Tl_curr = @view Tcurr[tank_begin:tank_end]
    Tl_prev = @view Tprev[tank_begin:tank_end]
    @assert length(Tp_curr) == length(Tf_curr) == Nc
    @assert length(Tl_curr) == Ns

    c1Δt = c1*Δt
    c2Δt = c2*Δt
    c4Δt = c4*Δt
    
    # update plate temperature
    Tp_curr .= Tp_prev .+ c1Δt .* (S .- hpf .* (Tp_prev .- Tf_prev) .- hpa .* (Tp_prev .- Ta) .- α .* (Tp_prev .^ 4 .- (T∞^4)))
    # update collector fluid temperature
    Tf_curr[begin] = Tf_prev[begin] + c2Δt *(c3*(Tp_prev[begin] - Tf_prev[begin]) - mdot* (Tf_prev[begin+1] - Tl_prev[end])/(2Δy))
    Tf_curr_m = @view Tf_curr[begin+1:end-1]  # exclude terminal entries
    Tf_prev_m = @view Tf_prev[begin+1:end-1]
    Tp_curr_m = @view Tp_curr[begin+1:end-1]
    Tp_prev_m = @view Tp_prev[begin+1:end-1]
    Tf_prev_m_sl = @view Tf_prev[begin:end-2]  # `_sl` ≡ shifted left
    Tf_prev_m_sr = @view Tf_prev[begin+2:end]  # `_sr` ≡ shifted right
    Tf_curr_m .= Tf_prev_m .+ c2Δt .* (c3.*(Tp_prev_m .- Tf_prev_m) .- mdot .* (Tf_prev_m_sr .- Tf_prev_m_sl)/(2Δy))
    Tf_curr[end] = Tf_prev[end] + c2Δt*(c3*(Tp_prev[end] - Tf_prev[end]) - mdot * (Tl_prev[end] - Tf_prev[end-1])/(2Δy))
    # update tank temperature
    Tl_curr[begin] = Tl_prev[begin] + c4Δt * (mdot * (Tf_prev[end] - Tl_prev[begin]) - c5*(Tl_prev[begin] - Ta))
    Tl_curr_m = @view Tl_curr[begin+1:end]
    Tl_prev_m = @view Tl_prev[begin+1:end]
    Tl_prev_m_sl = @view Tl_prev[begin:end-1]
    Tl_curr_m .= Tl_prev_m .+ c4Δt .* (mdot .* (Tl_prev_m_sl .- Tl_prev_m) - c5 .* (Tl_prev_m .- Ta))
    timeseries_data
end

function solve(prob::SWHSProblem; maxiter = prob.Nc^2)
    @unpack model, Nc, Tp0, Tf0 = prob
    @unpack Ns = model
    A = 4.0  # TODO: put this in the model
    N = 2Nc + Ns
    timeseries_data = zeros(N,maxiter)
    timeseries_data[begin:begin+Nc-1,:] .= Tp0
    timeseries_data[Nc+1:end] .= Tf0
    for i in 2:maxiter
        step!(timeseries_data, prob, i)
    end
    timeseries_data
end

function main()
    model = SWHSModel()
    prob = SWHSProblem(model)
    timeseries_data = solve(prob)
end

end # module
