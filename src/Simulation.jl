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
    c4::T
    c5::T
    Tf0::T
    Tp0::T
    function SWHSProblem(model::SWHSModel; Nc = 100, Δt = 0.1, Tf0 = 290.0, Tp0 = 300.0)
        T = value_type(model)
        @unpack Cpp, ρp, δ, L, Ns, mdot, V, Cpe = model
        @unpack W, hpf, ρf, Af, Cpl, ρt, hta = model
        Δy = L / Nc
        c1 = 1/(Cpp * ρp * δ)
        c2 = 1/(ρf * Af)
        c3 = W*hpf/Cpf
        c4 = 1/(ρe*V)
        c5 = hta * A / Cpe
        new{Float64}(model, Δy, Δt, Nc, c1, c2, c3, c4, Tf0, Tp0)
    end
end

function step!(timeseries_data, prob::SWHSProblem, i)
    @unpack Δt, Δy, Nc, c1, c2, c3, c4 = prob
    m = prob.model
    @unpack α = m
    Tprev = @view timeseries_data[:,i-1]
    Tcurr = @view timeseries_data[:,i]
    Tp_prev = @view Tprev[begin:begin+Nc-1]
    Tp_curr = @view Tcurr[begin:begin+Nc-1]
    @assert length(Tp_prev) == Nc
    Tf_prev = @view Tprev[Nc+1:2Nc]
    Tf_curr = @view Tcurr[Nc+1:2Nc]
    @assert length(Tf_prev) == Nc
    Tl_prev = @view Tprev[2Nc+1:end]
    Tl_curr = @view Tcurr[2Nc+1:end]
    @assert length(Tl_prev) == m.Ns
    Tl_prev_sl = @view Tprev[2Nc:end-1]     # `_sl` ≡ shifted left
    @assert length(Tl_prev_sl) == length(Tl_prev)
    Tf_prev_sl = @view Tprev[Nc:2Nc-1]
    Tf_prev_sr = @view Tprev[Nc+2:2Nc+1]
    @assert length(Tf_prev_sl) == length(Tf_prev_sr) == length(Tf_prev)
    # you should see a one-to-one correspondence with the fully discretized equations in the README.md
    @unpack S, hpf, hpa, Ta, T∞ = m
    Tp_curr .= Tp_prev .+ (Δt / c1) .* (S .- hpf .* (Tp_prev .- Tf_prev) .- hpa .* (Tp_prev .- Ta) .- α .* (Tp_prev.^4 .- T∞.^4))
    Tf_curr .= Tf_prev .+ Δt .* (c2 .* (Tp_prev .- Tf_prev) - (Tf_prev_sr .- Tf_prev_sl) ./ 2Δy)
    Tl_curr .= Tl_prev .+ Δt .* (c3 .* (Tl_prev_sl .- Tl_prev) .- c4 .* (Tl_prev .- Ta))
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
