module Simulation

using UnPack

export SWHSModel, SWHSProblem, solve

struct SWHSModel{T}
    S::T    # radiation flux
    Ta::T   # ambient temperature
    T∞::T   # sky temperature
    α::T    # radiation coefficient
    δ::T    # plate thickness
    W::T    # plate width
    L::T    # plate length
    Cpl::T  # specific heat of water
    Cpp::T  # specific heat of plate
    Cpt::T  # specific heat of tank
    ρl::T   # density of water
    ρp::T   # plate denisty
    ρt::T   # density of tank walls
    Vt::T   # volume of each stratified tank layer
    Ac::T   # transverse cross section area of risers
    ρAv::T  # mass flow rate of water
    hta::T   # loss coefficient (storage tank)
    hpf::T  # plate-fluid heat transfer coefficient
    hpa::T  # plate-atmosphere transfer coefficient
    Ns::Int # stratifications count

    # ctor with default arguments
    function SWHSModel(;S      = 800.0,
                       Ta      = 300,
                       T∞      = 295,
                       α       = 5.5e-8,
                       δ       = 0.1,
                       W       = 1.0,
                       L       = 2.0,
                       Cpl     = 4.2e3,
                       Cpp     = 450.0,
                       Cpt     = 5000.0,
                       ρl      = 1.0e3,
                       ρp      = 8.0e3,
                       ρt      = 1.5e3,
                       Vt      = 10.0,
                       Ac      = 0.6,
                       ρAv     = 800.0,
                       hta     = 500.0,
                       hpf     = 1000.0,
                       hpa     = 100.0,
                       Ns::Int = 10)
        new{Float64}(S,
                     Ta,
                     T∞,
                     α,
                     δ,
                     W,
                     L,
                     Cpl,
                     Cpp,
                     Cpt,
                     ρl,
                     ρp,
                     ρt,
                     Vt,
                     Ac,
                     ρAv,
                     hta,
                     hpf,
                     hpa,
                     Ns)
    end
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
    Tf0::T
    Tp0::T
    function SWHSProblem(model::SWHSModel; Nc = 100, Δt = 0.1, Tf0 = 290.0, Tp0 = 300.0)
        T = value_type(model)
        @unpack Cpp, ρp, δ, L, Ns, ρAv, Vt, Cpt = model
        @unpack W, hpf, ρl, Ac, Cpl, ρt, hta = model
        A = 4.0  # TODO: put this in the model
        Δy = L / Nc
        c1 = Cpp * ρp * δ
        c2 = W*hpf/(ρl * Ac * Cpl)
        c3 = ρAv / (ρt*Vt)
        c4 = hta * A/(ρt * Vt * Cpt)
        N = 2Nc + model.Ns
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
