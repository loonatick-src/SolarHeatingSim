using UnPack

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
    Ut::T   # loss coefficient (storage tank)
    Up::T   # loss coefficient (collector plate)
    hpf::T  # plate-fluid heat transfer coefficient
    hpa::T  # plate-atmosphere transfer coefficient
    Ns::Int # stratifications count
    function SWHSModel(;S = 800.0, Ta = 300, T∞ = 295, α = 5.5e-8,
                       δ = 0.1, W = 1.0, L = 2.0, Cpl = 4.2e3,
                       Cpp = 450.0, Cpt = 5000.0, ρl = 1.0e3, ρp = 8.0e3,
                       ρt = 1.5e3, Vt = 10.0, Ac = 0.6, ρAv = 800.0, Ut = 500.0,
                       hpf = 1000.0, hpa = 100.0, Ns = 10)
        new{Float64}(S, Ta, T∞, α, δ, W,
                     L, CPl, Cpp, Cpt,
                     ρl, ρp, ρt, Vt, ρAv,
                     Ut, Up, hpf, hpa, Ns)
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
    Tf0::T
    Tp0::T
    buf1::Vector{T}
    buf2::Vector{T}

    function SWHSProblem(model::SWHSModel; Nc = 100, Δt = 0.1, Tf0 = 290.0, Tp0 = 300.0)
        T = value_type(model)
        @unpack Cpp, ρp, δ, L, Nc, Ns = model
        @unpack W, hpf, ρf, Ac, Cpl = model
        Δy = L / Nc
        c1 = Cpp * ρp * δ
        c2 = W*hpf/(ρf * Ac * Cpl)
        N = 2Nc + model.Ns
        buf1 = Vector{T}(undef, N)
        buf2 = similar(buf1)
        buf1[begin:begin+Nc-1] .= Tp0
        buf1[Nc:end] = Tf0
        new{Float64}(Δy, Δt, Nc, c1, c2, Tf0, Tp0, buf1, buf2)
    end
end

function step!(timeseries_data, prob::SWHSProblem, i)
    # TODO
end

function solve(prob::SWHSProblem; maxiter = prob.Nc^2)
    # TODO
end
