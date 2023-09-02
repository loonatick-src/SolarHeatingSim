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
    ρAv::T  # mass flow rate of water
    Ut::T   # loss coefficient (storage tank)
    Up::T   # loss coefficient (collector plate)
    hpf::T  # plate-fluid heat transfer coefficient
    hpa::T  # plate-atmosphere transfer coefficient
    Ns::Int # stratifications count
    function SWHSModel(;S = 800.0, Ta = 300, T∞ = 295, α = 5.5e-8,
                       δ = 0.1, W = 1.0, L = 2.0, Cpl = 4.2e3,
                       Cpp = 450.0, Cpt = 5000.0, ρl = 1.0e3, ρp = 8.0e3,
                       ρt = 1.5e3, Vt = 10.0, ρAv = 800.0, Ut = 500.0,
                       hpf = 1000.0, hpa = 100.0, Ns = 10)
        new{Float64}(S, Ta, T∞, α, δ, W,
                     L, CPl, Cpp, Cpt,
                     ρl, ρp, ρt, Vt, ρAv,
                     Ut, Up, hpf, hpa, Ns)
    end
end

function semidiscretize(model::SWHSModel; Nc = 100)
    Δy = model.L / Nc
end
