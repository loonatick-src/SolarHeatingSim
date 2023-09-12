module Simulation

using UnPack
using LinearAlgebra: dot, mul!

export SWHSModel, SWHSProblem
export forward_euler, RK4, f!
export get_plate_view, get_fluid_view, get_tank_view
export pde_equations

include("utilities.jl")

@kwdef struct SWHSModel{T}
    Ns::Int= 10       # stratifications count for tank
    Nc::Int= 100      # number of points in the discretized mesh of collector    
    ρp::T  = 8.0e3    # density of plate
    ρf::T  = 1.0e3    # working fluid density
    ρe::T  = 1.5e3    # effective density (tank + wokring fluid)
    δ::T   = 0.1      # plate thickness
    W::T   = 1.0      # plate width
    L::T   = 2.0      # plate length
    H::T   = 10.0     # height of tank
    ΔH::T  = H/Ns     # height of each stratified layer
    Δy::T  = L/Nc
    dt::T  = 5.0      # tank diameter         
    dr::T  = 0.1      # riser pipe diameter
    dd::T  = dr       # down comer pipe diameter    
    kp::T  = 50.0     # thermal conductivity of plate
    Cpp::T = 450.0    # specific heat capacity of plate
    Cpf::T = 4.2e3    # specific heat capacity of working fluid
    Cpe::T = 5000.0   # effective specific heat of tank + working fluid
    S::T   = 800.0    # radiation flux
    hpf::T = 1000.0   # heat transfer coefficient (plate and working fluid)
    hpa::T = 100.0    # heat transfer coefficient (plate and atmosphere)
    hra::T = 100.0    # heat transfer coefficient (riser pipe and ambient)
    hda::T = hra      # heat transfer coefficient (down pipe and ambient)
    hta::T = 500.0    # loss coefficient (storage tank)    
    Ta::T  = 300.0    # ambient temperature
    T∞::T  = 295.0    # sky temperature
    α::T   = 5.5e-8   # radiation coefficient
    Af::T  = W*L*0.2  # transverse cross section area pipes in collector
    mdot::T = 0.5 * Af * ρf  # mass flow rate (assuming 0.5m per s in collector pipes
    A::T   = π*dt*ΔH  # ambient contact area per stratified layer
    V::T   = ΔH*π*(dt^2)/4  # tank volume per stratified layer
    ρCl::T = 1.0e3    # TODO: decide value
    mCr::T = 1.0e3    # TODO: decide value
    mCd::T = mCr
    plate_begin::Int = 1
    plate_end::Int = plate_begin + Nc - 1
    dc_idx::Int = plate_end+1
    fluid_begin::Int = dc_idx+1
    fluid_end::Int = fluid_begin + Nc - 1
    r_idx::Int = fluid_end+1
    tank_begin::Int = r_idx+1
    tank_end::Int = tank_begin+Ns-1
end

value_type(::SWHSModel{T}) where T = T

@kwdef struct SWHSCache
    model::SWHSModel    = SWHSModel()
    # finite difference (FD) matrices and the FD-like matrix
    Dp2::AbstractMatrix = plate_FD2_matrix(model.Δy, model.Nc)
    Df1::AbstractMatrix = fluid_FD1_matrix(model.Δy, model.Nc)
    Dl::AbstractMatrix  = tank_matrix(value_type(model), model.Ns)
    # pre-allocated buffers
    Dp2_Tp::AbstractVector = zeros(value_type(model), model.Nc)
    Df1_Tf::AbstractVector = zeros(value_type(model), model.Nc)
    Dl_Tl::AbstractVector = zeros(value_type(model), model.Ns)
    # coefficients used in the equations
    ρp_δ_Cpp_inv::Number    = 1/(model.ρp * model.δ * model.Cpp)
    δ_kp::Number        = model.δ * model.kp
    ρf_Af_Cpf_inv::Number   = 1/(model.ρf * model.Af * model.Cpf)
    mdot_Cpf::Number    = model.mdot * model.Cpf
    π_dr_hra::Number    = π*model.dr * model.hra
    π_dd_hda::Number    = π*model.dd * model.hda
end

struct SWHSProblem
    u::AbstractVector
    du::AbstractVector
    cache::SWHSCache
    function SWHSProblem(model::SWHSModel; Tf0 = 290.0, Tp0 = 300.0, Td0=Tf0, Tr0=Tf0)
        T = value_type(model)
        @unpack Nc, Ns, plate_begin, plate_end, fluid_begin, fluid_end, tank_begin, tank_end, dc_idx, r_idx = model
        N = 2Nc + Ns + 2
        u = zeros(T, N)
        du = zeros(T, N)
        @assert length(@view u[plate_begin:plate_end]) == Nc
        u[plate_begin:plate_end] .= Tp0
        u[dc_idx] = Td0
        @assert length(@view u[fluid_begin:fluid_end]) == Nc
        u[fluid_begin:fluid_end] .= Tf0
        u[r_idx] = Tr0
        @assert length(@view u[tank_begin:tank_end]) == Ns
        u[tank_begin:tank_end] .= Tf0
        cache = SWHSCache(;model=model)
        new(u, du, cache)
    end
end

function f!(du, u, cache::SWHSCache)
    @unpack model = cache
    @unpack plate_begin, plate_end, fluid_begin, fluid_end, tank_begin, tank_end, dc_idx, r_idx = model
    @unpack S, hpf, hpa, α, T∞, Ta, W, ρCl, mCr, mCd, mdot = model
    @unpack ρp_δ_Cpp_inv, ρf_Af_Cpf_inv, δ_kp, Dp2, Dp2_Tp, Df1, Df1_Tf, Dl, Dl_Tl = cache
    @unpack π_dd_hda, π_dr_hra, mdot_Cpf = cache
    u_p  = @view u[plate_begin:plate_end]
    u_f  = @view u[fluid_begin:fluid_end]
    u_l  = @view u[tank_begin:tank_end]
    du_p = @view du[plate_begin:plate_end]
    du_f = @view du[fluid_begin:fluid_end]
    du_l = @view du[tank_begin:tank_end]
    Td = u[dc_idx]
    Tr = u[r_idx]
    mul!(Dp2_Tp, Dp2, u_p)
    mul!(Df1_Tf, Df1, u_f)
    Df1_Tf[begin] -= Td
    Df1_Tf[end] += Tr
    mul!(Dl_Tl, Dl, u_l)
    Dl_Tl[begin] += Tr
    Tp⁴ = dot(u_p,u_p)^2
    Tr = u[r_idx]
    Td = u[dc_idx]
    W_hpf = W * hpf
    @. du_p = ρp_δ_Cpp_inv * (S + δ_kp*Dp2_Tp - hpf*(u_p - u_f) - hpa*(u_p - Ta) - α*(Tp⁴ - T∞))
    @. du_f = ρf_Af_Cpf_inv * (W_hpf * (u_p - u_f) - mdot_Cpf*Df1_Tf)
    @. du_l = 1/ρCl * Dl_Tl - π_dr_hra * (Tr - Ta)
    du[r_idx] = (mdot / mCr) * (u_f[end] - Tr) - (π_dr_hra/mCr) * (Tr - Ta)
    du[dc_idx] = (mdot / mCd) * (u_l[end] - Td) - (π_dd_hda/mCd) * (Td - Ta)
end

function reset_problem!(prob::SWHSProblem; Tf0 = 290.0, Tp0 = 300.0)
    @unpack Nc, Ns = prob.model
    @unpack u, du = prob
    u[1:Nc] .= Tp0
    u[Nc+1:end] .=   Tf0
    prob
end



function forward_euler(f!, tspan, dt, prob::SWHSProblem)
    @unpack u, du = prob
    t = first(tspan)
    while t < last(tspan)
        f!(du, u, prob.cache)
        u .+= du
        t += dt
    end
    u
end

function RK4(f!, tspan, dt, prob::SWHSProblem)
    @unpack u, du = prob
    @unpack cache = prob
    T = eltype(u)
    k1 = zeros(T, length(u))
    k2 = zeros(T, length(u))
    k3 = zeros(T, length(u))
    k4 = zeros(T, length(u))
    t = first(tspan)
    while t < last(tspan)
        f!(k1, u, cache)
        f!(k2, u .+ dt .* k1 ./ 2, cache)
        f!(k3, u .+ dt .* k2 ./ 2, cache)
        f!(k4, u .+ dt .* k3, cache)
        @. u += dt/6 * (k1 + 2k2 + 2k3 + k4)
        t += dt
    end
    u
end

function plate_view(prob::SWHSProblem)
    @unpack plate_begin, plate_end = prob.cache.model
    @view u[plate_begin:plate_end]
end


function fluid_view(prob::SWHSProblem)
    @unpack fluid_begin, fluid_end = prob.cache.model
    @view u[fluid_begin:fluid_end]
end

function tank_view(prob::SWHSProblem)
    @unpack tank_begin, tank_end = prob.cache.model
    @view u[tank_begin:tank_end]
end

include("example.jl")

end # module
