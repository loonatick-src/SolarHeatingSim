using UnPack
using ModelingToolkit, MethodOfLines, DomainSets

@kwdef struct AdvectionModel{T<:Number}
    ρp::T  = 8.0e3    # density of plate
    δ::T   = 0.1      # plate thickness
    W::T   = 1.0      # plate width
    L::T   = 2.0      # plate length
    H::T   = 5.0      # tank height
    Cpp::T = 450.0    # specific heat capacity of plate
    S::T   = 800.0    # radiation flux
    hpf::T = 1000.0   # heat transfer coefficient (plate and working fluid)
    hpa::T = 100.0    # heat transfer coefficient (plate and atmosphere)
    Ta::T  = 300.0    # ambient temperature
    T∞::T  = 295.0    # sky temperature
    α::T   = 5.5e-8   # radiation coefficient
    ρf::T  = 1.0e3    # working fluid density
    Af::T  = W*L*0.2      # transverse cross section area of risers
    Cpf::T = 4.2e3    # specific heat capacity of working fluid
    mdot::T = 0.5 * Af * ρf  # mass flow rate (assuming 0.5m per s in risers
    ρe::T  = 1.5e3    # effective density (tank + wokring fluid)
    Cpe::T = 5000.0   # effective specific heat of tank + working fluid
    # V::T   = 10.0     # volume of each stratified tank layer
    R::T   = 2.0      # inner radius of tank
    At::T  = π*R^2    # contact area between tank and atmosphere per stratification
    C::T = 2π*R       # circumference of tank
    hta::T = 500.0    # loss coefficient (storage tank)
    # Ns::Int= 10       # stratifications count for tank
    # Nc::Int= 100     # number of points in the discretized mesh of collector
end

value_type(::AdvectionModel{T}) where T = T

σ(x, a, b) = 1/(1 + exp(-(x-a)*b))

@kwdef struct AdvectionProblem{T}
    model::AdvectionModel{T}
    u::Vector{T}
    du::Vector{T}
    Tp0::T
    Tf0::T
    Nc::Int
    Nt::Int
    Δyc::T
    Δyt::T
    function AdvectionProblem(model::AdvectionModel, Nc::Int, Nt::Int, Tp0, Tf0)
        @unpack L, H = model
        T = value_type(model)
        N = 2Nc + Nt
        Δyc = L/Nc
        Δyt = H/Nt
        u = Vector{T}(undef, N)
        u[1:Nc] .= Tp0
        u[Nc+1:end] .= Tf0
        du = zeros(T, N)
        new{T}(u, du, Tp0, Tf0, Nc, Nt, Δyc, Δyt)
    end

    AdvectionProblem(model::AdvectionModel; Nc=10, Nt=10, Tp0=280.0, Tf0=280.0) = AdvectionProblem(model, Nc, Nt, Tp0, Tf0)
    AdvectionProblem() = AdvectionProblem(AdvectionModel())
end

function fluid_f(y, Tp, Tf;prob::AdvectionProblem = AdvectionProblem())
    @unpack L, H, W, hpf, C, hta, Ta = model
    y < L ? W*hpf*(Tp - Tf) : -C*hta*(Tf - Ta)
end

plate_f(Tp, Tf, W, hpf) = W*hpf*(Tp - Tf)
tank_f(Tf, C, hta, Ta) = -C*hta*(Tf - Ta)

fluid_f(y, Tp, Tf, L, H, W, hpf, C, hta, Ta) = (1-σ(y,L,1))*plate_f(Tp, Tf, W, hpf) + σ(y,L,1)*tank_f(Tf, C, hta, Ta) # y < L ? W*hpf*(Tp - Tf) : -C*hta*(Tf - Ta)

# function advection_f(y, Tp, Tf; prob::AdvectionProblem=AdvectionProblem())
#     @unpack model = prob
#     @unpack L, H, W, hpf, C, hta, Ta = model
#     @unpack ρf, Cpf, At, Af, mdot = model
#     1/(ρf*Cpf*A(y)) * (fluid_f(y, Tp, Tf) - mdot *Cpf)
# end

value_type(::AdvectionProblem{T}) where T = T

# """Area of cross section"""
# function A(y::T, model::AdvectionModel{T}) where T
#     @unpack L, H, Af, A = model
#     (y < 0 || y > L + H) && throw(DomainError(y, "Expected $(0) < y < $(L+H)."))
#     y < L ? Af : At
# end

A(y, A1, A2, L, H)  = (1-σ(y,L,1))*A1 + σ(y,L,1)*A2 # y < L ? A1 : A2


function pde_equations(tmin, tmax, Tf0, Tp0, dy; Lval=2.0, Hval=5.0)
    @parameters y t
    @parameters L H S δ Cpp hpf kp hpa Ta T∞ ρf Cpf mdot Af At ρp W C hta
    @variables Tp(..) Tf(..)
    Tp_init(y,t) = Tp0
    Tf_init(y,t) = Tf0
    Dt = Differential(t)
    Dy = Differential(y)
    domains = [t ∈ Interval(tmin, tmax),
               y ∈ Interval(0, Lval+Hval)]
    eq = [
        Dt(Tp(y,t)) ~ 1/(ρp*δ*Cpp) * (S - hpf*(Tp(y,t) - Tf(y,t)) - hpa*(Tp(y,t) - Ta) #=- α*(Tp(y,t)^4 - T∞^4)=#),
        Dt(Tf(y,t)) ~ 1/(ρf*Cpf*A(y,L,H,Af,At))*(fluid_f(y,Tp(y,t), Tf(y,t), L, H, W, hpf, C, hta, Ta) - mdot * Cpf * Dy(Tf(y,t)))
    ]
    bcs = [
        Tp(y,0) ~ Tp_init(y,0),  # initial values
        Tf(y,0) ~ Tf_init(y,0),  # initial values
        Tf(0,t) ~ Tf(Lval+Hval,t)      # periodic boundary conditions
    ]
    @named pdesys = PDESystem(eq, bcs, domains,
                              [y,t], [Tp(y,t), Tf(y,t)],
                              [ρp => 8.0e3, δ => 0.1, W => 1.0,
                               L => Lval, H => Hval, Cpp => 450.0,
                               S => 800.0, hpf => 1000.0,
                               hpa => 100.0, Ta => 300, T∞ => 295,
                               ρf => 1.0e3, Af => 1, Cpf=>4.2e3,
                               mdot=>500, At=>12.5, C=>12.5, hta=>500.0])
    discretization = MOLFiniteDifference([y=>dy], t)
    order=2
    prob = discretize(pdesys, discretization)
    # eq, domains, bcs
    prob
end
