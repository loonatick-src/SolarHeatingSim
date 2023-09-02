# Rudimentary Solar Water Heating System Simulation
We make the following assumptions and simplifications
- Flat plate collector with parallel riser tubes
- Temperature variations along the width of the solar collector (i.e. perpendicular to fluid flow direction) are not modelled, assume uniform mean temperature along with.
- Single phase regime -- the fluid does not undergo phase change
- Changes in specific heat capacity with temperature are negligible for all materials in the operating range
- Fixed ambient temperature, wind speed (i.e. fixed convection heat transfer coefficients), sky temperature, and solar irradiation
- Working fluid is incompressible, mass flow rate is constant throughout the system
- Heat exchanges occur only at the collector and storage tank, no heat transfer at pipes or pump
- The work done by the pump is completely balanced by e.g. the changes in potential energy in the storage tank and riser tubes, so it does not appear in the balance equations
- Radial temperature of the fluid is not modelled - assume an overall mean temperature across the cross section

## Model Parameters
```julia
struct SWHSModel{T}
    S::T    # radiation flux
    Ta::T   # ambient temperature
    T∞::T   # sky temperature
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
    Acs::T  # fluid flow cross section area
end
```
## Subsystems
### Collector
We consider a flat plate collector with parallel riser tubes (TODO: add figure) modelled as an absorber plate transfering heat to the fluid in the risers. This model is based on the work of [4]. Heat transfer between plate and fluid modelled as convection using Newton's law of cooling [1].
#### Plate
We have the following processes
- Solar radiation, $\propto T_p^4 - T_{\text{sky}}^4$ (Stefan-Boltzmann law)
- Convective transfer to fluid, $\propto  T_p - T_f$ (Newton's law of cooling)
- Convective loss to atmosphere, $\propto T_p - T_a$ (Newton's law of cooling)

$$\rho_p \delta \frac{\partial T_p}{\partial t} = S - h_{pf}(T_p - T_f) - h_{pa}(T_p - T_a) - \alpha(T_p^4 - T_{\text{sky}}^4)$$

where $\alpha$ is a constant that includes the Stefan-Boltzmann constant, emissivities etc. This equation equation models heat transfer rate per unit area of the collector. Since we assume no dynamics along the width of the collector, we multiply the equation by the width $W$ of the collector to get the dynamics along the length of the conductor ($T_p = T_p(t, y)$).

$$W\rho_p \delta \frac{\partial T_p}{\partial t} = WS - Wh_{pf}(T_p - T_f) - Wh_{pa}(T_p - T_a) - W\alpha(T_p^4 - T_{\text{sky}}^4)$$

#### Working Fluid
We model heat transfer through and by the fluid using a 1D advection equation with the source term being the plate-fluid convective transfer term.

$$\rho_fA_{\text{t}}C_{pf}\frac{\partial T_f}{\partial t} + \rho_f A_{\text{t}} C_{pf} V\frac{\partial T_f}{\partial y} = Wh_pf(T_p - T_f),$$

where $A_t$ is the transverse flow cross section.

### Storage Tank
We use a stratified, well-mixed tank model (TODO: add figure)
$$\rho_e V C_{pe} \frac{dT_{l,i}}{dt} = \dot{m} C_{pe} (T_{l,i-1} - T_{l,i}) - h_{ta}A(T_{l,i} - T_a),$$
for $i \in 1 \ldots N_s$, where $N_s$ is the number of stratification layers, (model parameter). We let $T_{l,0}(t) = T_f(t, L)$. The assumption here is that the pipes are short or transport from collector to tank is sufficiently fast.
Note that $T_l$ is the temperature of the same working fluid as in the collector, but we use a different subscript to distinguish between the two systems.
This model is based on lecture 29 of [3].

## Numerical Considerations
We use central differences for spatial derivatives in the collector equations.

$$\frac{\partial T_f}{\partial y} \approx \frac{T_{f,i+1} - T_{f,i-1}}{2\Delta y} \quad \forall i \in 1,\ldots N_c,$$

where $T_{f,i} = T_f(t, i\Delta y)$, and $N_c$ is the number of degrees of freedom in the discretized mesh. The spatial extremities of the collector
connect to the storage tank, so we take

$$T_{N_c+1} = T_{l,1}$$

$$T_{N_0} = T_{l, N_s}.$$

The semidiscretized equations are then
### Plate
$$W\rho_p \delta \frac{d T_{p,i}}{d t} = WS - Wh_{pf}(T_{p,i} - T_{f,i}) - Wh_{pa}(T_{p,i} - T_a) - W\alpha(T_{p,i}^4 - T_{\text{sky}}^4),\quad \forall i \in 1,\ldots N_c$$
### Fluid in collector

$$\rho_fA_{\text{t}}C_{pf}\frac{d T_{f,i}}{d t} + \rho_f A_{\text{t}} C_{pf} V\frac{T_{f,i+1} - T_{f,i-1}}{2\Delta y} = Wh_{pf}(T_{p,i} - T_{f,i}), \quad \forall i \in 1,\ldots N_c$$
with $T_{f,N_c+1} = T_{l,1}$ and $T_{f,0} = T_{l,N_s}$.

### Storage Tank

$$\rho V C_{pe} \frac{dT_{l,i}}{dt} = \dot{m} C_{pe} T_{l,i-1} - T_{l,i} - h_{ta}A(T_{l,i} - T_a), \quad \forall i \in 1,\ldots N_s$$
with $T_{l,0} = T_{f,N_c}$.

### Full discretization
For simplicity we choose to use forward difference approximation for the time derivatives.
This results in the explicit forward in time, central in space (FTCS) scheme. If we choose some other
finite difference scheme for temporal derivatives that results in an implicit scheme we will require
a nonlinear solver because of the radiative term. Alternatively we could neglect the radiative term and use
a linear solve instead.

The FTCS scheme then results in the following discretized forms.
#### Plate
$$T_{p,i}^{n+1} = \frac{\Delta t}{\rho_p \delta}\left(S - h_{pf}(T_{p,i}^n - T_{f,i}^n) - h_{pa}(T_{p,i}^n - T_a) - \alpha \left((T_{p,i}^n)^4 - T_{\text{sky}}^4\right)\right) + T_{p,i}^n$$

$\forall i \in 1,\ldots N_c$, $n \in 0,\ldots$ till convergence.
#### Fluid in collector
$$T_{f,i}^{n+1} = \Delta t\left( \frac{Wh_{pf}}{\rho_f A_t C_{pf}}(T_{p,i}^n - T_{f,i}^n) - \frac{T_{f,i+1}^n - T_{f,i-1}^n}{2\Delta y}\right) + T_{f,i}^n$$

$\forall i \in 1,\ldots N_c$ and $j \in 1,\ldots$ till convergence.

#### Fluid in storage tank
$$T_{l,i}^{n+1} = \Delta t \left( \frac{\dot m}{\rho_e V}(T_{l,i-1}^n - T_{l,i}^n) - \frac{h_{ta}A}{\rho_e V C_{pe}}(T_{l,i}^n - T_a)  \right) + T_{l,i}^n$$



## Possible Refinements


##  References
- [1] Incropera, Frank P., et al. Fundamentals of heat and mass transfer. Vol. 6. New York: Wiley, 1996.
- [2] Zeghib, I., & Chaker, A. (2011). Simulation of a solar domestic water heating system. Energy Procedia, 6, 292-301.
- [3] Kalita, K. (2020). Solar Energy Engineering and Technology [MOOC]. NPTEL. https://onlinecourses.nptel.ac.in/noc20_ph14/preview
- [4] Al-Tabbakh, A. A. (2022). Numerical transient modeling of a flat plate solar collector. Results in Engineering, 15, 100580