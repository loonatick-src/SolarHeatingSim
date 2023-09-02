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

## Collector
Heat transfer between plate and fluid modelled as convection using Newton's law of cooling [1].
### Plate
We have the following processes
- Solar radiation, $\propto T_p^4 - T_{\text{sky}}^4$ (Stefan-Boltzmann law)
- Convective transfer to fluid, $\propto  T_p - T_f$ (Newton's law of cooling)
- Convective loss to atmosphere, $\propto T_p - T_a$ (Newton's law of cooling)

$$\rho_p \delta \frac{\partial T_p}{\partial t} = S - h_{pf}(T_p - T_f) - h_{pa}(T_p - T_a) - \alpha(T_p^4 - T_{\text{sky}}^4)$$

where $\alpha$ is a constant that includes the Stefan-Boltzmann constant, emissivities etc. This equation equation models heat transfer rate per unit area of the collector. Since we assume no dynamics along the width of the collector, we multiply the equation by the width $W$ of the collector to get the dynamics along the length of the conductor ($T_p = T_p(t, y)$).

$$W\rho_p \delta \frac{\partial T_p}{\partial t} = WS - Wh_{pf}(T_p - T_f) - Wh_{pa}(T_p - T_a) - W\alpha(T_p^4 - T_{\text{sky}}^4)$$

### Working Fluid
We model heat transfer through and by the fluid using a 1D advection equation with the convective transfer
## Storage Tank
[2]
# References
- [1] Incropera, Frank P., et al. Fundamentals of heat and mass transfer. Vol. 6. New York: Wiley, 1996.
- [2] Zeghib, I., & Chaker, A. (2011). Simulation of a solar domestic water heating system. Energy Procedia, 6, 292-301.
- [3] Kalita, K. (2020). Solar Energy Engineering and Technology [MOOC]. NPTEL. https://onlinecourses.nptel.ac.in/noc20_ph14/preview
- [4] Al-Tabbakh, A. A. (2022). Numerical transient modeling of a flat plate solar collector. Results in Engineering, 15, 100580