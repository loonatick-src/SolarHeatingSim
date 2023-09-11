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
```
## Subsystems
### Collector
We consider a flat plate collector with parallel riser tubes (TODO: add figure) modelled as an absorber plate transfering heat to the fluid in the risers. This model is based on the work of [4]. Heat transfer between plate and fluid modelled as convection using Newton's law of cooling [1].
#### Plate
We have the following processes
- Solar radiation, $\propto T_p^4 - T_{\text{sky}}^4$ (Stefan-Boltzmann law)
- Convective transfer to fluid, $\propto  T_p - T_f$ (Newton's law of cooling)
- Convective loss to atmosphere, $\propto T_p - T_a$ (Newton's law of cooling)
$$\rho_p \delta C_{pp}\frac{\partial T_p}{\partial t} = S + \delta k_p \frac{\partial^2 T_p}{\partial y^2}- h_{pf}(T_p - T_f) - h_{pa}(T_p - T_a) - \alpha(T_p^4 - T_{\text{sky}}^4)$$

where $\alpha$ is a constant that includes the Stefan-Boltzmann constant, emissivities etc. This equation equation models heat transfer rate per unit area of the collector. Since we assume no dynamics along the width of the collector, we multiply the equation by the width $W$ of the collector to get the effective dynamics along the length of the conductor ($T_p = T_p(t, y)$). 

$$W\rho_p \delta C_{pp}\frac{\partial T_p}{\partial t} = WS - Wh_{pf}(T_p - T_f) - Wh_{pa}(T_p - T_a) - W\alpha(T_p^4 - T_{\text{sky}}^4)$$
Note that we only need the width factor to get the source term for the 1D advection equation used for simulating  the working fluid in the collector.

$$\frac{\partial T_p}{\partial t} = \frac{1}{\rho_p \delta C_{pp}}\left( S + \delta k_p \frac{\partial^2 T_p}{\partial y^2}- h_{pf}(T_p - T_f) - h_{pa}(T_p - T_a) - \alpha(T_p^4 - T_{\text{sky}}^4) \right)$$

We use Neumann boundary conditions:
$$\frac{\partial T_p}{\partial T_y} = 0\quad \text{at }y = 0, y = L$$
#### Working Fluid
We model heat transfer through and by the fluid using a 1D advection equation with the source term being the plate-fluid convective transfer term.

$$\rho_fA_f C_{pf}\frac{\partial T_f}{\partial t} + \rho_f A_f C_{pf} v_f\frac{\partial T_f}{\partial y} = Wh_{pf}(T_p - T_f),$$

where $A_f$ is the transverse flow cross section. Note that mass conservation requires $\rho_f A_f v_f = \dot m$ to be constant.

$$\rho_f A_c C_{pf}\frac{\partial T_f}{\partial t} + \dot{m} C_{pf}\frac{\partial T_f}{\partial y} = Wh_pf(T_p - T_f),$$
$$\implies \frac{\partial T_f}{\partial t} = \frac{1}{\rho_f A_c C_{pf}}\left(W h_{pf}(T_p - T_f) - \dot m C_{pf} \frac{\partial T_p}{\partial y}\right)$$

### Storage Tank
We use a stratified, well-mixed tank model (TODO: add figure)

$$(\rho C)_t V \frac{dT_{l,i}}{dt} = \dot{m} C_p(T_{l,i-1} - T_{l,i}) - h_{ta}A(T_{l,i} - T_a)$$
$$\implies \frac{dT_{l,i}}{dt} = \frac{1}{(\rho C)_t V}\left(\dot{m} C_p(T_{l,i-1} - T_{l,i}) - h_{ta}A(T_{l,i} - T_a)\right)$$

for $i \in 1 \ldots N_s$, where $N_s$ is the number of stratification layers, (model parameter). We let $T_{l,0}(t) = T_f(t, L)$. The assumption here is that the pipes are short or transport from collector to tank is sufficiently fast.
Note that $T_l$ is the temperature of the same working fluid as in the collector, but we use a different subscript to distinguish between the two systems.
This model is based on lecture 29 of [3].

### Connecting Pipes
We assume that the temperature drop across the connecting pipes is small and model them as single nodes, as done in (TODO: cite). Denote $\text{d}$ as the down comer pipe that leads flow into the collector, and $\text{r}$ as the up-riser pipe leading from the collector to the tank. $T_d$ = $T_f(0)$ and $T_r = T_{l,0}$.

$$(mC)_r\frac{dT_r}{dt} = \dot m(T_f(L) - T_r) - \pi d_r h_{ra}(T_r - T_a)$$

$$(mC)_d\frac{dT_d}{dt} = \dot m (T_{l,N_s} - T_d) - \pi d_dh_{da}(T_d - T_a)$$

## Numerical Considerations

$$\frac{\partial T}{\partial y} \approx \frac{T_{i+1} - T_{i-1}}{2\Delta y} \quad \forall i \in 1,\ldots N_c,$$

where $T_{i} = T(t, i\Delta y)$, and $N_c$ is the number of degrees of freedom in the discretized mesh.   We also need to take into account the boundary conditions.
### Plate
Consider first $T_p$.

```math
\begin{bmatrix}
T'_{p,1}\\
T_{p,2}'\\
T_{p,3}'\\
\vdots\\
T_{p,N_c-1}'\\
T_{p,N_c}'\end{bmatrix} \approx \frac{1}{2\Delta y}\begin{bmatrix}
0 & 0 & 0  &       & &\\
1 & 0 & -1 &       & &\\
  & 1 & 0 & -1     & &\\
  &   &   &        \ddots & &\\
  &   &   &        1      & 0 & -1\\
  &   &   &        0      & 0 & 0  
\end{bmatrix}
\begin{bmatrix}
T_{p,1}\\
T_{p,2}\\
T_{p,3}\\
\vdots\\
T_{p,N_c-1}\\
T_{p,N_c}
\end{bmatrix}
```
Where the top and bottom rows are zero because of the Neumann conditions on the plate boundaries. We need the second derivative for the conductive term as well.

```math
\begin{bmatrix}
T_{p,1}''\\
T_{p,2}''\\
T_{p,3}''\\
\vdots\\
T_{p,N_c-1}''\\
T_{p,N_c}''
\end{bmatrix} \approx \frac{1}{\Delta y ^2}
\begin{bmatrix}
-2 & 2  & 0\\
1  & -2 & 1\\
   & 1  & -2 &  1\\
   &    &    &  \ddots &\\
   &    &    &         & 1 & -2 & 1\\
   &    &    &         &   &  2 & -2
\end{bmatrix}
\begin{bmatrix}
T_{p,1}\\
T_{p,2}\\
T_{p,3}\\
\vdots\\
T_{p,N_c-1}\\
T_{p,N_c}
\end{bmatrix}
```
Where we estimated the second derivative at the boundaries using the Neumann conditions.
$$\frac{T_{p,N_c+1} - T_{p,N_c-1}}{2\Delta y} = 0 \implies T_{p,N_c+1} = T_{p,N_c-1}$$
$$\implies T''_{p,N_c} \approx \frac{T_{p,N_c+1} + T_{p,N_c-1} - 2T_{p,N_c}}{\Delta y ^2} = \frac{2T_{p,N_c-1} - 2T_{p,N_c}}{\Delta y^2}$$
and similarly for $T''_{p,1}$.
### Fluid in Collector
We have the boundary conditions $T_{f,1} = T_d \implies T'_{f,1} = T'_d$, and we can take $T_{f,N_c+1} = T_r$. Therefore, we have
```math
\begin{bmatrix}
T'_{f,1}\\
T'_{f,2}\\
T'_{f,3}\\
\vdots\\
T'_{f,N_c-1}\\
T'_{f,N_c}
\end{bmatrix} \approx \frac{1}{2\Delta y}
\begin{bmatrix}
0 & 0\\
1 & 0 & -1\\
  & 1 & 0 & -1\\
  &   &   &  \ddots\\
  &   &   &        & 1 & 0 & -1\\
  &   &   &        &   & -1 & 0\\
\end{bmatrix}
\begin{bmatrix}
T_{f,1}\\
T_{f,2}\\
T_{f,3}\\
\vdots\\
T_{f,N_c-1}\\
T_{f,N_c}
\end{bmatrix}
+
\begin{bmatrix}
T'_d\\
0\\
0\\
\vdots\\
0\\
\frac{T_r}{2\Delta y}
\end{bmatrix}
```

### Putting it all Together
**TODO**

##  References
- [1] Incropera, Frank P., et al. Fundamentals of heat and mass transfer. Vol. 6. New York: Wiley, 1996.
- [2] Zeghib, I., & Chaker, A. (2011). Simulation of a solar domestic water heating system. Energy Procedia, 6, 292-301.
- [3] Kalita, K. (2020). Solar Energy Engineering and Technology [MOOC]. NPTEL. https://onlinecourses.nptel.ac.in/noc20_ph14/preview
- [4] Al-Tabbakh, A. A. (2022). Numerical transient modeling of a flat plate solar collector. Results in Engineering, 15, 100580
- [5] Rosales, R. R. (2022) Notes: von Neumann Stability Analysis. [Course Notes] https://math.mit.edu/classes/18.300/Notes/Notes_vNSA.pdf
- [6] Rowell, D., & Wormley, D. N. (1997). System dynamics: an introduction (Vol. 635). Upper Saddle River: Prentice Hall.

# Appendix
(Everything beyond this point is one of the earlier attempts that missed a few details. Included for... proof of work? Not sure might delete later)

## Numerical Considerations
We use central differences for spatial derivatives in the collector equations.

$$\frac{\partial T_f}{\partial y} \approx \frac{T_{f,i+1} - T_{f,i-1}}{2\Delta y} \quad \forall i \in 1,\ldots N_c,$$

where $T_{f,i} = T_f(t, i\Delta y)$, and $N_c$ is the number of degrees of freedom in the discretized mesh. The spatial extremities of the collector
connect to the storage tank, so we take

$$T_{N_c+1} = T_{l,1}$$

$$T_{N_0} = T_{l, N_s}.$$

The semidiscretized equations are then
### Plate
$$\rho_p \delta C_{pp}\frac{d T_{p,i}}{d t} = S - h_{pf}(T_{p,i} - T_{f,i}) - h_{pa}(T_{p,i} - T_a) - \alpha(T_{p,i}^4 - T_{\text{sky}}^4),\quad \forall i \in 1,\ldots N_c$$
### Fluid in collector

$$\rho_fA_{f}C_{pf}\frac{d T_{f,i}}{d t} + \dot{m}C_{pf} \frac{T_{f,i+1} - T_{f,i-1}}{2\Delta y} = Wh_{pf}(T_{p,i} - T_{f,i}), \quad \forall i \in 1,\ldots N_c$$
with $T_{f,N_c+1} = T_{l,1}$ and $T_{f,0} = T_{l,N_s}$.

### Storage Tank
$$\rho_e V C_{pe} \frac{dT_{l,i}}{dt} = \dot{m} C_{pe} (T_{l,i-1} - T_{l,i}) - h_{ta}A(T_{l,i} - T_a), \quad \forall i \in 1,\ldots N_s$$
with $T_{l,0} = T_{f,N_c}$.
### The "Edge Cases"
To be explicit, we write enumerate the collector fluid and storage tank equations explicitly.
#### Fluid in Collector
$i = 1$
$$\rho_f A_f C_{pf} \frac{dT_{f,1}}{dt} + \dot m C_{pf}\frac{T_{f,2} - T_{l,N_s}}{2\Delta y} = Wh_{pf}(T_{p,1} - T_{f,1})$$
$i \in 2, \ldots N_c-1$
$$\rho_f A_f C_{pf} \frac{dT_{f,i}}{dt} + \dot m C_{pf}\frac{T_{f,i+1} - T_{f,i-1}}{2\Delta y} = Wh_{pf}(T_{p,i} - T_{f,i})$$
$i = N_c$
$$\rho_f A_f C_{pf} \frac{dT_{f,N_c}}{dt} + \dot m C_{pf}\frac{T_{l,0} - T_{f,N_c-1}}{2\Delta y} = Wh_{pf}(T_{p,N_c} - T_{f,N_c})$$
#### Storage Tank
$i = 1$
$$\rho_e V C_{pe} \frac{dT_{l,1}}{dt} = \dot{m} C_{pe} (T_{f,N_c} - T_{l,1}) - h_{ta}A(T_{l,1} - T_a)$$ 
$i \in 2,\ldots N_s$
$$\rho_e V C_{pe} \frac{dT_{l,i}}{dt} = \dot{m} C_{pe} (T_{l,i-1} - T_{l,i}) - h_{ta}A(T_{l,i} - T_a)$$

### Full discretization
For simplicity we choose to use forward difference approximation for the time derivatives.
This results in the explicit forward in time, central in space (FTCS) scheme. If we choose some other
finite difference scheme for temporal derivatives that results in an implicit scheme we will require
a nonlinear solver because of the radiative term. Alternatively we could neglect the radiative term and use
a linear solve instead.

The FTCS scheme then results in the following discretized forms.
#### Plate
$$T_{p,i}^{n+1} = T_{p,i}^n + \frac{\Delta t}{\rho_p C_{pp}\delta}\left(S - h_{pf}(T_{p,i}^n - T_{f,i}^n) - h_{pa}(T_{p,i}^n - T_a) - \alpha \left((T_{p,i}^n)^4 - T_{\text{sky}}^4\right)\right)$$

$\forall i \in 1,\ldots N_c$, $n \in 0,\ldots$ till convergence.
#### Fluid in collector
$$T_{f,i}^{n+1} = T_{f,i}^n + \frac{\Delta t}{\rho_f A_f} \left( \frac{Wh_{pf}}{C_{pf}}(T_{p,i}^n - T_{f,i}^n) - \dot m\frac{T_{f,i+1}^n - T_{f,i-1}^n}{2\Delta y}\right)$$

$\forall i \in 1,\ldots N_c$ and $j \in 1,\ldots$ till convergence.

#### Fluid in storage tank
$$T_{l,i}^{n+1} = T_{l,i}^n + \frac{\Delta t}{\rho_e V} \left( \dot m (T_{l,i-1}^n - T_{l,i}^n) - \frac{h_{ta}A}{C_{pe}}(T_{l,i}^n - T_a)  \right)$$

### Final Equations
We have equations of the form
$$\begin{bmatrix}
T^{n+1}_{p,1}\\
T^{n+1}_{p,2}\\
T^{n+1}_{p,3}\\
\vdots\\
T^{n+1}_{p,N_c}\\
T^{n+1}_{f,1}\\
T^{n+1}_{f,2}\\
\vdots\\
T^{n+1}_{f, N_c}\\
T^{n+1}_{l,1}\\
T^{n+1}_{l,2}\\
\vdots\\
T^{n+1}_{l,N_s}
\end{bmatrix} = f(\mathbf{T_p}^n, \mathbf{T_f}^n, \mathbf{T_l}^n)
$$
#### Plate
$$T_{p,i}^{n+1} = T_{p,i}^n + c_1\Delta t\left(S - h_{pf}(T_{p,i}^n - T_{f,i}^n) - h_{pa}(T_{p,i}^n - T_a) - \alpha \left((T_{p,i}^n)^4 - T_{\text{sky}}^4\right)\right)$$
- $c_1 = \frac{1}{\rho_p C_{pp}\delta}$
#### Fluid in Collector
$i=1$
$$T_{f,1}^{n+1} = T_{f,1}^n + c_2\Delta t \left( c_3(T_{p,1}^n - T_{f,1}^n) - \dot m\frac{T_{f,2}^n - T_{l,N_s}^n}{2\Delta y}\right)$$
$i \in 2,\ldots N_c-1$
$$T_{f,i}^{n+1} = T_{f,i}^n + c_2\Delta t \left( c_3(T_{p,i}^n - T_{f,i}^n) - \dot m\frac{T_{f,i+1}^n - T_{l,i-1}^n}{2\Delta y}\right)$$
$i = N_c$
$$T_{f,N_c}^{n+1} = T_{f,N_c}^n + c_2\Delta t \left( c_3(T_{p,N_c}^n - T_{f,N_c}^n) - \dot m\frac{T_{l,1}^n - T_{f,N_c-1}^n}{2\Delta y}\right)$$
- $c_2 = \frac{1}{\rho_f A_f}$
- $c_3 = \frac{Wh_{pf}}{C_{pf}}$
#### Storage Tank
$i=1$
$$T_{l,1}^{n+1} = T_{l,1}^n + c_4\Delta t \left( \dot m (T_{f,N_c}^n - T_{l,1}^n) -c_5(T_{l,1}^n - T_a)  \right)$$
$i \in 2,\ldots N_s$
$$T_{l,i}^{n+1} = T_{l,i}^n + c_4\Delta t \left( \dot m (T_{l,i-1}^n - T_{l,i}^n) - c_5(T_{l,i}^n - T_a)  \right)$$
- $c_4 = \frac{1}{\rho_e V}$
- $c_5 = \frac{h_{ta}A}{C_{pe}}$

## Limitations and Possible Refinements

- **Stability Analysis:** We can perform a e.g. von Neumann stability analysis to pick better values for the step size $\Delta t$ to strike a balance between convergence rate and tolerance while ensuring stability. The explicit scheme used is numerically unstable for larger values of $\Delta t$.
- **Using Implicit Schemes:** The simulation would benefit from using an implicit scheme such as the Crank-Nicolson method. This would require performing a non-linear solve on every iteration using e.g. a Newton-like method due to the radiative term in the plate equations. Alternatively, we can make do with a linear solve if we neglect losses due to radiation, which could be reasonable since the coefficient $\alpha \sim 1.0e-8$, many orders of magnitude smaller than typical convective coefficients.
- **A less complex model(?):** For starters I probably should have stuck with something like a lumped element model, if that is possible for thermal systems. Some model design strategies exist, such as modeling a "thermal circuit" with temperature $\equiv$ voltage and heat flux $\equiv$ electric current. This treatment is presented in [6]
- **A more complex model(?):** Both solar collectors and storage tanks have been modelled with a fair amount of complexity [2,4]. 



##  References
- [1] Incropera, Frank P., et al. Fundamentals of heat and mass transfer. Vol. 6. New York: Wiley, 1996.
- [2] Zeghib, I., & Chaker, A. (2011). Simulation of a solar domestic water heating system. Energy Procedia, 6, 292-301.
- [3] Kalita, K. (2020). Solar Energy Engineering and Technology [MOOC]. NPTEL. https://onlinecourses.nptel.ac.in/noc20_ph14/preview
- [4] Al-Tabbakh, A. A. (2022). Numerical transient modeling of a flat plate solar collector. Results in Engineering, 15, 100580
- [5] Rosales, R. R. (2022) Notes: von Neumann Stability Analysis. [Course Notes] https://math.mit.edu/classes/18.300/Notes/Notes_vNSA.pdf
- [6] Rowell, D., & Wormley, D. N. (1997). System dynamics: an introduction (Vol. 635). Upper Saddle River: Prentice Hall.