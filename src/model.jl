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
    Ns::Int # striations count
    function SWHSModel()
        
    end
end

function semidiscretize(model::SWHSModel; Nc = 100)
    
end
