"""
Previous attempts that were not particularly well architected.
They weren't correct either, so RIP.
"""

function step!(timeseries_data, prob::SWHSProblem, i)
    @unpack Δt, Δy, Nc, c1, c2, c3, c4, c5 = prob
    @unpack Ns, S, mdot, hpf, hpa, hta, Ta, T∞, α = prob.model
    Tcurr = @view timeseries_data[:,i]
    Tprev = @view timeseries_data[:,i-1]
    
    plate_begin = firstindex(Tcurr)
    plate_end = Nc
    fluid_begin = plate_end+1
    fluid_end = 2Nc
    tank_begin = fluid_end + 1
    tank_end = lastindex(Tcurr)
    
    Tp_curr = @view Tcurr[plate_begin:plate_end]
    Tp_prev = @view Tprev[plate_begin:plate_end]
    Tf_curr = @view Tcurr[fluid_begin:fluid_end]
    Tf_prev = @view Tprev[fluid_begin:fluid_end]
    Tl_curr = @view Tcurr[tank_begin:tank_end]
    Tl_prev = @view Tprev[tank_begin:tank_end]
    @assert length(Tp_curr) == length(Tf_curr) == Nc
    @assert length(Tl_curr) == Ns

    c1Δt = c1*Δt
    c2Δt = c2*Δt
    c4Δt = c4*Δt
    
    # update plate temperature
    Tp_curr .= Tp_prev .+ c1Δt .* (S .- hpf .* (Tp_prev .- Tf_prev) .- hpa .* (Tp_prev .- Ta) .- α .* (Tp_prev .^ 4 .- (T∞^4)))
    # update collector fluid temperature
    Tf_curr[begin] = Tf_prev[begin] + c2Δt *(c3*(Tp_prev[begin] - Tf_prev[begin]) - mdot* (Tf_prev[begin+1] - Tl_prev[end])/(2Δy))
    Tf_curr_m = @view Tf_curr[begin+1:end-1]  # exclude terminal entries
    Tf_prev_m = @view Tf_prev[begin+1:end-1]
    Tp_curr_m = @view Tp_curr[begin+1:end-1]
    Tp_prev_m = @view Tp_prev[begin+1:end-1]
    Tf_prev_m_sl = @view Tf_prev[begin:end-2]  # `_sl` ≡ shifted left
    Tf_prev_m_sr = @view Tf_prev[begin+2:end]  # `_sr` ≡ shifted right
    Tf_curr_m .= Tf_prev_m .+ c2Δt .* (c3.*(Tp_prev_m .- Tf_prev_m) .- mdot .* (Tf_prev_m_sr .- Tf_prev_m_sl)/(2Δy))
    Tf_curr[end] = Tf_prev[end] + c2Δt*(c3*(Tp_prev[end] - Tf_prev[end]) - mdot * (Tl_prev[end] - Tf_prev[end-1])/(2Δy))
    # update tank temperature
    Tl_curr[begin] = Tl_prev[begin] + c4Δt * (mdot * (Tf_prev[end] - Tl_prev[begin]) - c5*(Tl_prev[begin] - Ta))
    Tl_curr_m = @view Tl_curr[begin+1:end]
    Tl_prev_m = @view Tl_prev[begin+1:end]
    Tl_prev_m_sl = @view Tl_prev[begin:end-1]
    Tl_curr_m .= Tl_prev_m .+ c4Δt .* (mdot .* (Tl_prev_m_sl .- Tl_prev_m) - c5 .* (Tl_prev_m .- Ta))
    timeseries_data
end

function solve(prob::SWHSProblem; T=0.1)
    @unpack model, Nc, Tp0, Tf0, Δt = prob
    @unpack Ns = model
    maxiter = convert(Int, T ÷ Δt)
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

function f!(du, u, prob::SWHSProblem)
    @unpack ρf, Af, Cpf, mdot, W, hpf, Ns, Nc, hpa, T∞ = prob.model
    @unpack Cpp, δ, ρp, ρe, Cpe, V, A, hta, S, Ta, α, L = prob.model
    Δy = L / Nc
    # TODO: consider precomputing these values and storing in a cache struct
    plate_begin = firstindex(u)
    plate_end = Nc
    fluid_begin = plate_end+1
    fluid_end = 2Nc
    tank_begin = fluid_end + 1
    tank_end = lastindex(u)

    Tp = @view u[plate_begin:plate_end]
    dTp = @view du[plate_begin:plate_end]
    Tf = @view u[fluid_begin:fluid_end]
    dTf = @view du[fluid_begin:fluid_end]
    @assert length(Tp) == length(Tf) == Nc
    Tl = @view u[tank_begin:tank_end];
    dTl = @view du[tank_begin:tank_end]
    Tf_m = @view u[fluid_begin+1:fluid_end-1]
    Tp_m = @view u[plate_begin+1:plate_end-1]
    @assert length(Tf_m) == length(Tp_m) == Nc-2
    dTf_m = @view du[fluid_begin+1:fluid_end-1]
    Tf_sl = @view u[fluid_begin:fluid_end-2]
    Tf_sr = @view u[fluid_begin+2:fluid_end]
    Tl_r = @view u[tank_begin+1:end]
    dTl_r = @view du[tank_begin+1:end]
    Tl_l = @view u[tank_begin:end-1]

    @assert length(Tl) == Ns

    ρp_δ_Cpp_inv = 1/(ρp*δ*Cpp)
    ρf_Af_Cpf_inv = 1/(ρf * Af * Cpf)
    mdot_Cpf = mdot * Cpf
    W_hpf = W*hpf
    ρe_V_Cpe_inv = 1/(ρe * V * Cpe)
    mdot_Cpe = mdot * Cpe
    hta_A = hta * A
    
    # compute dTₚ/dt
    @. dTp = ρp_δ_Cpp_inv * (S - hpf * (Tp - Tf) - hpa * (Tp - Ta) - α * (Tp^4 - T∞^4))
    # compute dTf/dt
    # i = 1
    dTf[begin] = ρf_Af_Cpf_inv * (W_hpf * (Tp[begin] - Tf[begin]) - mdot_Cpf * (Tf[begin+1] - Tl[end])/(2Δy))
    @. dTf_m = ρf_Af_Cpf_inv * (W_hpf * (Tp_m - Tf_m) - mdot_Cpf * (Tf_sr - Tf_sl)/(2Δy))
    dTf[end] = ρf_Af_Cpf_inv * (W_hpf * (Tp[end] - Tf[end]) - mdot_Cpf * (Tl[begin] - Tf[end-1])/(2Δy))

    # compute dTₗ/dt
    dTl[begin] = ρe_V_Cpe_inv * (mdot_Cpe * (Tf[end] - Tl[begin]) - hta_A*(Tl[begin] - Ta))
    @. dTl_r = ρe_V_Cpe_inv * (mdot_Cpe * (Tl_r - Tl_l) - hta_A*(Tl_r - Ta))
    
    du
end
