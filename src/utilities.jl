using LinearAlgebra: Tridiagonal

function plate_FD1_matrix(Δy::T, Nc) where {T <: Number}
    a = 1/(2Δy)
    dl = fill(a, Nc-1)
    dl[end] = 0
    du = fill!(-a, Nc-1)
    du[begin] = 0
    d = zeros(T,Nc)
    Tridiagonal(dl, d, du)
end

function plate_FD2_matrix(Δy::T, Nc) where {T<:Number}
    a = 1/(Δy^2)
    dl = fill(a, Nc-1)
    dl[end] = 2a
    du = fill(a, Nc-1)
    du[begin] = 2a
    d = fill(-2a, Nc)
    Tridiagonal(dl, d, du)
end

function fluid_FD1_matrix(Δy::T, Nc) where {T<:Number}
    a = 1/(2Δy)
    d = zeros(T, Nc)
    dl = fill(a, Nc-1)
    dl[end] = -a
    du = fill(-a, Nc-1)
    du[begin] = 0
    Tridiagonal(dl, d, du)
end
