"""
Abstract problem.
"""
abstract type AbstractProblem{T} end

"""
    BTPDE(; model, matrices)

Bloch-Torrey PDE problem.
"""
Base.@kwdef struct BTPDE{T,N} <: AbstractProblem{T}
    model::Model{T,N}
    matrices::NamedTuple
end

"""
    HADC(; model, matrices)

HADC problem.
"""
Base.@kwdef struct HADC{T,N} <: AbstractProblem{T}
    model::Model{T,N}
    matrices::NamedTuple
end

"""
    Karger(; model, difftensors)

Karger problem.
"""
Base.@kwdef struct Karger{T,N} <: AbstractProblem{T}
    model::Model{T,N}
    difftensors::Vector{Matrix{T}}
end

"""
    Laplace(; model, matrices, neig_max)

Laplace eigenvalue problem.
"""
Base.@kwdef struct Laplace{T,N} <: AbstractProblem{T}
    model::Model{T,N}
    matrices::NamedTuple
    neig_max::Int
end

"""
    MatrixFormalism(; model, matrices, lap_eig)

Matrix formalism problem. Given a Laplace eigendecomposition `lap_eig`, this problem
consists of computing the MF magnetization.
"""
Base.@kwdef struct MatrixFormalism{T,N} <: AbstractProblem{T}
    model::Model{T,N}
    matrices::NamedTuple
    lap_eig::NamedTuple
end

"""
    AnalyticalLaplace(; ρ, r, D, W, T₂, γ, dim, eiglim, eigstep)

Analytical radial Laplace eigenvalue problem.
"""
Base.@kwdef struct AnalyticalLaplace{T} <: AbstractProblem{T}
    ρ::Vector{Complex{T}}
    r::Vector{T}
    D::Vector{T}
    W::Vector{T}
    T₂::Vector{T}
    γ::T
    dim::Int
    eiglim::T
    eigstep::T
end

"""
    AnalyticalMatrixFormalism(; analytical_laplace, lap_mat, volumes)

Analytical Matrix formalism problem. Given a radial Laplace eigendecomposition
`analytical_laplace`, this problem consists of computing the MF compartment signals.
"""
Base.@kwdef struct AnalyticalMatrixFormalism{T} <: AbstractProblem{T}
    analytical_laplace::AnalyticalLaplace{T}
    lap_mat::NamedTuple
    volumes::Vector{T}
end
