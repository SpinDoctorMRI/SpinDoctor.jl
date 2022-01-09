"""
Abstract problem.
"""
abstract type AbstractProblem{T} end

"""
    GeneralBTPDE(;
        model,
        matrices,
        abstol = 1e-6,
        reltol = 1e-4,
        odesolver = QNDF(autodiff = false),
    )

General BTPDE problem.
"""
Base.@kwdef struct GeneralBTPDE{T,S} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    abstol::T = 1e-6
    reltol::T = 1e-4
    odesolver::S = QNDF(autodiff = false)
end

"""
    IntervalConstanBTPDE(; model, matrices, θ = 0.5, timestep)

BTPDE problem specialized on intervalwise constant `ScalarGradient`s, e.g [`PGSE`](@ref),
[`DoublePGSE`](@ref).
"""
Base.@kwdef struct IntervalConstanBTPDE{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    θ::T = 0.5
    timestep::T
end

"""
    HADC(;
        model,
        matrices,
        abstol = 1e-6,
        reltol = 1e-4,
        odesolver = QNDF(autodiff = false),
    )

HADC problem.
"""
Base.@kwdef struct HADC{T,S} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    abstol::T = 1e-6
    reltol::T = 1e-4
    odesolver::S = QNDF(autodiff = false)
end

"""
    Karger(; model, difftensors, odesolver, timestep)

Karger problem.
"""
Base.@kwdef struct Karger{T,S} <: AbstractProblem{T}
    model::Model{T}
    difftensors::Vector{SMatrix{3,3,T,9}}
    odesolver::S = MagnusGL6()
    timestep::T
end

"""
    Laplace(; model, matrices, neig_max)

Laplace eigenvalue problem.
"""
Base.@kwdef struct Laplace{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    neig_max::Int
end

"""
    MatrixFormalism(; model, matrices, lap_eig, ninterval)

Matrix formalism problem. Given a Laplace eigendecomposition `lap_eig`, this problem
consists of computing the MF magnetization.
"""
Base.@kwdef struct MatrixFormalism{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    lap_eig::NamedTuple
    ninterval::Int
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
