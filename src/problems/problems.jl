abstract type AbstractProblem{T} end

Base.@kwdef struct GeneralBTPDE{T,S} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    abstol::T
    reltol::T
    odesolver::S
end

Base.@kwdef struct IntervalConstanBTPDE{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    θ::T
    timestep::T
end

Base.@kwdef struct HADC{T,S} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    abstol::T
    reltol::T
    odesolver::S
end

Base.@kwdef struct Karger{T,S} <: AbstractProblem{T}
    model::Model{T}
    difftensors::Vector{SMatrix{3,3,T,9}}
    odesolver::S
    timestep::T
end

Base.@kwdef struct Laplace{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    neig_max::Int
end

Base.@kwdef struct MatrixFormalism{T} <: AbstractProblem{T}
    model::Model{T}
    matrices::NamedTuple
    lap_eig::NamedTuple
    ninterval::Int
end

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

Base.@kwdef struct AnalyticalMatrixFormalism{T} <: AbstractProblem{T}
    analytical_laplace::AnalyticalLaplace{T}
    lap_mat::NamedTuple
    volumes::Vector{T}
end
