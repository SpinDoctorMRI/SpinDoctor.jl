abstract type AbstractSetup{T} end

Base.@kwdef struct SphereSetup{T} <: AbstractSetup{T}
    name::String
    ncell::Int
    rmin::T
    rmax::T
    dmin::T
    dmax::T
    include_in::Bool
    in_ratio::T
    ecs_shape::String
    ecs_ratio::T
    refinement::T = Inf
end

Base.@kwdef struct CylinderSetup{T} <: AbstractSetup{T}
    name::String
    ncell::Int
    rmin::T
    rmax::T
    dmin::T
    dmax::T
    height::T
    bend::T
    twist::T
    include_in::Bool
    in_ratio::T
    ecs_shape::String
    ecs_ratio::T
    refinement::T = Inf
end

Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    name::String
    width::T
    depth::T
    heights::Vector{T}
    bend::T
    twist::T
    refinement::T = Inf
end

Base.@kwdef struct NeuronSetup{T} <: AbstractSetup{T}
    name::String
    ecs_shape::Symbol
    ecs_ratio::T
    refinement::T = Inf
end
