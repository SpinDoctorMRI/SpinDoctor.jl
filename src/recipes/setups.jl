"""
Abstract setup recipe.
"""
abstract type AbstractSetup{T} end

"""
Setup recipe for a set of stacked plates.
"""
Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    name::String
    width::T
    depth::T
    heights::Vector{T}
    bend::T = 0
    twist::T = 0
    refinement::T = Inf
end

"""
Setup recipe for a set of cylinders immersed in a ECS.
"""
Base.@kwdef struct CylinderSetup{T} <: AbstractSetup{T}
    name::String
    ncell::Int
    rmin::T
    rmax::T
    dmin::T
    dmax::T
    height::T
    bend::T = 0
    twist::T = 0
    include_in::Bool = false
    in_ratio::T = 0.5
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
end

"""
Setup recipe for a set of spheres immersed in a ECS.
"""
Base.@kwdef struct SphereSetup{T} <: AbstractSetup{T}
    name::String
    ncell::Int
    rmin::T
    rmax::T
    dmin::T
    dmax::T
    include_in::Bool = false
    in_ratio::T = 0.5
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
end

"""
Setup recipe for a neuron, possibly wrapped in a ECS.
"""
Base.@kwdef struct NeuronSetup{T} <: AbstractSetup{T}
    name::String
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
end
