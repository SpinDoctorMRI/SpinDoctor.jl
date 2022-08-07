"""
Abstract setup recipe.
"""
abstract type AbstractSetup{T} end

"""
Setup recipe for a set of stacked plates.
"""
Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    name::String = "plates"
    meshdir = "meshfiles"
    depth::T
    widths::Vector{T} = [1, 1, 1, 1, 1]
    refinement::T = Inf
end

"""
    DiskSetup

Setup recipe for a set of 2D disks immersed in a ECS.
"""
Base.@kwdef struct DiskSetup{T} <: AbstractSetup{T}
    name::String = "disks"
    meshdir = "meshfiles"
    ncell::Int
    nlayer::Int = 1
    layersizes::Vector{T} = [1]
    rmin::T = 2
    rmax::T = 6
    dmin::T = 2//5
    dmax::T = 3//5
    nsidewall::Int = 12 # Average number of edges of circle defining the cylinder
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 1//2
    refinement::T = Inf
end

"""
    ExtrusionSetup

Setup recipe for an "extruded" 2D geometry, based on a 2D setup `groundsetup`.
"""
Base.@kwdef struct ExtrusionSetup{T,S} <: AbstractSetup{T}
    name::String = "extrusion"
    meshdir::String = "meshfiles"
    groundsetup::S
    refinement::T = Inf
    height::T = 1
    bend::T = 0
    twist::T = 0
end

"""
Setup recipe for a set of spheres immersed in a ECS.
"""
Base.@kwdef struct SphereSetup{T} <: AbstractSetup{T}
    name::String = "spheres"
    meshdir = "meshfiles"
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
    name::String = "neuron"
    meshdir = "meshfiles"
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
end
