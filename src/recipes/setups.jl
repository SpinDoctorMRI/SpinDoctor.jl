"""
Abstract setup recipe.
"""
abstract type AbstractSetup{T} end

"""
Setup recipe for a set of stacked plates.
"""
Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    depth::T = 1
    widths::Vector{T} = [1, 1, 1, 1, 1]
    refinement::T = Inf
end

"""
    DiskSetup

Setup recipe for a set of 2D disks immersed in a ECS.
"""
Base.@kwdef struct DiskSetup{T} <: AbstractSetup{T}
    ncell::Int
    layersizes::Vector{T} = [1]
    rmin::T = 2
    rmax::T = 6
    dmin::T = 2 // 5
    dmax::T = 3 // 5
    nsidewall::Int = 12 # Average number of edges of circle defining the cylinder
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 1 // 2
    refinement::T = Inf
end

"""
Setup recipe for a set of spheres immersed in a ECS.
"""
Base.@kwdef struct SphereSetup{T} <: AbstractSetup{T}
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
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
end

"""
    ExtrusionSetup

Setup recipe for an "extruded" 2D geometry, based on a 2D setup `groundsetup`.
"""
Base.@kwdef struct ExtrusionSetup{T,S} <: AbstractSetup{T}
    groundsetup::S
    refinement::T = Inf
    height::T = 1
    bend::T = 0
    twist::T = 0
end

CylinderSetup{T} = ExtrusionSetup{T,DiskSetup{T}}
CylinderSetup{T}(;
    ncell = 1,
    layersizes = [1],
    rmin = 2,
    rmax = 6,
    dmin = 2 // 5,
    dmax = 3 // 5,
    nsidewall = 12,
    ecs_shape = :no_ecs,
    ecs_ratio = 1 // 2,
    height = 1,
    bend = 0,
    twist = 0,
    refinement = Inf,
) where {T} = ExtrusionSetup(;
    groundsetup = DiskSetup{T}(;
        ncell,
        layersizes,
        rmin,
        rmax,
        dmin,
        dmax,
        nsidewall,
        ecs_shape,
        ecs_ratio,
    ),
    height = T(height),
    bend = T(bend),
    twist = T(twist),
    refinement = T(refinement),
)

SlabSetup{T} = ExtrusionSetup{T,PlateSetup{T}}
SlabSetup{T}(;
    widths = [1, 1, 1, 1, 1],
    depth = 5,
    height = 1,
    bend = 0,
    twist = 0,
    refinement = Inf,
) where {T} = ExtrusionSetup(;
    groundsetup = PlateSetup{T}(;
        widths,
        depth,
    ),
    height = T(height),
    bend = T(bend),
    twist = T(twist),
    refinement = T(refinement),
)
