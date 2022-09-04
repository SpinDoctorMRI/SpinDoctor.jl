"""
    AbstractSetup

Abstract setup recipe.
"""
abstract type AbstractSetup{T} end

"""
    PlateSetup

Setup recipe for a set of stacked plates.
"""
@with_kw struct PlateSetup{T} <: AbstractSetup{T}
    depth::T = 1.0
    @assert depth > 0
    widths::Vector{T} = [1.0, 1.0, 1.0, 1.0, 1.0]
    @assert all(>(0), widths)
    refinement::T = Inf
    @assert refinement > 0
end

"""
    DiskSetup

Setup recipe for a set of 2D disks immersed in a ECS.
"""
@with_kw struct DiskSetup{T,E} <: AbstractSetup{T}
    ncell::Int = 3
    @assert ncell > 0
    layersizes::Vector{T} = [1.0]
    @assert all(s -> 0 < s ≤ 1, layersizes)
    @assert issorted(layersizes)
    rmin::T = 2.0
    rmax::T = 6.0
    @assert 0 < rmin ≤ rmax
    dmin::T = 0.2
    dmax::T = 0.6
    @assert 0 < dmin ≤ dmax
    nsidewall::Int = 12
    @assert nsidewall > 0
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    SphereSetup

Setup recipe for a set of spheres immersed in a ECS.
"""
@with_kw struct SphereSetup{T,E} <: AbstractSetup{T}
    ncell::Int = 3
    @assert ncell > 0
    layersizes::Vector{T} = [1.0]
    @assert all(s -> 0 < s ≤ 1, layersizes)
    @assert issorted(layersizes)
    rmin::T = 2.0
    rmax::T = 6.0
    @assert 0 < rmin ≤ rmax
    dmin::T = 0.2
    dmax::T = 0.6
    @assert 0 < dmin ≤ dmax
    nsidewall::Int = 200
    @assert nsidewall > 0
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    NeuronSetup

Setup recipe for a neuron, possibly wrapped in a ECS.
"""
@with_kw struct NeuronSetup{T,E} <: AbstractSetup{T}
    ecs::E = NoECS()
    refinement::T = Inf
    @assert refinement > 0
end

"""
    ExtrusionSetup

Setup recipe for an "extruded" 2D geometry, based on a 2D setup `groundsetup`.
"""
@with_kw struct ExtrusionSetup{T,S} <: AbstractSetup{T}
    groundsetup::S
    refinement::T = Inf
    @assert refinement > 0
    height::T = 1.0
    @assert height > 0
    bend::T = 0.0
    twist::T = 0.0
end

"""
    CylinderSetup(;
        ncell = 1,
        layersizes = [1.0],
        rmin = 2.0,
        rmax = 6.0,
        dmin = 0.2,
        dmax = 0.6,
        nsidewall = 12,
        ecs = NoECS(),
        height = 1,
        bend = 0.0,
        twist = 0.0,
        refinement = Inf,
    )

Cylinder setup.
"""
CylinderSetup{T,E} = ExtrusionSetup{T,DiskSetup{T,E}}
CylinderSetup(;
    ncell = 1,
    layersizes = [1.0],
    rmin = 2.0,
    rmax = 6.0,
    dmin = 0.2,
    dmax = 0.6,
    nsidewall = 12,
    ecs = NoECS(),
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
    refinement = Inf,
) = ExtrusionSetup(;
    groundsetup = DiskSetup(; ncell, layersizes, rmin, rmax, dmin, dmax, nsidewall, ecs),
    height,
    bend,
    twist,
    refinement,
)

"""
    SlabSetup

Slab setup.
"""
SlabSetup{T} = ExtrusionSetup{T,PlateSetup{T}}
SlabSetup(;
    widths = [1.0, 1.0, 1.0, 1.0, 1.0],
    depth = 5.0,
    height = 1.0,
    bend = 0.0,
    twist = 0.0,
    refinement = Inf,
) = ExtrusionSetup(;
    groundsetup = PlateSetup(; widths, depth),
    height,
    bend,
    twist,
    refinement,
)
