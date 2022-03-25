"""
Abstract setup recipe.
"""
abstract type AbstractSetup{T} end

"""
Setup recipe for a set of stacked plates.
"""
Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    name::String
    meshdir = "meshfiles"
    width::T
    depth::T
    heights::Vector{T}
    deform_angle::NamedTuple{(:bend, :twist), Tuple{T, T}} # [bend; twist]
    refinement::T = Inf
    tetgen_option::String = ""
    D::Vector{Matrix{T}}
    ρ::Vector{Complex{T}}
    T₂::Vector{T}
    κ::NamedTuple{(:interfaces, :boundaries), Tuple{Vector{T}, Vector{T}}}
    κ_interfaces::Vector{T}
    κ_boundaries::Vector{T}
    γ::T = 2.6752218744e-7
end

"""
Setup recipe for a set of cylinders immersed in a ECS.
"""
# add inner constructor to make sure the checksum is the same when include_in = false or ecs=:no_ecs
@with_kw struct CylinderSetup{T} <: AbstractSetup{T}
    name::String
    meshdir = "meshfiles"
    ncell::Int
    r_range::NamedTuple{(:rmin, :rmax), Tuple{T, T}} # (; rmin, rmax)
    d_range::NamedTuple{(:dmin, :dmax), Tuple{T, T}} # (; dmin, dmax)
    height::T
    deform_angle::NamedTuple{(:bend, :twist), Tuple{T, T}} # (; bend, twist)
    include_in::Bool = false
    in_ratio::T = 0.5
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
    tetgen_option::String = ""
    D::NamedTuple{(:in, :out, :ecs), Tuple{Matrix{T}, Matrix{T}, Matrix{T} } } # (; in, out, ecs)
    ρ::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=1.0, out=1.0, ecs=1.0)  # (; in, out, ecs)
    T₂::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=Inf, out=Inf, ecs=Inf)  # (; in, out, ecs)
    κ::NamedTuple{(:in_out, :out_ecs, :in, :out, :ecs), Tuple{T, T, T, T, T}} = (; in_out=1e-4, out_ecs=1e-4, in=0.0, out=0.0, ecs=0.0) # (; in_out, out_ecs, in, out, ecs)
    γ::T = 2.6752218744e-7

    CylinderSetup{T}(name,meshdir,ncell,r_range,d_range, height,
        deform_angle,include_in,in_ratio,ecs_shape,ecs_ratio,
        refinement,tetgen_option,D,ρ,T₂,κ,γ ) where{T} = ( 
            @assert ecs_shape∈[:no_ecs, :box, :convex_hull, :tight_wrap];
            new(name,meshdir,ncell,r_range,d_range,height,
                deform_angle,include_in,include_in*in_ratio,ecs_shape,(ecs_shape != :no_ecs)*ecs_ratio,
                refinement,tetgen_option,D,ρ,T₂,κ,γ)
    )
end

"""
Setup recipe for a set of spheres immersed in a ECS.
"""
Base.@kwdef struct SphereSetup{T} <: AbstractSetup{T}
    name::String
    meshdir = "meshfiles"
    ncell::Int
    r_range::NamedTuple{(:rmin, :rmax), Tuple{T, T}} # (; rmin, rmax)
    d_range::NamedTuple{(:dmin, :dmax), Tuple{T, T}} # (; dmin, dmax)
    include_in::Bool = false
    in_ratio::T = 0.5
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
    tetgen_option::String = ""
    D::NamedTuple{(:in, :out, :ecs), Tuple{Matrix{T}, Matrix{T}, Matrix{T} }} # (; in, out, ecs)
    ρ::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=1.0, out=1.0, ecs=1.0)  # (; in, out, ecs)
    T₂::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=Inf, out=Inf, ecs=Inf)  # (; in, out, ecs)
    κ::NamedTuple{(:in_out, :out_ecs, :in, :out, :ecs), Tuple{T, T, T, T, T}} = (; in_out=1e-4, out_ecs=1e-4, in=0.0, out=0.0, ecs=0.0) # (; in_out, out_ecs, in, out, ecs)
    γ::T = 2.6752218744e-7

    SphereSetup{T}(name,meshdir,ncell,r_range,d_range,
        include_in,in_ratio,ecs_shape,ecs_ratio,
        refinement,tetgen_option,D,ρ,T₂,κ,γ) where{T} = (
            @assert ecs_shape∈[:no_ecs, :box, :convex_hull, :tight_wrap];
            new(name,meshdir,ncell,r_range,d_range,
            include_in,include_in*in_ratio,ecs_shape,(ecs_shape != :no_ecs)*ecs_ratio,
            refinement,tetgen_option,D,ρ,T₂,κ,γ) 
    )
end

"""
Setup recipe for a neuron, possibly wrapped in a ECS.
"""
Base.@kwdef struct NeuronSetup{T} <: AbstractSetup{T}
    name::String
    meshdir = "meshfiles"
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.5
    refinement::T = Inf
    tetgen_option::String = ""
    D::NamedTuple{(:neuron, :ecs), Tuple{Matrix{T}, Matrix{T} } } # (; neuron, ecs)
    ρ::NamedTuple{(:neuron, :ecs), Tuple{T, T}} = (; neuron=1.0, ecs=1.0)  # (; neuron, ecs)
    T₂::NamedTuple{(:neuron, :ecs), Tuple{T, T}} = (; neuron=Inf, ecs=Inf)  # (; neuron, ecs)
    κ::NamedTuple{(:neuron_ecs, :neuron, :ecs), Tuple{T, T, T}} = (; neuron_ecs=1e-4, neuron=0.0, ecs=0.0) # (; neuron_ecs, neuron, ecs)
    γ::T = 2.6752218744e-7

    NeuronSetup{T}(name,meshdir,ecs_shape,ecs_ratio,
        refinement,tetgen_option,D,ρ,T₂,κ,γ ) where {T} = (@assert ecs_shape∈[:no_ecs, :box, :convex_hull, :tight_wrap];
        new(name,meshdir,ecs_shape,(ecs_shape != :no_ecs)*ecs_ratio,refinement,tetgen_option,D,ρ,T₂,κ,γ)
    )
end

"""
Setup recipe for a set of fibers wrapped in a ECS.
"""
@with_kw struct FiberSetup{T} <: AbstractSetup{T}
    name::String
    meshdir::String = "meshfiles"
    ncell::Int
    deform_angle::NamedTuple{(:bend, :twist), Tuple{T, T}} # (; bend, twist)
    include_in::Bool = false
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.0
    refinement::T = Inf
    tetgen_option::String = ""
    D::NamedTuple{(:in, :out, :ecs), Tuple{Matrix{T}, Matrix{T}, Matrix{T} }} # (; in, out, ecs)
    ρ::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=1.0, out=1.0, ecs=1.0)  # (; in, out, ecs)
    T₂::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=Inf, out=Inf, ecs=Inf)  # (; in, out, ecs)
    κ::NamedTuple{(:in_out, :out_ecs, :in, :out, :ecs), Tuple{T, T, T, T, T}} = (; in_out=1e-4, out_ecs=1e-4, in=0.0, out=0.0, ecs=0.0) # (; in_out, out_ecs, in, out, ecs)
    γ::T = 2.6752218744e-7

    FiberSetup{T}(name,meshdir,ncell,deform_angle,include_in,ecs_shape,ecs_ratio,
        refinement,tetgen_option,D,ρ,T₂,κ,γ ) where{T} = ( 
            @assert ecs_shape ∈ [:no_ecs, :box, :convex_hull, :tight_wrap];
            @assert !include_in;
            new(name,meshdir,ncell,deform_angle,include_in,
                ecs_shape,(ecs_shape != :no_ecs)*ecs_ratio,refinement,tetgen_option,D,ρ,T₂,κ,γ)
    )
end

"""
Setup recipe for one EM image wrapped in a ECS.
"""
@with_kw struct EMSetup{T} <: AbstractSetup{T}
    name::String
    meshdir::String = "meshfiles"
    ncell::Int
    height::T
    deform_angle::NamedTuple{(:bend, :twist), Tuple{T, T}} # (; bend, twist)
    include_in::Bool = false
    ecs_shape::Symbol = :no_ecs
    ecs_ratio::T = 0.0
    refinement::T = Inf
    tetgen_option::String = ""
    resolution::T
    D::NamedTuple{(:in, :out, :ecs), Tuple{Matrix{T}, Matrix{T}, Matrix{T} }} # (; in, out, ecs)
    ρ::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=1.0, out=1.0, ecs=1.0)  # (; in, out, ecs)
    T₂::NamedTuple{(:in, :out, :ecs), Tuple{T, T, T}} = (; in=Inf, out=Inf, ecs=Inf)  # (; in, out, ecs)
    κ::NamedTuple{(:in_out, :out_ecs, :in, :out, :ecs), Tuple{T, T, T, T, T}} = (; in_out=1e-4, out_ecs=1e-4, in=0.0, out=0.0, ecs=0.0) # (; in_out, out_ecs, in, out, ecs)
    γ::T = 2.6752218744e-7

    EMSetup{T}(name,meshdir,ncell,height,deform_angle,include_in,ecs_shape,ecs_ratio,
        refinement,tetgen_option,resolution,D,ρ,T₂,κ,γ ) where{T} = ( 
            @assert ecs_shape ∈ [:no_ecs, :box, :convex_hull, :tight_wrap];
            new(name,meshdir,ncell,height,deform_angle,include_in,ecs_shape,ecs_ratio,
                refinement,tetgen_option,resolution,D,ρ,T₂,κ,γ)
    )
end