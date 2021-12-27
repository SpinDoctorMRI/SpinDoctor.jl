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
    refinement::Union{T,Nothing} = nothing
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
    refinement::Union{T,Nothing} = nothing
end

Base.@kwdef struct PlateSetup{T} <: AbstractSetup{T}
    name::String
    width::T
    depth::T
    heights::Vector{T}
    bend::T
    twist::T
    refinement::Union{T,Nothing} = nothing
end

Base.@kwdef struct NeuronSetup{T} <: AbstractSetup{T}
    name::String
    ncell::Int = 1
    include_in::Bool = false
    in_ratio::T = 0
    ecs_shape::String
    ecs_ratio::T
    refinement::Union{T,Nothing} = nothing
end
