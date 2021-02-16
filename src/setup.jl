@with_kw struct CellSetup
    name::String
    shape::String = "cylinder"
    ncell::Int = 1
    rmin::Union{Float64, Nothing} = nothing
    rmax::Union{Float64, Nothing} = nothing
    dmin::Union{Float64, Nothing} = nothing
    dmax::Union{Float64, Nothing} = nothing
    deformation::Tuple{Float64, Float64} = (0., 0.)
    height::Union{Float64, Nothing} = nothing
    include_nucleus::Bool = false
    nucleus_radiusratio::Float64 = 0.
    include_ecs::Bool = false
    ecs_shape::String = "box"
    ecs_gap::Float64 = 0.
    refinement::Union{Float64, Nothing} = nothing
end

function check_consistency(setup::CellSetup)
    if setup.shape == "neuron"
        @assert setup.ncell == 1
        @assert setup.deformation[1] ≈ 0.
        @assert setup.deformation[2] ≈ 0.
    elseif setup.shape == "cylinder"
        @assert setup.ncell ≥ 1
        @assert  setup.height > 0.
        @assert setup.deformation[1] ≥ 0.
        @assert setup.deformation[2] ≥ 0.
        @assert 0. < setup.rmin ≤ setup.rmax < Inf
        @assert 0. < setup.dmin ≤ setup.dmax < Inf
    elseif setup.shape == "sphere"
        @assert setup.ncell ≥ 1
        @assert setup.deformation[1] ≈ 0.
        @assert setup.deformation[2] ≈ 0.
        @assert 0. < setup.rmin ≤ setup.rmax < Inf
        @assert 0. < setup.dmin ≤ setup.dmax < Inf
    else
        error("""CellSetup shape must be "sphere", "cylinder" or "neuron".""")
    end
    setup.include_nucleus ? (@assert setup.nucleus_radiusratio > 0.) : nothing
    @assert setup.ecs_shape ∈ ["no_ecs", "box", "convexhull", "tightwrap"]
end

@with_kw struct DomainSetup
    diffusivity_in::Float64
    diffusivity_out::Float64
    diffusivity_ecs::Float64
    relaxation_in::Float64 = Inf
    relaxation_out::Float64 = Inf
    relaxation_ecs::Float64 = Inf
    initial_density_in::Float64
    initial_density_out::Float64
    initial_density_ecs::Float64
    permeability_in_out::Float64
    permeability_out_ecs::Float64
end

@with_kw struct BTPDE
    odesolver::DiffEqBase.DEAlgorithm = Trapezoid()
    reltol::Float64 = 1e-04
    abstol::Float64 = 1e-06
    nsave::Int = 1
end

@with_kw struct HADC
    odesolver::DiffEqBase.DEAlgorithm = Trapezoid()
    reltol::Float64 = 1e-04
    abstol::Float64 = 1e-06
end

@with_kw struct MF
    length_scale::Float64 = 0.
    neig_max::Int = Inf
    ninterval::Int = 100
end

@with_kw struct ExperimentSetup
    ndirection::Int = 1
    flat_dirs::Bool = false
    direction::Array{Float64} = [1.; 0.; 0.]
    sequences::Array{TimeProfile}
    values::Array{Float64}
    values_type::Char
    btpde::Union{BTPDE, Nothing} = nothing
    hadc::Union{HADC, Nothing} = nothing
    mf::Union{MF, Nothing} = nothing
end
