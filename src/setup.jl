@with_kw mutable struct Setup
    name::String
    geometry::Union{Dict{Symbol,Any},Nothing} = nothing
    pde::Union{Dict{Symbol,Any},Nothing} = nothing
    gradient::Union{Dict{Symbol,Any},Nothing} = nothing
    btpde::Union{Dict{Symbol,Any},Nothing} = nothing
    hadc::Union{Dict{Symbol,Any},Nothing} = nothing
    mf::Union{Dict{Symbol,Any},Nothing} = nothing
    analytical::Union{Dict{Symbol,Any},Nothing} = nothing
end

function check_consistency(setup::Setup)
    if setup.shape == "neuron"
        @assert setup.ncell == 1
        @assert setup.deformation[1] ≈ 0.0
        @assert setup.deformation[2] ≈ 0.0
    elseif setup.shape == "cylinder"
        @assert setup.ncell ≥ 1
        @assert setup.height > 0.0
        @assert setup.deformation[1] ≥ 0.0
        @assert setup.deformation[2] ≥ 0.0
        @assert 0.0 < setup.rmin ≤ setup.rmax < Inf
        @assert 0.0 < setup.dmin ≤ setup.dmax < Inf
    elseif setup.shape == "sphere"
        @assert setup.ncell ≥ 1
        @assert setup.deformation[1] ≈ 0.0
        @assert setup.deformation[2] ≈ 0.0
        @assert 0.0 < setup.rmin ≤ setup.rmax < Inf
        @assert 0.0 < setup.dmin ≤ setup.dmax < Inf
    else
        error("""CellSetup shape must be "sphere", "cylinder" or "neuron".""")
    end
    setup.include_in ? (@assert setup.in_ratio > 0.0) : nothing
    @assert setup.ecs_shape ∈ ["no_ecs", "box", "convexhull", "tightwrap"]
end
