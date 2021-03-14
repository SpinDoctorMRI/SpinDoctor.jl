"""
    prepare_pde!(setup)

Prepare PDE compartments.
"""
function prepare_pde!(setup)
    @unpack cell_shape, ncell, include_in, in_ratio, ecs_shape, ecs_ratio = setup.geometry

    include_ecs = ecs_shape != "no_ecs"

    # Check for correct radius ratios and that neurons do not have in-compartments
    @assert !include_in || 0 < in_ratio && in_ratio < 1 && cell_shape != "neuron"
    @assert !include_ecs || 0 < ecs_ratio

    # Determine number of compartments and boundaries
    ncompartment = (1 + include_in) * ncell + include_ecs
    if cell_shape == "sphere"
        # For a sphere, there is one interface
        nboundary = (include_in + 1) * ncell + include_ecs
    elseif cell_shape == "cylinder"
        # An axon has a side interface, and a top-bottom boundary
        nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs
    elseif cell_shape == "neuron"
        # For a neuron, there is one interface
        nboundary = 1 + include_ecs
    end

    boundaries = String[]
    boundary_markers = fill(false, ncompartment, nboundary)
    if include_in
        # Add in-compartments and in-out-interfaces
        compartments = repeat(["in"], ncell)
        boundaries = repeat(["in,out"], ncell)
        boundary_markers[CartesianIndex.(1:2*ncell, repeat(1:ncell, 2))] .= true
        ncompartment_old = ncell
        nboundary_old = ncell
    else
        compartments = []
        boundaries = []
        ncompartment_old = 0
        nboundary_old = 0
    end

    # Add out-compartments
    append!(compartments, repeat(["out"], ncell))

    if include_ecs
        # Add ecs-compartment and out-ecs interfaces
        push!(compartments, "ecs")
        append!(boundaries, repeat(["out,ecs"], ncell))
        boundary_markers[CartesianIndex.(
            [ncompartment_old+1:ncompartment_old+ncell; fill(ncompartment, ncell)],
            repeat(nboundary_old+1:nboundary_old+ncell, 2),
        )] .= true
        nboundary_old = nboundary_old + ncell
    end

    if cell_shape == "cylinder"
        if include_in
            # Add in boundary
            append!(boundaries, repeat(["in"], ncell))
            ncompartment_old = ncell
            boundary_markers[CartesianIndex.(
                1:ncell,
                nboundary_old+1:nboundary_old+ncell,
            )] .= true
            nboundary_old = nboundary_old + ncell
        else
            ncompartment_old = 0
        end

        # Add out boundary
        append!(boundaries, repeat(["out"], ncell))
        boundary_markers[CartesianIndex.(
            ncompartment_old+1:ncompartment_old+ncell,
            nboundary_old+1:nboundary_old+ncell,
        )] .= true

        if include_ecs
            # Add ecs boundary
            push!(boundaries, "ecs")
            boundary_markers[end, end] = true
        end
    elseif include_ecs
        # Add ecs boundary
        push!(boundaries, "ecs")
        boundary_markers[end, end] = true
    else
        # Add out boundary
        push!(boundaries, "out")
        boundary_markers[end, end] = true
    end

    # Initialize output arrays
    σ = zeros(ncompartment)
    T₂ = zeros(ncompartment)
    κ = zeros(nboundary)
    ρ = zeros(ncompartment)

    # Distribute material properties to compartments and boundaries
    ρ[compartments.=="in"] .= setup.pde[:ρ_in]
    ρ[compartments.=="out"] .= setup.pde[:ρ_out]
    ρ[compartments.=="ecs"] .= setup.pde[:ρ_ecs]
    σ[compartments.=="in"] .= setup.pde[:σ_in]
    σ[compartments.=="out"] .= setup.pde[:σ_out]
    σ[compartments.=="ecs"] .= setup.pde[:σ_ecs]
    T₂[compartments.=="in"] .= setup.pde[:T₂_in]
    T₂[compartments.=="out"] .= setup.pde[:T₂_out]
    T₂[compartments.=="ecs"] .= setup.pde[:T₂_ecs]
    κ[boundaries.=="in,out"] .= setup.pde[:κ_in_out]
    κ[boundaries.=="out,ecs"] .= setup.pde[:κ_out_ecs]
    κ[boundaries.=="in"] .= setup.pde[:κ_in]
    κ[boundaries.=="out"] .= setup.pde[:κ_out]
    κ[boundaries.=="ecs"] .= setup.pde[:κ_ecs]

    # Add fields to dictinary
    setup.pde[:ρ] = ρ
    setup.pde[:σ] = σ
    setup.pde[:T₂] = T₂
    setup.pde[:κ] = κ
    setup.pde[:compartments] = compartments
    setup.pde[:boundaries] = boundaries
    setup.pde[:boundary_markers] = boundary_markers

end
