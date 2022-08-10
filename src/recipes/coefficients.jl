"""
    coefficients(setup; D, T₂, ρ, κ, γ)

Order coefficients compartment arrays.
"""
function coefficients end

function coefficients(setup::PlateSetup{T}; D, T₂, ρ, κ, γ) where {T}
    n = length(setup.widths)
    coeffs = (;
        D = [zeros(T, 0, 0) for _ = 1:n],
        T₂ = zeros(T, n),
        κ = zeros(T, n * (n + 1) ÷ 2),
        ρ = zeros(Complex{T}, n),
        γ = T(γ),
    )
    for i = 1:n-1
        coeffs.κ[(i-1)*(2n-i)÷2+1] = κ.interfaces[i]
    end
    for i = 1:n
        coeffs.D[i] = D[i]
        coeffs.T₂[i] = T₂[i]
        coeffs.ρ[i] = ρ[i]
        coeffs.κ[(n-1)*n÷2+i] = κ.boundaries[i]
    end
    coeffs
end

function coefficients(setup::NeuronSetup{T}; D, T₂, ρ, κ, γ) where {T}
    (; ecs) = setup

    include_ecs = !(ecs isa NoECS)

    # Determine number of compartments and boundaries
    ncompartment = 1 + include_ecs
    nboundary = 1 + include_ecs
    if include_ecs
        compartments = ["neuron", "ecs"]
        boundaries = ["neuron,ecs", "ecs"]
    else
        compartments = ["neuron"]
        boundaries = ["neuron"]
    end

    # Initialize output arrays
    coeffs = (;
        D = [zeros(T, 3, 3) for _ = 1:ncompartment],
        T₂ = zeros(T, ncompartment),
        κ = zeros(T, nboundary),
        ρ = zeros(Complex{T}, ncompartment),
        γ = T(γ),
    )

    # Distribute material properties to compartments and boundaries
    coeffs.ρ[compartments.=="neuron"] .= ρ.neuron
    coeffs.ρ[compartments.=="ecs"] .= ρ.ecs
    coeffs.D[compartments.=="neuron"] .= [D.neuron]
    coeffs.D[compartments.=="ecs"] .= [D.ecs]
    coeffs.T₂[compartments.=="neuron"] .= T₂.neuron
    coeffs.T₂[compartments.=="ecs"] .= T₂.ecs
    coeffs.κ[boundaries.=="neuron,ecs"] .= κ.neuron_ecs
    coeffs.κ[boundaries.=="neuron"] .= κ.neuron
    coeffs.κ[boundaries.=="ecs"] .= κ.ecs

    coeffs
end

function coefficients(setup::ExtrusionSetup; D, T₂, ρ, κ, γ)
    (; groundsetup) = setup
    coefficients(groundsetup; D, T₂, ρ, κ, γ)
end

function coefficients(setup::Union{DiskSetup{T},SphereSetup{T}}; D, T₂, ρ, κ, γ) where {T}
    (; ncell, layersizes, ecs) = setup

    # Each disk cell may have multiple layers
    include_ecs = !(ecs isa NoECS)
    nlayer = length(layersizes)
    ncompartment = nlayer * ncell + include_ecs

    # All possible interfaces (n-1)n/2 + all outer boundaries n
    nboundary = ncompartment * (ncompartment + 1) ÷ 2

    # Initialize output arrays
    coeffs = (;
        D = [zeros(T, size(D.cell[1])) for _ = 1:ncompartment],
        T₂ = zeros(T, ncompartment),
        κ = zeros(T, nboundary),
        ρ = zeros(Complex{T}, ncompartment),
        γ = T(γ),
    )

    edgenumber(i, j) = (i - 1) * (2 * ncompartment - i) ÷ 2 + (j - i)

    # Distribute material properties to compartments and boundaries
    k = 1
    for i = 1:nlayer, j = 1:ncell
        coeffs.ρ[k] = ρ.cell[i]
        coeffs.D[k] .= D.cell[i]
        coeffs.T₂[k] = T₂.cell[i]

        # Interfaces
        if i < nlayer
            # (i, j) is linked to (i + 1, j)
            coeffs.κ[edgenumber(k, k + ncell)] = κ.cell_interfaces[i]
        elseif include_ecs
            # The outer layer is linked to the ECS (last region)
            coeffs.κ[edgenumber(k, ncompartment)] = κ.cell_ecs
        end

        # Outer boundaries Γ_k: may be of measure zero in the 2D case, but
        # still assigned in the case of an extrusion setup
        coeffs.κ[(ncompartment-1)*ncompartment÷2+k] = κ.cell_boundaries[i]

        k += 1
    end

    # ECS
    if include_ecs
        coeffs.ρ[end] = ρ.ecs
        coeffs.D[end] .= D.ecs
        coeffs.T₂[end] = T₂.ecs
        coeffs.κ[end] = κ.ecs
    end

    coeffs
end
