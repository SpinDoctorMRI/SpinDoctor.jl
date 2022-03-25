"""
    couple_flux_matrix(model, Q_blocks[, symmetrical])

Create global flux matrix with coupling between compartments.
"""
function couple_flux_matrix(model, Q_blocks, symmetrical = false)

    # Extract mesh fields
    (; mesh, ρ, κ) = model
    (; point_map, points, facets) = mesh

    # Sizes
    nboundary = size(facets, 2)
    inds_cmpts = cumsum([0; size.(points, 2)])
    npoint = inds_cmpts[end]
    get_inds(icmpt) = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]


    # Assemble sparse block diagonal flux matrix
    Q = spzeros(npoint, npoint)

    # Couple flux matrix at boundaries
    for iboundary = 1:nboundary

        # Check if boundary is an interface between two compartments
        cmpts_touch = findall(.!isempty.(facets[:, iboundary]))
        ntouch = length(cmpts_touch)
        if ntouch == 1
            # Only one compartment touches boundary. It is thus an outer boundary,
            # and may possibly have a boundary relaxation coefficient
            cmpt = cmpts_touch[1]

            # Global indices of the boundary
            inds = get_inds(cmpt)

            # Add boundary contribution to global flux matrix for compartment
            Q[inds, inds] += κ[iboundary] * Q_blocks[cmpt, iboundary]
        elseif ntouch == 2
            # Extract flux matrices from corresponding compartments
            cmpt₁, cmpt₂ = cmpts_touch[1], cmpts_touch[2]
            Q₁₁ = Q_blocks[cmpt₁, iboundary]
            Q₂₂ = Q_blocks[cmpt₂, iboundary]
            Q₁₂ = spzeros(size(Q₁₁, 1), size(Q₂₂, 2))
            Q₂₁ = spzeros(size(Q₂₂, 1), size(Q₁₁, 2))

            # Identify pairs of points using global numbering
            inds₁ = unique(facets[cmpt₁, iboundary])
            inds₂ = unique(facets[cmpt₂, iboundary])

            # Check if local indices are already sorted with a one-to-one correspondence
            if all(point_map[cmpt₁][inds₁] .== point_map[cmpt₂][inds₂])
                indinds₁ = 1:length(inds₁)
                indinds₂ = 1:length(inds₂)
            else
                pairmap = point_map[cmpt₁][inds₁] .== point_map[cmpt₂][inds₂]'
                indinds₁ = 1:length(inds₁)
                indinds₂ = [findfirst(pairmap[i, :]) for i ∈ indinds₁]
            end

            # Create coupling blocks
            Q₁₂[:, inds₂[indinds₂]] = Q₁₁[:, inds₁[indinds₁]]
            Q₂₁[:, inds₁[indinds₁]] = Q₂₂[:, inds₂[indinds₂]]

            # Check whether to use same permeability coefficients on both sides
            if symmetrical
                # Use the same permeability coefficient on each side of the boundary
                c₁₂ = 1
                c₂₁ = 1
            else
                # Weigh permeability coefficients on each side of the boundary with
                # initial spin density equilibrium
                ρ₁ = ρ[cmpt₁]
                ρ₂ = ρ[cmpt₂]
                c₂₁ = 2 * ρ₂ / (ρ₁ + ρ₂)
                c₁₂ = 2 * ρ₁ / (ρ₁ + ρ₂)
            end

            # Adjust permeability coefficients
            κ₁ = c₂₁ * κ[iboundary]
            κ₂ = c₁₂ * κ[iboundary]

            # Global indices of the boundary in each of the compartments
            inds₁ = get_inds(cmpt₁)
            inds₂ = get_inds(cmpt₂)

            # Add interface contribution to global flux matrix for compartment 1
            Q[inds₁, inds₁] += κ₁ * Q₁₁
            Q[inds₁, inds₂] -= κ₂ * Q₁₂

            # Add interface contribution to global flux matrix for compartment 2
            Q[inds₂, inds₁] -= κ₁ * Q₂₁
            Q[inds₂, inds₂] += κ₂ * Q₂₂
        else
            error("Each interface only touches 1 or 2 compartments")
        end
    end

    Q
end
