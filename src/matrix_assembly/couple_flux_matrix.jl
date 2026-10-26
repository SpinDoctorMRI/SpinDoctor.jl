"""
    couple_flux_matrix(model[, symmetrical])

Assemble the global flux matrix `Q` of the permeable interfaces between compartments and of
the relaxing outer boundaries.

At an interface with permeability `κ`, the coupling coefficients on each side are weighted
by the initial spin densities `ρ₁` and `ρ₂` (`2ρ₂ / (ρ₁ + ρ₂) κ` and `2ρ₁ / (ρ₁ + ρ₂) κ`), so
that the equilibrium magnetization `ρ` has no net flux. With `symmetrical = true`, `κ` is
used on both sides.
"""
function couple_flux_matrix(model, symmetrical = false)
    (; mesh, ρ, κ) = model
    (; point_map, points, facets) = mesh
    ncompartment, nboundary = size(facets)
    offsets = cumsum([0; size.(points, 2)])
    npoint = offsets[end]

    I, J, V = Int[], Int[], eltype(points[1])[]
    for iboundary = 1:nboundary
        cmpts = findall(!isempty, facets[:, iboundary])
        if length(cmpts) == 1
            # Outer boundary, possibly with a boundary relaxation coefficient
            c = only(cmpts)
            iszero(κ[iboundary]) && continue
            f = facets[c, iboundary]
            g = f .+ offsets[c]
            facet_mass_triplets!(I, J, V, points[c], f, g, g, κ[iboundary])
        elseif length(cmpts) == 2
            # Interface: the facets in the local numbering of both compartments
            c₁, c₂ = cmpts
            f₁ = facets[c₁, iboundary]
            f₂ = local_indices(point_map[c₂], point_map[c₁][f₁])
            g₁ = f₁ .+ offsets[c₁]
            g₂ = f₂ .+ offsets[c₂]

            if symmetrical
                κ₁ = κ₂ = κ[iboundary]
            else
                ρ₁, ρ₂ = ρ[c₁], ρ[c₂]
                κ₁ = real(2ρ₂ / (ρ₁ + ρ₂)) * κ[iboundary]
                κ₂ = real(2ρ₁ / (ρ₁ + ρ₂)) * κ[iboundary]
            end

            # Q₁₁ += κ₁ M, Q₁₂ -= κ₂ M, Q₂₁ -= κ₁ M, Q₂₂ += κ₂ M (M: interface mass matrix)
            facet_mass_triplets!(I, J, V, points[c₁], f₁, g₁, g₁, κ₁)
            facet_mass_triplets!(I, J, V, points[c₁], f₁, g₁, g₂, -κ₂)
            facet_mass_triplets!(I, J, V, points[c₁], f₁, g₂, g₁, -κ₁)
            facet_mass_triplets!(I, J, V, points[c₁], f₁, g₂, g₂, κ₂)
        elseif length(cmpts) > 2
            error("Boundary $iboundary touches more than two compartments")
        end
    end

    sparse(I, J, V, npoint, npoint)
end

"""
    local_indices(point_map, globals)

Local indices of the global points `globals` in a compartment with (sorted) `point_map`.
"""
function local_indices(point_map, globals)
    issorted(point_map) || error("The point map of the compartment must be sorted")
    map(globals) do g
        i = searchsortedfirst(point_map, g)
        i ≤ length(point_map) && point_map[i] == g || error("Point $g is not in the compartment")
        i
    end
end
