"""
    assemble_matrices(model)

Assemble finite element matrices.

Returns a named tuple with the global (block diagonal) mass matrix `M`, stiffness matrix
`S`, relaxation matrix `R` (`M` divided by `T₂`), first order moment matrices `Mx`
(`∫ xᵢ φⱼ φₖ`, one per axis) and flux matrix `Q` (permeable interfaces and outer
boundaries), the compartment blocks `M_cmpts`, `S_cmpts` and `Mx_cmpts`, the surface
integrals `G` (`∫_Γ φᵢ n dS` over the boundary `Γ` of each compartment, one row per point and
one column per axis) and the compartment `volumes`.
"""
function assemble_matrices(model::Model{T,dim}) where {T,dim}
    (; mesh, D, T₂) = model
    ncompartment = length(mesh.points)

    M_cmpts = SparseMatrixCSC{T,Int}[]
    S_cmpts = SparseMatrixCSC{T,Int}[]
    Mx_cmpts = [SparseMatrixCSC{T,Int}[] for _ = 1:dim]
    G = Matrix{T}[]
    volumes = zeros(T, ncompartment)
    for icmpt = 1:ncompartment
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        pattern = FEPattern(elements, size(points, 2))
        volumes[icmpt] = sum(first(get_mesh_volumes(points, elements)))

        push!(M_cmpts, assemble_mass_matrix(elements, points; pattern))
        push!(S_cmpts, assemble_stiffness_matrix(elements, points, D[icmpt]; pattern))
        for d = 1:dim
            Mx = assemble_mass_matrix(elements, points; nodeweights = points[d, :], pattern)
            push!(Mx_cmpts[d], Mx)
        end

        # Surface integrals ∫ φᵢ n dS with outward normals n
        facets = reduce(hcat, mesh.facets[icmpt, :])
        areas, _, normals = get_mesh_surfacenormals(points, elements, facets)
        Gᵢ = zeros(T, size(points, 2), dim)
        for f ∈ axes(facets, 2), a ∈ axes(facets, 1), d = 1:dim
            Gᵢ[facets[a, f], d] += areas[f] / size(facets, 1) * normals[d, f]
        end
        push!(G, Gᵢ)
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T₂)...)
    Mx = [blockdiag(Mx_cmpts[d]...) for d = 1:dim]
    Q = couple_flux_matrix(model)

    (; M, S, R, Mx, Q, M_cmpts, S_cmpts, Mx_cmpts, G, volumes)
end

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
        i ≤ length(point_map) && point_map[i] == g ||
            error("Point $g is not in the compartment")
        i
    end
end
