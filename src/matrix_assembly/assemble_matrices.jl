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
