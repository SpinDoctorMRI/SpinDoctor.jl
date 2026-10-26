# Finite element matrices of a model

"""
    assemble_matrices(model; periodic = false)

Assemble the finite element matrices of `model`.

Returns a named tuple with

- the global (block diagonal) mass matrix `M`, stiffness matrix `S` (with the diffusion
  tensors), relaxation matrix `R` (`M` divided by `T2`), first order moment matrices `Mx`
  (``\\int x_d \\varphi_i \\varphi_j``, one per axis `d`) and flux matrix `Q` (permeable
  interfaces and outer boundaries, see [`assemble_flux_matrix`](@ref));
- the compartment blocks `M_cmpts`, `S_cmpts` and `Mx_cmpts`;
- the surface integrals `G`, one matrix per compartment with the integrals
  ``\\int_{\\partial \\Omega} \\varphi_i n_d \\, \\mathrm{d}s`` over the boundary of the
  compartment with outward normal `n` (one row per point, one column per axis);
- the compartment `volumes`.

With `periodic = true` (or one `Bool` per axis), the mesh is a sample of periodic tissue
in its bounding box, such as the meshes of [`PeriodicBoxECS`](@ref) and
[`VoronoiSetup`](@ref), and each compartment is periodic along the given axes: the points
on opposite faces of the box are identified. The solvers then use periodic boundary
conditions instead of impermeable box faces. The matrices are still those of the whole
mesh, and the named tuple also contains

- the prolongation `P` from the periodic points to all the points (see
  [`periodic_restriction`](@ref)), its compartment blocks `P_cmpts`, and the index of one
  point for each periodic point (`masters`);
- the matrices `Cx` (``\\int \\varphi_i (D \\nabla \\varphi_j)_d``, one per axis `d`, see
  [`assemble_derivative_matrix`](@ref)).

Without periodicity, these fields are `nothing`.
"""
function assemble_matrices(model::Model{T,dim}; periodic = false) where {T,dim}
    (; mesh, D, T2) = model
    ncompartment = length(mesh.points)
    periodic isa Bool && (periodic = fill(periodic, dim))
    length(periodic) == dim || error("Give `periodic` as a `Bool` or one per axis")
    dims = findall(periodic)

    M_cmpts = SparseMatrixCSC{T,Int}[]
    S_cmpts = SparseMatrixCSC{T,Int}[]
    Mx_cmpts = [SparseMatrixCSC{T,Int}[] for _ = 1:dim]
    G = Matrix{T}[]
    volumes = zeros(T, ncompartment)
    for icmpt = 1:ncompartment
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        pattern = FEPattern(elements, size(points, 2))
        volumes[icmpt] = sum(first(element_geometry(points, elements)))

        push!(M_cmpts, assemble_mass_matrix(elements, points; pattern))
        push!(S_cmpts, assemble_stiffness_matrix(elements, points, D[icmpt]; pattern))
        for d = 1:dim
            Mx = assemble_mass_matrix(elements, points; nodeweights = points[d, :], pattern)
            push!(Mx_cmpts[d], Mx)
        end

        # Surface integrals of the basis functions times the outward normal
        facets = reduce(hcat, mesh.facets[icmpt, :])
        areas, _, normals = facet_geometry(points, facets, elements)
        G_i = zeros(T, size(points, 2), dim)
        for f in axes(facets, 2), a in axes(facets, 1), d = 1:dim
            G_i[facets[a, f], d] += areas[f] / size(facets, 1) * normals[d, f]
        end
        push!(G, G_i)
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T2)...)
    Mx = [blockdiag(Mx_cmpts[d]...) for d = 1:dim]
    Q = assemble_flux_matrix(model)

    if isempty(dims)
        P = P_cmpts = masters = Cx = nothing
    else
        # The outer boundaries are the faces of the box, which are not boundaries with
        # periodicity
        outer = [j for j in axes(mesh.facets, 2) if count(!isempty, mesh.facets[:, j]) == 1]
        all(iszero, model.kappa[outer]) || error(
            "Outer boundaries with relaxation (nonzero `kappa`) are not available with periodic boundary conditions",
        )

        # The box is the bounding box of all the compartments
        allpoints = reduce(hcat, mesh.points)
        box = vec(minimum(allpoints; dims = 2)), vec(maximum(allpoints; dims = 2))
        atol = 1e-8 * norm(box[2] - box[1])
        restrictions = map(1:ncompartment) do i
            facets = reduce(hcat, mesh.facets[i, :])
            periodic_restriction(mesh.points[i], facets, box, dims; atol)
        end
        P_cmpts = first.(restrictions)
        P = blockdiag(P_cmpts...)
        offsets = cumsum([0; size.(mesh.points, 2)])
        masters =
            reduce(vcat, [m .+ offsets[i] for (i, m) in enumerate(last.(restrictions))])
        Cx = map(1:dim) do d
            blockdiag(
                (
                    assemble_derivative_matrix(mesh.elements[i], mesh.points[i], d, D[i]) for i = 1:ncompartment
                )...,
            )
        end
    end

    (; M, S, R, Mx, Q, M_cmpts, S_cmpts, Mx_cmpts, G, volumes, P, P_cmpts, masters, Cx)
end

"""
    assemble_flux_matrix(model; symmetrical = false)

Assemble the global flux matrix `Q` of the permeable interfaces between compartments and of
the relaxing outer boundaries.

At an interface with permeability `kappa` between compartments with initial spin densities
`rho1` and `rho2`, the coefficients on each side are `2 rho2 / (rho1 + rho2) kappa` and
`2 rho1 / (rho1 + rho2) kappa`, so that the equilibrium magnetization has no net flux. With
`symmetrical = true`, `kappa` is used on both sides.
"""
function assemble_flux_matrix(model; symmetrical = false)
    (; mesh, rho, kappa) = model
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
            iszero(kappa[iboundary]) && continue
            f = facets[c, iboundary]
            g = f .+ offsets[c]
            facet_mass_triplets!(I, J, V, points[c], f, g, g, kappa[iboundary])
        elseif length(cmpts) == 2
            # Interface: the facets in the local numbering of both compartments
            c1, c2 = cmpts
            f1 = facets[c1, iboundary]
            f2 = local_indices(point_map[c2], point_map[c1][f1])
            g1 = f1 .+ offsets[c1]
            g2 = f2 .+ offsets[c2]

            if symmetrical
                kappa1 = kappa2 = kappa[iboundary]
            else
                rho1, rho2 = rho[c1], rho[c2]
                kappa1 = real(2rho2 / (rho1 + rho2)) * kappa[iboundary]
                kappa2 = real(2rho1 / (rho1 + rho2)) * kappa[iboundary]
            end

            # With the mass matrix M of the interface: Q11 += kappa1 M, Q12 -= kappa2 M,
            # Q21 -= kappa1 M, Q22 += kappa2 M
            facet_mass_triplets!(I, J, V, points[c1], f1, g1, g1, kappa1)
            facet_mass_triplets!(I, J, V, points[c1], f1, g1, g2, -kappa2)
            facet_mass_triplets!(I, J, V, points[c1], f1, g2, g1, -kappa1)
            facet_mass_triplets!(I, J, V, points[c1], f1, g2, g2, kappa2)
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
        i <= length(point_map) && point_map[i] == g ||
            error("Point $g is not in the compartment")
        i
    end
end
