"""
    assemble_matrices(model)

Assemble finite element matrices.
"""
function assemble_matrices(model)
    @unpack mesh, D, T₂, ρ = model

    # Deduce sizes
    ncompartment = length(ρ)

    # Assemble finite element matrices compartment-wise
    M_cmpts = []
    S_cmpts = []
    Mx_cmpts = [[] for _ = 1:3]
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        volumes, = get_mesh_volumes(points, elements)

        # Assemble mass, stiffness and flux matrices
        push!(M_cmpts, assemble_mass_matrix(elements', volumes))
        push!(S_cmpts, assemble_stiffness_matrix(elements', points', D[icmpt]))

        # Assemble first order product moment matrices
        for dim = 1:3
            push!(Mx_cmpts[dim], assemble_mass_matrix(elements', volumes, points[dim, :]))
        end
    end

    # Assemble global finite element matrices
    M = blockdiag(M_cmpts...)
    S = blockdiag(S_cmpts...)
    R = blockdiag((M_cmpts ./ T₂)...)
    Mx = [blockdiag(Mx_cmpts[dim]...) for dim = 1:3]
    Q_blocks = assemble_flux_matrices(mesh.points, mesh.facets)
    Q = couple_flux_matrix(model, Q_blocks, false)

    (; M, S, R, Mx, Q, M_cmpts, S_cmpts, Mx_cmpts)
end
