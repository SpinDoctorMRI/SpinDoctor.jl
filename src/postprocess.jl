"""
    compute_signal(M, ξ)

Compute signal from magnetization `ξ`, using the mass matrix `M` for integration.

Given a mesh `mesh` and a vector of compartment mass matrices `M_cmpts`, a vector of
compartment signals may be obtained by `compute_signal.(M_cmpts, split_field(mesh, ξ))`.
"""
compute_signal(M, ξ) = sum(M * ξ)

"""
    fit_adc(bvalues, signals)

Fit the apparent diffusion coefficient (ADC) using a polynomial logfit of the normalized signals
`signals` against the b-values `bvalues`.
"""
function fit_adc(bvalues, signals)
    # Fit log of signal
    p = fit(bvalues, log.(abs.(signals)))

    # ADC is the negative linear term
    -p[1]
end

"""
    fit_tensors(directions, adcs)

Fit effective diffusion tensors to directionalized ADCs.
The six components of the symmetric diffusion tensors are fitted by least
squares to the gradient directions and resulting ADCs.
"""
function fit_tensors(directions, adcs)
    T = eltype(directions)
    dim = size(directions, 1)

    adc = reduce(hcat, adcs)
    ncompartment, ndirection = size(adc)

    # Get vectors of gradient interactions
    d = directions
    if dim == 2
        G = zeros(ndirection, 3)
        G[:, 1] = d[1, :] .^ 2
        G[:, 2] = d[2, :] .^ 2
        G[:, 3] = 2 .* d[1, :] .* d[2, :]
    elseif dim == 3
        G = zeros(ndirection, 6)
        G[:, 1] = d[1, :] .^ 2
        G[:, 2] = d[2, :] .^ 2
        G[:, 3] = d[3, :] .^ 2
        G[:, 4] = 2 .* d[1, :] .* d[2, :]
        G[:, 5] = 2 .* d[1, :] .* d[3, :]
        G[:, 6] = 2 .* d[2, :] .* d[3, :]
    end

    tensors = [zeros(T, 0, 0) for _ = 1:ncompartment]
    for icmpt = 1:ncompartment
        # Deduce vector of effective diffusion tensor components using a
        # least squares fit of the directions to the computed ADCs
        Dvec = G \ adc[icmpt, :]

        # Tensorize component vectors (symmetric)
        if dim == 2
            tensors[icmpt] = [
                Dvec[1, :] Dvec[3, :]
                Dvec[3, :] Dvec[2, :]
            ]
        elseif dim == 3
            tensors[icmpt] = [
                Dvec[1, :] Dvec[4, :] Dvec[5, :]
                Dvec[4, :] Dvec[2, :] Dvec[6, :]
                Dvec[5, :] Dvec[6, :] Dvec[3, :]
            ]
        end
    end

    tensors
end

"""
    savefield(mesh, ξ, filename; fieldname = "Magnetization")

Save field `ξ` to a VTK file. It may then be visualized in ParaView.
"""
function savefield(
    mesh::FEMesh{T,dim},
    ξ,
    filename;
    fieldname = "Magnetization",
) where {T,dim}
    if dim == 2
        cell_type = VTKCellTypes.VTK_TRIANGLE
    elseif dim == 3
        cell_type = VTKCellTypes.VTK_TETRA
    end
    ξ_cmpts = split_field(mesh, ξ)
    isdir(dirname(filename)) || mkpath(dirname(filename))
    vtmfile = vtk_multiblock(filename)
    for icmpt = 1:length(mesh.points)
        ξᵢ = ξ_cmpts[icmpt]
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        cells = [MeshCell(cell_type, elements[:, i]) for i = 1:size(elements, 2)]
        vtkfile = vtk_grid(vtmfile, points, cells)
        vtkfile[fieldname*" (real part)"] = real(ξᵢ)
        vtkfile[fieldname*" (imaginary part)"] = imag(ξᵢ)
    end
    vtk_save(vtmfile)
end
