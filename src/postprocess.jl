# Signals, ADC and tensor fits, and output of fields

"""
    compute_signal(M, u)

Compute the signal from the magnetization `u`, using the mass matrix `M` for integration.

Given a mesh `mesh` and the compartment mass matrices `M_cmpts`, the compartment signals
are `compute_signal.(M_cmpts, split_field(mesh, u))`.
"""
compute_signal(M, u) = sum(M * u)

"""
    fit_adc(bvalues, signals; degree = length(bvalues) - 1)

Compute the apparent diffusion coefficient (ADC) as the initial slope of the logarithm of
the `signals` with respect to the `bvalues`, by fitting a polynomial of degree `degree`
(by default interpolating all the points; use `degree = 1` for a log-linear least squares
fit).
"""
function fit_adc(bvalues, signals; degree = length(bvalues) - 1)
    p = fit(bvalues, log.(abs.(signals)), degree)
    -p[1]
end

"""
    fit_tensors(directions, adcs)

Fit effective diffusion tensors to the ADCs `adcs[i]` (one value per compartment) in the
`directions[:, i]`. The components of the symmetric tensors are fitted by least squares.
"""
function fit_tensors(directions, adcs)
    dim, ndirection = size(directions)
    adc = reduce(hcat, adcs)
    d = directions

    # The ADC in the direction d is d' * D * d, linear in the tensor components
    if dim == 2
        pairs = [(1, 1), (2, 2), (1, 2)]
    else
        pairs = [(1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)]
    end
    G = [(i == j ? 1 : 2) * d[i, k] * d[j, k] for k = 1:ndirection, (i, j) in pairs]

    map(eachrow(adc)) do a
        c = G \ a
        D = zeros(eltype(c), dim, dim)
        for ((i, j), cij) in zip(pairs, c)
            D[i, j] = D[j, i] = cij
        end
        D
    end
end

"""
    save_field(filename, mesh, u; name = "Magnetization")

Save the field `u` (real and imaginary parts) to the VTK multiblock file `filename.vtm`,
with one block per compartment. It may then be visualized in ParaView.
"""
function save_field(filename, mesh, u; name = "Magnetization")
    mkpath(dirname(abspath(filename)))
    vtm = vtk_multiblock(filename)
    add_vtk_field!(vtm, mesh, u, name)
    vtk_save(vtm)
end

# Add one VTK grid per compartment to the multiblock file `vtm`, with the real and
# imaginary parts of the field `u`, and return the grids
function add_vtk_field!(vtm, mesh::FEMesh{T,dim}, u, name) where {T,dim}
    celltype = dim == 2 ? VTKCellTypes.VTK_TRIANGLE : VTKCellTypes.VTK_TETRA
    map(mesh.points, mesh.elements, split_field(mesh, u)) do points, elements, u_i
        cells = [MeshCell(celltype, e) for e in eachcol(elements)]
        vtk = vtk_grid(vtm, points, cells)
        vtk[name*" (real part)", VTKPointData()] = real(u_i)
        vtk[name*" (imaginary part)", VTKPointData()] = imag(u_i)
        vtk
    end
end
