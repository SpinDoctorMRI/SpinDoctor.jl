"""
    savefield(mesh, ξ, filename; fieldname = "Magnetization")

Save field `ξ` to a VTK file. It may then be visualized in ParaView.
"""
function savefield(mesh::FEMesh{T,dim}, ξ, filename; fieldname = "Magnetization") where {T,dim}
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
        cells =
            [MeshCell(cell_type, elements[:, i]) for i = 1:size(elements, 2)]
        vtkfile = vtk_grid(vtmfile, points, cells)
        vtkfile[fieldname*" (real part)"] = real(ξᵢ)
        vtkfile[fieldname*" (imaginary part)"] = imag(ξᵢ)
    end
    vtk_save(vtmfile)
end
