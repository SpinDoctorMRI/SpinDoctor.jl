function savefield(mesh, ξ, filename::String; fieldname = "Magnetization")
    ξ_cmpts = split_field(mesh, ξ)
    isdir(dirname(filename)) || mkpath(dirname(filename))
    vtmfile = vtk_multiblock(filename)
    for icmpt = 1:length(mesh.points)
        ξᵢ = ξ_cmpts[icmpt]
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]
        cells =
            [MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)]
        vtkfile = vtk_grid(vtmfile, points, cells)
        vtkfile[fieldname * " (real part)"] = real(ξᵢ)
        vtkfile[fieldname * " (imaginary part)"] = imag(ξᵢ)
    end
    vtk_save(vtmfile)
end
