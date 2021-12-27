function savefield(mesh, field, filename::String, fieldname = "Magnetization")

    # Make sure that directory exists
    isdir(dirname(filename)) || mkpath(dirname(filename))

    vtmfile = vtk_multiblock(filename)

    for icmpt = 1:length(mesh.points)
        f = field[icmpt]
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]

        cells =
            [MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)]

        vtkfile = vtk_grid(vtmfile, points, cells)
        if size(f, 2) == 1
            vtkfile[fieldname*" (real part)"] = real(f)
            vtkfile[fieldname*" (imaginary part)"] = imag(f)
        else
            for ifield = 1:min(size(f, 2), 99)
                fieldfieldname = "$fieldname $(@sprintf("%2d",ifield))"
                vtkfile[fieldfieldname*" (real part)"] = real(f[:, ifield])
                vtkfile[fieldfieldname*" (imaginary part)"] = imag(f[:, ifield])
            end
        end
    end

    vtk_save(vtmfile)
end
