function savefield(mesh, field, filename::String, fieldname="Magnetization")

    vtmfile = vtk_multiblock(filename)

    for icmpt = 1:mesh.ncmpt

        f = field[icmpt]
        points = mesh.points[icmpt]
        elements = mesh.elements[icmpt]

        cells = [MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)]

        vtkfile = vtk_grid(vtmfile, points, cells)
        if size(f, 2) == 1
            vtkfile[fieldname] = real(f)
        else
            for ifield = 1:min(size(f, 2), 99)
                vtkfile["$(fieldname) $(@sprintf("%2d",ifield))"] = real(f[:, ifield])
            end
        end

    end

    outfiles = vtk_save(vtmfile)
end


function savefield_time(mesh, time, field, filename::String, fieldname="Magnetization")

    ncmpt = mesh.ncmpt

    pvd = paraview_collection(filename)

    for (it, t) ∈ enumerate(time)

        vtmfile = vtk_multiblock(filename * "_$it")


        for icmpt = 1:ncmpt
            points = mesh.points[icmpt]
            elements = mesh.elements[icmpt]

            cells = [MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)]

            # vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile[fieldname] = real(field[icmpt][:, it])
        end

        pvd[t] = vtmfile

    end

    outfiles = vtk_save(pvd)

end


function save_btpde_results(mesh, btpde, experiment, directions, filename)

    ncmpt = mesh.ncmpt
    g = directions[:, 1]
    f = experiment.sequences[1]

    time = btpde.time[1, 1, 1]
    magnetization = btpde.magnetization[:, 1, 1, 1]
    signal_allcmpts = btpde.signal_allcmpts[1, 1, 1]

    pvd = paraview_collection(filename)

    for (it, t) ∈ enumerate(time)

        vtmfile = vtk_multiblock(filename * "_$it")

        for icmpt = 1:ncmpt
            points = mesh.points[icmpt]
            elements = mesh.elements[icmpt]

            cells = [MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)]

            # vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile["Magnetization", VTKPointData()] = real(magnetization[icmpt][:, it])
            vtkfile["Gradient", VTKFieldData()] = f(t) * g
            vtkfile["Sequence", VTKFieldData()] = f(t)
            vtkfile["Signal", VTKFieldData()] = real(signal_allcmpts[it])
        end

        pvd[t] = vtmfile

    end

    outfiles = vtk_save(pvd)

end
