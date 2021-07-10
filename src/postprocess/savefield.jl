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


function savefield_time(mesh, time, field, filename::String, fieldname = "Magnetization")

    # Make sure that directory exists
    isdir(dirname(filename)) || mkpath(dirname(filename))

    ncompartment = length(mesh.points)

    pvd = paraview_collection(filename)

    for (it, t) ∈ enumerate(time)
        vtmfile = vtk_multiblock(filename * "_$it")

        for icmpt = 1:ncompartment
            points = mesh.points[icmpt]
            elements = mesh.elements[icmpt]

            cells = [
                MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for
                i = 1:size(elements, 2)
            ]

            # vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile[fieldname*" (real part)"] = real(field[icmpt][:, it])
            vtkfile[fieldname*" (imaginary part)"] = imag(field[icmpt][:, it])
        end

        pvd[t] = vtmfile
    end

    vtk_save(pvd)
end


function save_btpde_results(mesh, experiment, btpde, filename)

    ncompartment = length(mesh.points)
    g = experiment.gradient.directions[:, 1]
    f = experiment.gradient.sequences[1]

    time = btpde.time[1, 1, 1]
    magnetization = btpde.magnetization[:, 1, 1, 1]

    pvd = paraview_collection(filename)

    for (it, t) ∈ enumerate(time)
        vtmfile = vtk_multiblock(filename * "_$it")

        for icmpt = 1:ncompartment
            points = mesh.points[icmpt]
            elements = mesh.elements[icmpt]

            cells = [
                MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for
                i = 1:size(elements, 2)
            ]

            # vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile["Magnetization (real part)", VTKPointData()] =
                real(magnetization[icmpt][:, it])
            vtkfile["Magnetization (imaginary part)", VTKPointData()] =
                imag(magnetization[icmpt][:, it])
            vtkfile["Gradient", VTKFieldData()] = f(t) * g
            vtkfile["Sequence", VTKFieldData()] = f(t)
        end

        pvd[t] = vtmfile
    end

    vtk_save(pvd)
end
