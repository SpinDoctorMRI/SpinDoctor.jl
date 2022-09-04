update!(p::Printer, problem, gradient, ξ, t) = p.verbosity ≥ 2 && @info "t = $t"

function update!(writer::VTKWriter, problem, gradient, ξ, t)
    if writer.n % writer.nupdate == 0
        (; pvd, dir, filename, ifile) = writer
        (; model) = problem

        ncompartment = length(model.mesh.points)
        npoint_cmpts = size.(model.mesh.points, 2)
        inds_cmpts = [0; cumsum(npoint_cmpts[1:end])]

        vtmfile = vtk_multiblock("$dir/$(filename)_$ifile")

        for icmpt = 1:ncompartment
            points = model.mesh.points[icmpt]
            elements = model.mesh.elements[icmpt]

            cells = [
                MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for i = 1:size(elements, 2)
            ]
            ξ_cmpt = ξ[1+inds_cmpts[icmpt]:inds_cmpts[icmpt+1]]

            vtkfile = vtk_grid(vtmfile, points, cells)
            vtkfile["Magnetization (real part)", VTKPointData()] = real(ξ_cmpt)
            vtkfile["Magnetization (imaginary part)", VTKPointData()] = imag(ξ_cmpt)
            vtkfile["Gradient", VTKFieldData()] = gradient(t)
        end

        pvd[t] = vtmfile
        writer.ifile += 1
    end
    writer.n += 1
end

function update!(p::Plotter, problem, gradient, ξ, t)
    if p.n % p.nupdate == 0
        # ξ updates first: important!
        p.ξ[] = ξ
        p.t[] = t

        sleep(0.01)
    end
    p.n += 1
end
