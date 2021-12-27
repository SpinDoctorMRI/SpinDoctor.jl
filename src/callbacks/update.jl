function update!(writer::VTKWriter, simulation, gradient, ξ, t)
    if writer.n % writer.nupdate == 0 
        (; pvd, dir, filename, ifile) = writer
        (; model) = simulation

        ncompartment = length(model.mesh.points)
        npoint_cmpts = size.(model.mesh.points, 2)
        inds_cmpts = [0; cumsum(npoint_cmpts[1:end])]

        vtmfile = vtk_multiblock("$dir/$(filename)_$ifile")
        
        for icmpt = 1:ncompartment
            points = model.mesh.points[icmpt]
            elements = model.mesh.elements[icmpt]

            cells = [
                MeshCell(VTKCellTypes.VTK_TETRA, elements[:, i]) for
                i = 1:size(elements, 2)
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

function update!(p::Plotter, simulation, gradient, ξ, t)
    if p.n % p.nupdate == 0 
        femesh = simulation.model.mesh
        M = simulation.matrices.M
        ncompartment, nboundary = size(femesh.facets)
        ξ_cmpts = split_field(femesh, ξ)

        push!(p.t[], t)
        p.f[] = push!(p.f[], gradient.profile(t))
        p.attenuation[] = push!(p.attenuation[], sum(abs, M * ξ) / p.S₀)
        p.ξ[] = ξ
        for icmpt = 1:ncompartment, iboundary = 1:nboundary
            p.magnitude[icmpt, iboundary][] = abs.(ξ_cmpts[icmpt])
            p.phase[icmpt, iboundary][] = angle.(ξ_cmpts[icmpt])
        end

        sleep(0.05)
    end
    p.n += 1
end
