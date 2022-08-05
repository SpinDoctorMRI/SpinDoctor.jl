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

function update!(p::Plotter{T,dim}, problem, gradient, ξ, t) where {T,dim}
    if p.n % p.nupdate == 0
        femesh = problem.model.mesh
        M = problem.matrices.M
        ncompartment, nboundary = size(femesh.facets)
        ξ_cmpts = split_field(femesh, ξ)
        push!(p.t[], t)
        if gradient isa ScalarGradient
            p.f[] = push!(p.f[], gradient.profile(t))
        else
            grad = Vec{dim,Float32}(gradient(t))
            p.g⃗[] = [grad]
            p.g⃗_hist[] = push!(p.g⃗_hist[], grad)
        end
        p.attenuation[] = push!(p.attenuation[], sum(abs, M * ξ) / p.S₀)
        p.ξ[] = ξ
        if dim == 2
            for icmpt = 1:ncompartment
                p.magnitude[icmpt][] = abs.(ξ_cmpts[icmpt])
                p.phase[icmpt][] = angle.(ξ_cmpts[icmpt])
            end
        elseif dim == 3
            for icmpt = 1:ncompartment, iboundary = 1:nboundary
                p.magnitude[icmpt, iboundary][] = abs.(ξ_cmpts[icmpt])
                p.phase[icmpt, iboundary][] = angle.(ξ_cmpts[icmpt])
            end
        end

        sleep(0.01)
    end
    p.n += 1
end
