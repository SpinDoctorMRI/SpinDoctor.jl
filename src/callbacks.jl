"""
    AbstractCallback

Callback to be performed after every time step when solving the BTPDE. The callback is
initialized before time stepping using the `initialize!` function, updated after every step
using `update!`, and finalized after time stepping using `finalize!`.
"""
abstract type AbstractCallback end

"""
    Printer(; nupdate = 1)

Print time stepping information to the console.
"""
Base.@kwdef mutable struct Printer <: AbstractCallback
    verbosity::Int = 2
    nupdate::Int = 1
    n::Int = 1
end

"""
    VTKWriter(; nupdate = 1, dir = "output", filename = "solution")

Write magnetization field to a VTK file after every `nupdate` timestep. The files are stored in a
ParaView data collection file (PVD). The magnetization time series may be visualized in
ParaView by opening the file `"\$dir/\$filename.pvd"`. The compartments are labeled.
"""
Base.@kwdef mutable struct VTKWriter <: AbstractCallback
    nupdate::Int = 1
    n::Int = 1
    ifile::Int = 1
    dir::String = "output/time_stepping"
    filename::String = "solution"
    pvd::WriteVTK.CollectionFile = paraview_collection("")
end

"""
    Plotter(; nupdate = 1)

Plot the evolution of the BTPDE during time stepping. This requires loading a Makie
plotting backend (`using GLMakie` for an interactive window, or `using CairoMakie`). The
plot is updated every `nupdate` time step. The resulting figure contains a plot of the time
profile, total signal attenuation, and magnetization field (complex magnitude and phase
shift), and is available as `plotter.fig` after solving.
"""
Base.@kwdef mutable struct Plotter <: AbstractCallback
    nupdate::Int = 1
    n::Int = 1
    t::Any = nothing
    ξ::Any = nothing
    fig::Any = nothing
end

function initialize!(p::Printer, problem, gradient, ξ, t)
    @info "Solving the BTPDE" gradient
end

function initialize!(writer::VTKWriter, problem, gradient, ξ, t)
    (; dir, filename) = writer
    mkpath(dir)
    writer.pvd = paraview_collection(joinpath(dir, filename))
    update!(writer, problem, gradient, ξ, t)
end

initialize!(p::Plotter, problem, gradient, ξ, t) =
    initialize_plotter!(p, problem, gradient, ξ, t)

update!(p::Printer, problem, gradient, ξ, t) = p.verbosity ≥ 2 && @info "t = $t"

function update!(writer::VTKWriter, problem, gradient, ξ, t)
    if writer.n % writer.nupdate == 0
        (; pvd, dir, filename, ifile) = writer
        (; model) = problem

        ncompartment = length(model.mesh.points)
        npoint_cmpts = size.(model.mesh.points, 2)
        inds_cmpts = [0; cumsum(npoint_cmpts[1:end])]

        vtmfile = vtk_multiblock("$dir/$(filename)_$ifile")

        dim = size(model.mesh.points[1], 1)
        celltype = dim == 2 ? VTKCellTypes.VTK_TRIANGLE : VTKCellTypes.VTK_TETRA

        for icmpt = 1:ncompartment
            points = model.mesh.points[icmpt]
            elements = model.mesh.elements[icmpt]

            cells = [MeshCell(celltype, elements[:, i]) for i = 1:size(elements, 2)]
            ξ_cmpt = ξ[(1+inds_cmpts[icmpt]):inds_cmpts[icmpt+1]]

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

update!(p::Plotter, problem, gradient, ξ, t) = update_plotter!(p, problem, gradient, ξ, t)

finalize!(::Printer) = nothing
finalize!(::Plotter) = nothing
finalize!(writer::VTKWriter) = vtk_save(writer.pvd)
