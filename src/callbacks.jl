# Callbacks for the BTPDE time stepping

"""
    AbstractCallback

Callback to be performed during the time stepping of the BTPDE. The callback is initialized
before time stepping with `initialize!(callback, problem, gradient, u, t)`, updated after
every step with `update!(callback, problem, gradient, u, t)`, and finalized after time
stepping with `finalize!(callback)`, where `problem = (; model, matrices)` and `u` is the
magnetization at time `t`.
"""
abstract type AbstractCallback end

"""
    Printer(; nupdate = 1, verbosity = 2)

Print time stepping information to the console every `nupdate` time step (if
`verbosity >= 2`).
"""
Base.@kwdef mutable struct Printer <: AbstractCallback
    nupdate::Int = 1
    verbosity::Int = 2
    n::Int = 1
end

initialize!(p::Printer, problem, gradient, u, t) =
    p.verbosity >= 1 && @info "Solving the BTPDE" gradient

function update!(p::Printer, problem, gradient, u, t)
    p.verbosity >= 2 && p.n % p.nupdate == 0 && @info "t = $t"
    p.n += 1
end

finalize!(::Printer) = nothing

"""
    VTKWriter(; nupdate = 1, dir = "output/time_stepping", filename = "solution")

Write the magnetization to VTK files every `nupdate` time step. The files are collected in
the ParaView data collection file `dir/filename.pvd`, which may be opened in ParaView. The
compartments are separate blocks.
"""
Base.@kwdef mutable struct VTKWriter <: AbstractCallback
    nupdate::Int = 1
    dir::String = "output/time_stepping"
    filename::String = "solution"
    n::Int = 1
    ifile::Int = 1
    pvd::Any = nothing
end

function initialize!(writer::VTKWriter, problem, gradient, u, t)
    mkpath(writer.dir)
    writer.n = 1
    writer.ifile = 1
    writer.pvd = paraview_collection(joinpath(writer.dir, writer.filename))
    update!(writer, problem, gradient, u, t)
end

function update!(writer::VTKWriter, problem, gradient, u, t)
    if writer.n % writer.nupdate == 0
        (; dir, filename, ifile, pvd) = writer
        vtm = vtk_multiblock(joinpath(dir, "$(filename)_$ifile"))
        for vtk in add_vtk_field!(vtm, problem.model.mesh, u, "Magnetization")
            vtk["Gradient", VTKFieldData()] = gradient(t)
        end
        pvd[t] = vtm
        writer.ifile += 1
    end
    writer.n += 1
end

finalize!(writer::VTKWriter) = vtk_save(writer.pvd)

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
    u::Any = nothing
    fig::Any = nothing
end

initialize!(p::Plotter, problem, gradient, u, t) =
    initialize_plotter!(p, problem, gradient, u, t)
update!(p::Plotter, problem, gradient, u, t) = update_plotter!(p, problem, gradient, u, t)
finalize!(::Plotter) = nothing
