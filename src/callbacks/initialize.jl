function initialize!(p::Printer, problem, gradient, ξ, t)
    @info "Solving for" typeof(problem) gradient
end

function initialize!(writer::VTKWriter, problem, gradient, ξ, t)
    (; dir, filename) = writer
    mkpath(dir)
    writer.pvd = paraview_collection(joinpath(dir, filename))
    update!(writer, problem, gradient, ξ, t)
end

initialize!(p::Plotter, problem, gradient, ξ, t) =
    initialize_plotter!(p, problem, gradient, ξ, t)
