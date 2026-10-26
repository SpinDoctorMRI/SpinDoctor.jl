finalize!(::Printer) = nothing
finalize!(::Plotter) = nothing
finalize!(writer::VTKWriter) = vtk_save(writer.pvd)
