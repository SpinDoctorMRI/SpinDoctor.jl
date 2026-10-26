finalize!(::Printer) = nothing
finalize!(::AbstractPlotter) = nothing
finalize!(writer::VTKWriter) = vtk_save(writer.pvd)
