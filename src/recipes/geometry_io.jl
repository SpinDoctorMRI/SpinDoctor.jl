# Saving and loading geometries in JLD2 files (implemented in the JLD2 extension, loaded
# with `using JLD2`)

"""
    save_geometry(file, setup, mesh, surfaces, cells)

Save a geometry created by [`create_geometry`](@ref) to the JLD2 file `file`, together with
a description of `setup`. Only plain arrays are stored, so the file does not depend on the
SpinDoctor version. Requires `using JLD2`.
"""
function save_geometry end

"""
    (; setup, mesh, surfaces, cells) = load_geometry(file)

Load a geometry saved with [`save_geometry`](@ref). The `setup` is returned as the
description string that was saved. Requires `using JLD2`.
"""
function load_geometry end
