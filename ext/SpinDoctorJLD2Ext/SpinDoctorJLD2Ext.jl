"""
Saving and loading geometries in JLD2 files.
"""
module SpinDoctorJLD2Ext

using JLD2
using SpinDoctor
using SpinDoctor: FEMesh
import SpinDoctor: save_geometry, load_geometry

function save_geometry(file, setup, mesh, surfaces, cells)
    mkpath(dirname(abspath(file)))
    jldsave(
        file;
        format = "SpinDoctor geometry 1",
        setup = repr(setup),
        point_map = mesh.point_map,
        points = mesh.points,
        facets = mesh.facets,
        elements = mesh.elements,
        surfaces,
        cells,
    )
end

function load_geometry(file)
    jldopen(file, "r") do f
        f["format"] == "SpinDoctor geometry 1" || error("Unknown geometry file format")
        mesh = FEMesh(;
            point_map = f["point_map"],
            points = f["points"],
            facets = f["facets"],
            elements = f["elements"],
        )
        (; setup = f["setup"], mesh, surfaces = f["surfaces"], cells = f["cells"])
    end
end

end
