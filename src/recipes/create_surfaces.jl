"""
    create_surfaces(setup, cells)

Create surface triangulation from a setup and cells.

Returns a named tuple containing:

  - `points` (size `dim × npoint`),
  - `facets` (edges in 2D) (size `dim × nfacet`),
  - `facetmarkers` (edge markers in 2D) (size `nfacet`),
  - `regions` (one point inside each region) (size `dim × nregion`)
"""
function create_surfaces end

include("create_surfaces_extrusion.jl")
include("create_surfaces_plates.jl")
include("create_surfaces_neuron.jl")
include("create_surfaces_disk_sphere.jl")
