"""
    create_surfaces(setup::ExtrusionSetup, cells)

Create surface triangulation of extruded setup.

The ground surface is triangulated first, before the walls are "extruded" and
the top surface is copied from the ground surface.
"""
function create_surfaces(setup::ExtrusionSetup, cells)
    (; name, groundsetup, height, refinement) = setup

    # using SpinDoctor: create_surfaces, create_mesh

    # Create ground surface triangulation using 2D disk setup
    groundfacets = create_surfaces(groundsetup, cells)
    groundmesh = create_mesh(groundfacets, groundsetup.refinement)

    # groundfacets.facetmarkers
    # plot_surfaces(groundfacets, 1:10)
    # plot_mesh(split_mesh(groundmesh), 1:4)

    groundregions = groundfacets.regions
    nregion = size(groundregions, 2)

    # 2D facets are 3D edges, 2D elements are 3D facets
    points = groundmesh.points
    edges = groundmesh.facets
    edgemarkers = groundmesh.facetmarkers
    facets = groundmesh.elements
    facetmarkers = groundmesh.elementmarkers

    # Copy ground points to the two 3D planes top and bottom
    npoint = size(points, 2)
    z = fill(height / 2, 1, npoint)
    points = [points points; -z z]

    # Create two side facet triangles from each ground edge. They form a
    # rectangle representing the extruded edge. The top layer points are
    # numbered `npoint` higher than the bottom layer
    le = edges[1:1, :]
    ri = edges[2:2, :]
    sidefacets = [
        le le
        ri (ri.+npoint)
        (ri.+npoint) (le.+npoint)
    ]

    # The boundaries are ordered as follows:
    #   Γ_1_2 (we skip Γ_2_1)
    #   Γ_1_3 (we skip Γ_3_1)
    #   ⋮
    #   Γ_nregion_nregion
    #   Γ_1
    #   ⋮
    #   Γ_nregion
    sidefacetmarkers = [edgemarkers; edgemarkers]

    # Bottom facets, top facets and side facets
    ninterface = (nregion - 1) * nregion ÷ 2
    facets = [sidefacets facets (facets .+ npoint)]
    facetmarkers =
        [sidefacetmarkers; facetmarkers .+ ninterface; facetmarkers .+ ninterface]

    # Put a point in the middle of each region (z = 0 is in the middle of
    # bottom and top layer)
    regions = [groundregions[1:2, :]; zeros(1, nregion)]

    (; points, facets, facetmarkers, regions)
end
