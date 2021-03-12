"""
    call_tetgen(setup, surfaces)

Call Tetgen on surface geometry.
"""
function call_tetgen(surfaces, refinement = nothing)
    @unpack points, facets, facetmarkers, regions = surfaces

    nregion = size(regions, 2)

    input = TetGen.RawTetGenIO{Cdouble}()
    input.pointlist = points

    TetGen.facetlist!(input, facets)

    input.facetmarkerlist = facetmarkers

    input.regionlist = [
        regions
        collect(1:nregion)'
        repeat([0.1], 1, nregion)
    ]

    tetgen = isnothing(refinement) ? tetrahedralize(input, "pqA") : tetrahedralize(input, "pqAa$(refinement)")

    points = tetgen.pointlist
    facets = Int.(tetgen.trifacelist)
    facetmarkers = Int.(tetgen.trifacemarkerlist)
    elements = Int.(tetgen.tetrahedronlist)
    elementmarkers = Int.(tetgen.tetrahedronattributelist[:])

    # Return named tuple
    (; points, facets, facetmarkers, elements, elementmarkers)
end
