"""
    call_tetgen(setup, surfaces)

Call Tetgen on surface geometry.
"""
function call_tetgen(surfaces, refinement = nothing)
    (; points, facets, facetmarkers, regions) = surfaces

    nregion = size(regions, 2)

    input = RawTetGenIO{Cdouble}()
    input.pointlist = points

    facetlist!(input, facets)

    input.facetmarkerlist = facetmarkers

    input.regionlist = [
        regions
        collect(1:nregion)'
        fill(0.1, 1, nregion)
    ]

    if isnothing(refinement)
        tetgen = tetrahedralize(input, "pqA")
    else
        tetgen = tetrahedralize(input, "pqAa$(refinement)")
    end

    points = tetgen.pointlist
    facets = Int.(tetgen.trifacelist)
    facetmarkers = Int.(tetgen.trifacemarkerlist)
    elements = Int.(tetgen.tetrahedronlist)
    elementmarkers = Int.(tetgen.tetrahedronattributelist[:])

    (; points, facets, facetmarkers, elements, elementmarkers)
end
