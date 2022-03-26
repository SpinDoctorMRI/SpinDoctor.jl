"""
    call_tetgen(surfaces, refinement)

Call Tetgen on surface geometry. A refinement is applied if `refinement < Inf`.
"""
function call_tetgen(surfaces, refinement, tetgen_option)
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

    if isempty(tetgen_option)
        if isinf(refinement)
            tetgen = tetrahedralize(input, "pqACR")
        else
            tetgen = tetrahedralize(input, "pqACRa$(refinement)")
        end
    else
        if isinf(refinement)
            tetgen = tetrahedralize(input, "$(tetgen_option)")
        else
            tetgen = tetrahedralize(input, "$(tetgen_option)a$(refinement)")
        end
    end

    points = tetgen.pointlist
    facets = Int.(tetgen.trifacelist)
    facetmarkers = Int.(tetgen.trifacemarkerlist)
    elements = Int.(tetgen.tetrahedronlist)
    elementmarkers = Int.(tetgen.tetrahedronattributelist[:])

    (; points, facets, facetmarkers, elements, elementmarkers)
end
