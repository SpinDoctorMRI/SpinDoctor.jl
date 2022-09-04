"""
    create_mesh(surfaces, refinement)

Create mesh from surface geometry. A refinement is applied if `refinement < Inf`.
For 3D, TetGen is called. For 2D, Triangle is called.
"""
function create_mesh(surfaces, refinement)
    (; points, facets, facetmarkers, regions) = surfaces

    dim = size(points, 1)
    nregion = size(regions, 2)

    if dim == 2
        regionlist = [regions; (1:nregion)'; fill(refinement, 1, nregion)]

        triin = TriangulateIO()
        triin.pointlist = points
        triin.segmentlist = facets
        triin.segmentmarkerlist = facetmarkers
        triin.regionlist = regionlist
        triout = first(triangulate("pcaAQ", triin))

        points = triout.pointlist
        facets = Int.(triout.segmentlist)
        facetmarkers = Int.(triout.segmentmarkerlist)
        elements = Int.(triout.trianglelist)
        elementmarkers = round.(Int, triout.triangleattributelist[:])
    elseif dim == 3
        input = RawTetGenIO{Cdouble}()
        input.pointlist = points

        facetlist!(input, facets)

        input.facetmarkerlist = facetmarkers

        input.regionlist = [
            regions
            collect(1:nregion)'
            fill(0.1, 1, nregion)
        ]

        if isinf(refinement)
            tetgen = tetrahedralize(input, "pqA")
        else
            tetgen = tetrahedralize(input, "pqAa$(refinement)")
        end

        points = tetgen.pointlist
        facets = Int.(tetgen.trifacelist)
        facetmarkers = Int.(tetgen.trifacemarkerlist)
        elements = Int.(tetgen.tetrahedronlist)
        elementmarkers = Int.(tetgen.tetrahedronattributelist[:])
    end

    (; points, facets, facetmarkers, elements, elementmarkers)
end
