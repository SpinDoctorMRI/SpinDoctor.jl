function create_mesh(cellsetup::CellSetup, domain, tri)
    @unpack ncell, refinement = cellsetup
    @unpack ncompartment = domain
    @unpack points, facets, facetmarkers, regions = tri

    input = TetGen.RawTetGenIO{Cdouble}()
    input.pointlist = points

    TetGen.facetlist!(input, facets)

    input.facetmarkerlist = facetmarkers

    input.regionlist = [
        regions
        collect(1:ncompartment)'
        repeat([0.1], 1, ncompartment)
    ]

    isnothing(refinement) ? tetrahedralize(input, "pqA") : tetrahedralize(input, "pqAa$(refinement)")
end
