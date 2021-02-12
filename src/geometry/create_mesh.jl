function create_mesh(cellsetup::CellSetup, domain, tri)
    @unpack ncell, refinement = cellsetup
    @unpack ncmpt = domain
    @unpack points, facets, facetmarkers, regions = tri

    input = TetGen.RawTetGenIO{Cdouble}()
    input.pointlist = points'

    TetGen.facetlist!(input, facets')

    input.facetmarkerlist = facetmarkers

    input.regionlist = [regions 1:ncmpt repeat([0.1], ncmpt)]'

    isnothing(refinement) ? tetrahedralize(input, "pqA") : tetrahedralize(input, "pqAa$(refinement)")
end
