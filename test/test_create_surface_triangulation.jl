using SpinDoctor, Test

# Test create cylinders
cellsetup = CellSetup(
    name = "test",
    shape = "cylinder",
    ncell = 10,
    rmin = 1.,
    rmax = 5.,
    dmin = 0.2,
    dmax = 0.3,
    height = 3.,
    include_nucleus = true,
    nucleus_radiusratio = 0.7,
    include_ecs = true,
    ecs_shape = "convexhull",
    ecs_gap = 0.3,
)
domainsetup = DomainSetup(
    diffusivity_in = 0.002,
    diffusivity_out = 0.002,
    diffusivity_ecs = 0.002,
    relaxation_in = Inf,
    relaxation_out = Inf,
    relaxation_ecs = Inf,
    initial_density_in = 1.,
    initial_density_out = 1.,
    initial_density_ecs = 1.,
    permeability_in_out = 1e-4,
    permeability_out_ecs = 1e-4,
)

cells = create_cells(cellsetup)
domain = prepare_pde(cellsetup, domainsetup)
tri = create_surface_triangulation(cellsetup, cells, domain)
@test size(tri.points)[2] == 3
@test size(tri.facets)[2] == 3
@test length(tri.facetmarkers) == size(tri.facets)[1]
@test size(tri.regions) == (domain.ncompartment, 3)
