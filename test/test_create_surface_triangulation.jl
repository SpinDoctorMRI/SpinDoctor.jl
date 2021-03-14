using SpinDoctor, Test

# Test create cylinders
cellsetup = CellSetup(
    name = "test",
    shape = "cylinder",
    ncell = 10,
    rmin = 1.0,
    rmax = 5.0,
    dmin = 0.2,
    dmax = 0.3,
    height = 3.0,
    include_in = true,
    in_ratio = 0.7,
    include_ecs = true,
    ecs_shape = "convexhull",
    ecs_ratio = 0.3,
)
domainsetup = DomainSetup(
    σ_in = 0.002,
    σ_out = 0.002,
    σ_ecs = 0.002,
    T₂_in = Inf,
    T₂_out = Inf,
    T₂_ecs = Inf,
    ρ_in = 1.0,
    ρ_out = 1.0,
    ρ_ecs = 1.0,
    κ_in_out = 1e-4,
    κ_out_ecs = 1e-4,
    κ_in = 0.0,
    κ_out = 0.0,
    κ_ecs = 0.0,
)

cells = create_cells(cellsetup)
domain = prepare_pde(cellsetup, domainsetup)
tri = create_surface_triangulation(cellsetup, cells, domain)
@test size(tri.points)[2] == 3
@test size(tri.facets)[2] == 3
@test length(tri.facetmarkers) == size(tri.facets)[1]
@test size(tri.regions) == (domain.ncompartment, 3)
