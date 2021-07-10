@testset "create_surfaces_cylinder.jl" begin
    # Test create cylinders
    setup = CylinderSetup(
        name = "test",
        ncell = 10,
        rmin = 1.0,
        rmax = 5.0,
        dmin = 0.2,
        dmax = 0.3,
        bend = 0.0,
        twist = 0.0,
        height = 3.0,
        include_in = true,
        in_ratio = 0.6,
        ecs_shape = "convex_hull",
        ecs_ratio = 0.5,
        refinement = 0.5,
        D_in = 0.002 * I(3),
        D_out = 0.002 * I(3),
        D_ecs = 0.002 * I(3),
        T₂_in = Inf,
        T₂_out = Inf,
        T₂_ecs = Inf,
        ρ_in = 1.0,
        ρ_out = 1.0,
        ρ_ecs = 1.0,
        κ_in_out = 1e-5,
        κ_out_ecs = 1e-5,
        κ_in = 0.0,
        κ_out = 0.0,
        κ_ecs = 0.0,
    )

    cells = create_cells(setup)
    surfaces = create_surfaces_cylinder(cells, setup)

    ncompartment = 2 * setup.ncell + 1

    @test size(surfaces.points)[1] == 3
    @test size(surfaces.facets)[1] == 3
    @test length(surfaces.facetmarkers) == size(surfaces.facets)[2]
    @test size(surfaces.regions) == (3, ncompartment)
end
