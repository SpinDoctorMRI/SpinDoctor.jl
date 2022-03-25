@testset "create_cells.jl" begin
    # Test create cylinders
    cylsetup = CylinderSetup(
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
        ecs_shape = :convex_hull,
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

    cells = create_cells(cylsetup)
    @test cells isa NamedTuple
    @test size(cells.radii) == (1, cylsetup.ncell)
    @test size(cells.centers) == (2, cylsetup.ncell)
    @test all(cylsetup.rmin .≤ cells.radii .≤ cylsetup.rmax)


    # Test create spheres
    spheresetup = SphereSetup(
        name = "test",
        ncell = 4,
        rmin = 3.0,
        rmax = 5.0,
        dmin = 0.2,
        dmax = 0.3,
        include_in = true,
        in_ratio = 0.7,
        ecs_shape = :box,
        ecs_ratio = 0.3,
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
    cells = create_cells(spheresetup)
    @test cells isa NamedTuple
    @test size(cells.radii) == (1, spheresetup.ncell)
    @test size(cells.centers) == (3, spheresetup.ncell)
    @test all(spheresetup.rmin .≤ cells.radii .≤ spheresetup.rmax)
end
