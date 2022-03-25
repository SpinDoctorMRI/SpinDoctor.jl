@testset "create_cells.jl" begin
    # Test create cylinders
    cylsetup = CylinderSetup(
        name = "test",
        ncell = 10,
        r_range = (; rmin = 3.0, rmax =  5.0),
        d_range = (; dmin = 0.2, dmax =  0.3),
        deform_angle = (; bend = 0.0, twist =  0.0),
        height = 3.0,
        include_in = true,
        in_ratio = 0.6,
        ecs_shape = :convex_hull,
        ecs_ratio = 0.5,
        refinement = 0.5,
        D = (; in = 2e-3 * I(3), out = 2e-3 * I(3), ecs = 2e-3 * I(3)),
        T₂ = (; in = Inf, out = Inf, ecs = Inf),
        ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
        κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
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
        r_range = (; rmin = 3.0, rmax =  5.0),
        d_range = (; dmin = 0.2, dmax =  0.3),
        include_in = true,
        in_ratio = 0.7,
        ecs_shape = :box,
        ecs_ratio = 0.3,
        D = (; in = 2e-3 * I(3), out = 2e-3 * I(3), ecs = 2e-3 * I(3)),
        T₂ = (; in = Inf, out = Inf, ecs = Inf),
        ρ = (; in = 1.0, out = 1.0, ecs = 1.0),
        κ = (; in_out = 1e-4, out_ecs = 1e-4, in = 0.0, out = 0.0, ecs = 0.0),
    )
    cells = create_cells(spheresetup)
    @test cells isa NamedTuple
    @test size(cells.radii) == (1, spheresetup.ncell)
    @test size(cells.centers) == (3, spheresetup.ncell)
    @test all(spheresetup.rmin .≤ cells.radii .≤ spheresetup.rmax)
end
