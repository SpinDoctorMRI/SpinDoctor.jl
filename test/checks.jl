# Checks of computed values against exact identities and analytical solutions

γ = 2.67513e-4

@testset "Compartment boundaries" begin
    # Every boundary facet of a compartment is a face of one of its elements. The
    # saved meshes of bent and layered cylinders used to violate this (and could
    # not be assembled), since facets with all their points on an interface were
    # also assigned to the neighboring compartment.
    for (name, setup) ∈ [
        "cylinders_bent_hull" => CylinderSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = π / 4,
        ),
        "cylinders_layers_hull" => CylinderSetup(;
            ncell = 2,
            layersizes = [0.6, 1.0],
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = π / 4,
        ),
    ]
        savedir = joinpath(@__DIR__, "regression", "meshes", name)
        mesh, = create_geometry(setup; savedir, recreate = false)
        for i ∈ eachindex(mesh.points), f ∈ mesh.facets[i, :]
            isempty(f) || @test all(>(0), SpinDoctor.facet_elements(mesh.elements[i], f))
        end
        model = Model(; mesh, get_coeffs(setup)...)
        @test assemble_matrices(model) isa NamedTuple
    end
end

@testset "Outward normals" begin
    # Divergence theorem: ∫_∂Ω x_d n_d dS = |Ω| for each compartment Ω and axis d.
    # P1 elements represent x_d exactly, so the surface integrals G must reproduce
    # the compartment volumes, also for non-convex compartments such as the ECS.
    setups = [
        DiskSetup(; ncell = 3, rmin = 2.0, rmax = 4.0, ecs = BoxECS(; margin = 1.0)),
        SphereSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            nsidewall = 60,
            ecs = BoxECS(; margin = 1.0),
        ),
        CylinderSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = π / 4,
        ),
        CylinderSetup(;
            ncell = 2,
            layersizes = [0.6, 1.0],
            rmin = 2.0,
            rmax = 4.0,
            ecs = BoxECS(; margin = 1.0),
            height = 10.0,
        ),
    ]
    for setup ∈ setups
        Random.seed!(1)
        mesh, = create_geometry(setup)
        model = Model(; mesh, get_coeffs(setup)...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        for i ∈ eachindex(mesh.points), d = 1:dim
            @test mesh.points[i][d, :]' * matrices.G[i][:, d] ≈ matrices.volumes[i] rtol =
                1e-8
        end

        # With equal spin densities, permeable interfaces conserve magnetization
        # and the flux matrix is symmetric
        @test norm(matrices.Q * ones(size(matrices.Q, 1))) ≤ 1e-10 * norm(matrices.Q, 1)
        @test matrices.Q ≈ matrices.Q'
    end
end

@testset "Permeable plates vs analytical solution" begin
    # For a gradient across the plates, the 2D plates and the 3D slabs reduce to
    # the one-dimensional problem of the analytical matrix formalism
    profile = PGSE(2000.0, 6000.0)
    for setup ∈ [
        PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement = 0.5),
        SlabSetup(; depth = 20.0, widths = fill(5.0, 4), height = 5.0, refinement = 2.0),
    ]
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        dir = [1.0; zeros(dim - 1)]
        gradient = ScalarGradient(dir, profile, √(1000 / int_F²(profile)) / γ)
        ξ = solve(
            BTPDE(; model, matrices),
            gradient,
            IntervalConstantSolver(; timestep = 5.0),
        )
        volumes = get_cmpt_volumes(mesh)
        signal = compute_signal(matrices.M, ξ) / sum(volumes)
        acoeffs = analytical_coefficients(setup, coeffs)
        laplace = AnalyticalLaplace(; acoeffs..., eiglim = length2eig(0.3, 0.002), eigstep = 1e-8)
        amf = AnalyticalMatrixFormalism(;
            analytical_laplace = laplace,
            lap_mat = solve(laplace),
            volumes,
        )
        signal_analytical = solve(amf, gradient) / sum(volumes)
        @test real(signal) ≈ real(signal_analytical) rtol = 1e-2
    end
end

@testset "Utilities" begin
    # Directions on a hemisphere
    d = unitsphere(20; half = true, normal = [1.0, 1.0, 0.0])
    @test size(d) == (3, 20)
    @test all(≥(0), [1, 1, 0]' * d)
    @test all(≈(1), sum(abs2, d; dims = 1))

    # Reproducible cell configurations
    setup = SphereSetup(; ncell = 3)
    @test SpinDoctor.create_cells(setup; rng = Xoshiro(1)) ==
          SpinDoctor.create_cells(setup; rng = Xoshiro(1))

    # VTK output of 2D (triangle) and 3D (tetrahedron) meshes
    for setup ∈ [
        PlateSetup(; depth = 2.0, widths = [1.0, 1.0]),
        SlabSetup(; depth = 2.0, widths = [1.0, 1.0], height = 1.0),
    ]
        mesh, = create_geometry(setup)
        model = Model(; mesh, get_coeffs(setup)...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        gradient = ScalarGradient([1.0; zeros(dim - 1)], PGSE(10.0, 20.0), 0.01)
        mktempdir() do dir
            writer = VTKWriter(; dir = joinpath(dir, "nested", "output"), nupdate = 2)
            solver = IntervalConstantSolver(; timestep = 5.0)
            ξ = solve(BTPDE(; model, matrices), gradient, solver; callbacks = [writer])
            @test isfile(joinpath(dir, "nested", "output", "solution.pvd"))
            savefield(mesh, ξ, joinpath(dir, "field"))
            @test isfile(joinpath(dir, "field.vtm"))
        end
    end
end
