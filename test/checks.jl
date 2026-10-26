# Checks of computed values against exact identities and analytical solutions

gamma = 2.67513e-4

@testset "Finite element matrices" begin
    # Exact identities for P1 elements, which represent the coordinates exactly:
    # 1ᵀM1 = |Omega|, xᵀMy = ∫ x y, xᵀSy = ∫ ∇x ⋅ D ∇y = D_xy |Omega|, 1ᵀMₓ1 = ∫ x
    D = [3.0 1.0 0.0; 1.0 2.0 0.5; 0.0 0.5 1.0]
    for (setup, D_i) in [
        (DiskSetup(; ncell = 2, ecs = BoxECS(; margin = 1.0)), D[1:2, 1:2]),
        (SphereSetup(; ncell = 2, nsidewall = 60, ecs = BoxECS(; margin = 1.0)), D),
    ]
        Random.seed!(2)
        mesh, = create_geometry(setup)
        dim = size(mesh.points[1], 1)
        for i in eachindex(mesh.points)
            p, e = mesh.points[i], mesh.elements[i]
            V = sum(first(SpinDoctor.element_geometry(p, e)))
            M = SpinDoctor.assemble_mass_matrix(e, p)
            S = SpinDoctor.assemble_stiffness_matrix(e, p, D_i)
            @test M == M' && S == S'
            @test sum(M) ≈ V
            for a = 1:dim, b = 1:dim
                xa, xb = p[a, :], p[b, :]
                @test xa' * S * xb ≈ D_i[a, b] * V atol = 1e-10 * norm(D_i) * V
            end
            # ∫ x y = (xᵀMy), and with x as a P1 weight: ∫ x phi_i phi_j
            Mx = SpinDoctor.assemble_mass_matrix(e, p; nodeweights = p[1, :])
            @test sum(Mx) ≈ ones(size(p, 2))' * M * p[1, :]
            @test p[2, :]' * Mx * ones(size(p, 2)) ≈ p[1, :]' * M * p[2, :]
            # Piecewise constant weights
            w = fill(2.0, size(e, 2))
            @test SpinDoctor.assemble_mass_matrix(e, p; elementweights = w) ≈ 2M
        end
    end
end

@testset "Compartment boundaries" begin
    # Every boundary facet of a compartment is a face of one of its elements. The
    # saved meshes of bent and layered cylinders used to violate this (and could
    # not be assembled), since facets with all their points on an interface were
    # also assigned to the neighboring compartment.
    for (name, setup) in [
        "cylinders_bent_hull" => CylinderSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = pi / 4,
        ),
        "cylinders_layers_hull" => CylinderSetup(;
            ncell = 2,
            layersizes = [0.6, 1.0],
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = pi / 4,
        ),
    ]
        savedir = joinpath(@__DIR__, "regression", "meshes", name)
        mesh, = create_geometry(setup; savedir, recreate = false)
        for i in eachindex(mesh.points), f in mesh.facets[i, :]
            isempty(f) || @test all(>(0), SpinDoctor.facet_elements(mesh.elements[i], f))
        end
        model = Model(; mesh, get_coeffs(setup)...)
        @test assemble_matrices(model) isa NamedTuple
    end
end

@testset "Outward normals" begin
    # Divergence theorem: ∫_∂Omega x_d n_d dS = |Omega| for each compartment Omega and axis d.
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
            twist = pi / 4,
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
    for setup in setups
        Random.seed!(1)
        mesh, = create_geometry(setup)
        model = Model(; mesh, get_coeffs(setup)...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        for i in eachindex(mesh.points), d = 1:dim
            @test mesh.points[i][d, :]' * matrices.G[i][:, d] ≈ matrices.volumes[i] rtol =
                1e-8
        end

        # With equal spin densities, permeable interfaces conserve magnetization
        # and the flux matrix is symmetric
        @test norm(matrices.Q * ones(size(matrices.Q, 1))) <= 1e-10 * norm(matrices.Q, 1)
        @test matrices.Q ≈ matrices.Q'
    end
end

@testset "Permeable plates vs analytical solution" begin
    # For a gradient across the plates, the 2D plates and the 3D slabs reduce to
    # the one-dimensional problem of the analytical matrix formalism
    profile = PGSE(2000.0, 6000.0)
    for setup in [
        PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement = 0.5),
        SlabSetup(; depth = 20.0, widths = fill(5.0, 4), height = 5.0, refinement = 2.0),
    ]
        coeffs = get_coeffs(setup)
        mesh, = create_geometry(setup)
        model = Model(; mesh, coeffs...)
        matrices = assemble_matrices(model)
        dim = size(mesh.points[1], 1)
        dir = [1.0; zeros(dim - 1)]
        gradient = ScalarGradient(dir, profile; bvalue = 1000, gamma)
        u = solve_btpde(model, matrices, gradient, IntervalConstantSolver(; timestep = 5.0))
        volumes = compartment_volumes(mesh)
        signal = compute_signal(matrices.M, u) / sum(volumes)
        acoeffs = analytical_coefficients(setup, coeffs)
        laplace = solve_analytical_laplace(;
            acoeffs...,
            eiglim = length2eig(0.3, 0.002),
            eigstep = 1e-8,
        )
        signal_analytical = solve_analytical_mf(laplace, volumes, gradient) / sum(volumes)
        @test real(signal) ≈ real(signal_analytical) rtol = 1e-2
    end
end

@testset "Compartment meshes and box faces" begin
    for setup in [
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
            ecs = BoxECS(; margin = 1.0),
            height = 10.0,
        ),
    ]
        Random.seed!(3)
        mesh, = create_geometry(setup)
        dim = size(mesh.points[1], 1)
        volumes = compartment_volumes(mesh)

        # The ECS as a mesh and model of its own
        ecs = compartment_mesh(mesh, length(mesh.points))
        @test compartment_volumes(ecs) ≈ volumes[end:end]
        model = Model(;
            mesh = ecs,
            rho = [1.0 + 0im],
            D = [0.002 * Matrix(I, dim, dim)],
            T2 = [Inf],
            kappa = zeros(size(ecs.facets, 2)),
            gamma = gamma,
        )
        matrices = assemble_matrices(model)
        @test sum(matrices.M) ≈ volumes[end]
        @test compartment_weights(model) ≈ [1.0]

        # Faces of the box: the faces crossed by no cell are full rectangles
        faces = box_faces(mesh)
        p = mesh.points[end]
        L = vec(maximum(p; dims = 2) - minimum(p; dims = 2))
        area(f) = sum(first(SpinDoctor.facet_geometry(p, f)))
        ncomplete = setup isa ExtrusionSetup ? 2 : dim
        for d = 1:ncomplete, face in (faces[2d-1], faces[2d])
            @test area(face) ≈ prod(L[[1:(d-1); (d+1):dim]])
        end
        if setup isa ExtrusionSetup
            # The cylinders cross the top and bottom faces
            @test 0 < area(faces.zmin) < L[1] * L[2]
            @test area(faces.zmin) ≈ area(faces.zmax)
        end
    end

    # Given cell configuration
    setup = DiskSetup(; ncell = 2, nsidewall = 40, ecs = BoxECS(; margin = 1.0))
    cells = (; centers = [-4.0 4.0; 0.0 0.0], radii = [3.0 2.0])
    mesh, = create_geometry(setup; cells)
    volumes = compartment_volumes(mesh)
    @test volumes[1:2] ≈ pi * [3.0, 2.0] .^ 2 rtol = 0.02
    @test_throws ErrorException create_geometry(SphereSetup(; ncell = 3); cells)
end

@testset "Geometry files" begin
    setup = CylinderSetup(; ncell = 2, rmin = 2.0, rmax = 3.0, ecs = BoxECS(), height = 5.0)
    other = CylinderSetup(; ncell = 2, rmin = 2.0, rmax = 3.0, ecs = BoxECS(), height = 6.0)
    same(m1, m2) =
        m1.points == m2.points && m1.elements == m2.elements && m1.facets == m2.facets
    mktempdir() do dir
        mesh, _, cells = create_geometry(setup; rng = Xoshiro(1))

        # JLD2 cache: created, then loaded (also with another random generator)
        file = joinpath(dir, "geometry.jld2")
        @test same(create_geometry(setup; cache = file, rng = Xoshiro(1))[1], mesh)
        @test isfile(file)
        mesh_cached, _, cells_cached =
            create_geometry(setup; cache = file, rng = Xoshiro(2))
        @test same(mesh_cached, mesh)
        @test cells_cached == cells

        # A different setup is not taken from the cache
        @test !same(create_geometry(other; cache = file, rng = Xoshiro(1))[1], mesh)
        @test load_geometry(file).setup == repr(other)

        # Text files: created, then reused (with full precision)
        savedir = joinpath(dir, "text")
        @test same(create_geometry(setup; savedir, rng = Xoshiro(1))[1], mesh)
        @test same(create_geometry(setup; savedir, rng = Xoshiro(2))[1], mesh)
        @test !same(create_geometry(other; savedir, rng = Xoshiro(1))[1], mesh)
    end
end

@testset "Utilities" begin
    # Directions on a hemisphere
    d = unitsphere(20; half = true, normal = [1.0, 1.0, 0.0])
    @test size(d) == (3, 20)
    @test all(>=(0), [1, 1, 0]' * d)
    @test all(≈(1), sum(abs2, d; dims = 1))

    # Reproducible cell configurations
    setup = SphereSetup(; ncell = 3)
    @test SpinDoctor.create_cells(setup; rng = Xoshiro(1)) ==
          SpinDoctor.create_cells(setup; rng = Xoshiro(1))

    # VTK output of 2D (triangle) and 3D (tetrahedron) meshes
    for setup in [
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
            u = solve_btpde(model, matrices, gradient, solver; callbacks = [writer])
            @test isfile(joinpath(dir, "nested", "output", "solution.pvd"))
            save_field(joinpath(dir, "field"), mesh, u)
            @test isfile(joinpath(dir, "field.vtm"))
        end
    end
end
