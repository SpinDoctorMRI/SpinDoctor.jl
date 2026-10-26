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
        u = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 5.0))
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

        # The analytical solution converges with the number of eigenmodes (the closely
        # spaced eigenvalues of weakly coupled layers used to be missed)
        coarse = solve_analytical_laplace(;
            acoeffs...,
            eiglim = length2eig(1.0, 0.002),
            eigstep = 1e-8,
        )
        @test length(coarse.alpha) < length(laplace.alpha)
        @test solve_analytical_mf(coarse, volumes, gradient) / sum(volumes) ≈
              signal_analytical rtol = 1e-4
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

@testset "Packings and periodic samples" begin
    # Periodic packing: requested fraction, gaps between all pairs (minimum image)
    rng = Xoshiro(5)
    radii = 0.5 .+ rand(rng, 20)
    cells = periodic_packing(radii; fraction = 0.75, gap = 0.05, rng)
    L = cells.box[2][1]
    @test sum(pi .* cells.radii .^ 2) / L^2 ≈ 0.75
    gaps = [
        norm(mod.(cells.centers[:, i] - cells.centers[:, j] .+ L / 2, L) .- L / 2) -
        cells.radii[i] - cells.radii[j] for j = 2:20 for i = 1:(j-1)
    ]
    @test minimum(gaps) >= 0.9 * 0.05

    # Disks cut by the box: the pieces of each cell (with its periodic images) make up the
    # whole cell polygon (up to vertices moved onto the faces), and the opposite faces of
    # the box match
    setup = DiskSetup(;
        ncell = 20,
        rmin = 0.5,
        rmax = 1.5,
        nsidewall = 16,
        layersizes = [0.6, 1.0],
        ecs = PeriodicBoxECS(),
    )
    mesh, = create_geometry(setup; cells)
    V = compartment_volumes(mesh)
    nside(r) = max(8, round(Int, 16 * r))
    polyarea(r) = nside(r) / 2 * r^2 * sin(2pi / nside(r))
    r = vec(cells.radii)
    @test sum(V) ≈ L^2
    @test V[1:20] ≈ polyarea.(0.6 .* r) rtol = 1e-2
    @test V[21:40] ≈ polyarea.(r) - polyarea.(0.6 .* r) rtol = 1e-2
    faces = box_faces(mesh)
    len(f) = sum(first(SpinDoctor.facet_geometry(mesh.points[end], f)))
    @test len(faces.xmin) ≈ len(faces.xmax)
    @test len(faces.ymin) ≈ len(faces.ymax)

    # Without a given configuration, a packing is created
    setup = CylinderSetup(;
        ncell = 10,
        rmin = 1.0,
        rmax = 2.0,
        nsidewall = 16,
        ecs = PeriodicBoxECS(; fraction = 0.6),
        height = 2.0,
    )
    mesh, _, cells = create_geometry(setup; rng = Xoshiro(6))
    V = compartment_volumes(mesh)
    r = vec(cells.radii)
    nside16(r) = max(8, round(Int, 16 * r / 1.5))
    @test sum(V) ≈ 2 * cells.box[2][1]^2
    @test V[1:10] ≈ 2 * nside16.(r) / 2 .* r .^ 2 .* sin.(2pi ./ nside16.(r)) rtol = 1e-2
    @test sum(pi * r .^ 2) / cells.box[2][1]^2 ≈ 0.6

    # Voronoi cells: ECS fraction, periodic faces, and outward normals
    setup = VoronoiSetup(; ncell = 8, boxsize = 10.0, fraction = 0.8)
    mesh, = create_geometry(setup; rng = Xoshiro(1))
    V = compartment_volumes(mesh)
    @test sum(V) ≈ 1000
    @test V[end] / sum(V) ≈ 0.2 rtol = 1e-6
    faces = box_faces(mesh)
    area(f) = sum(first(SpinDoctor.facet_geometry(mesh.points[end], f)))
    @test area(faces.xmin) ≈ area(faces.xmax)
    @test area(faces.ymin) ≈ area(faces.ymax)
    @test area(faces.zmin) ≈ area(faces.zmax)
    model = Model(; mesh, get_coeffs(SphereSetup(; ncell = 8, ecs = BoxECS()))...)
    matrices = assemble_matrices(model)
    for i in eachindex(mesh.points), d = 1:3
        @test mesh.points[i][d, :]' * matrices.G[i][:, d] ≈ matrices.volumes[i] rtol = 1e-8
    end
end

@testset "Periodic boundary conditions" begin
    # Free diffusion in a periodic rectangle, with boundary points every 0.25: the periodic
    # part of the magnetization stays constant, and the signal is exp(-b u' D u)
    Lx, Ly, h = 4.0, 3.0, 0.25
    corners = [[0.0, 0.0], [Lx, 0.0], [Lx, Ly], [0.0, Ly]]
    boundary = reduce(vcat, map(zip(corners, circshift(corners, -1))) do (a, b)
        n = round(Int, norm(b - a) / h)
        [a + (b - a) * i / n for i = 0:(n-1)]
    end)
    n = length(boundary)
    surfaces = (;
        points = reduce(hcat, boundary),
        facets = [1:n circshift(1:n, -1)]',
        facetmarkers = ones(Int, n),
        regions = [1.0; 1.0;;],
        box = ([0.0, 0.0], [Lx, Ly]),
    )
    mesh = SpinDoctor.split_mesh(SpinDoctor.create_mesh(surfaces, 0.05))
    D = [0.002 0.0005; 0.0005 0.001]
    model = Model(; mesh, rho = [1.0 + 0im], D = [D], T2 = [Inf], kappa = [0.0], gamma)
    matrices = assemble_matrices(model; periodic = true)
    @test size(matrices.P, 2) == size(mesh.points[1], 2) - size(boundary, 1) / 2 - 1
    @test norm(matrices.P' * matrices.G[1]) < 1e-12
    for d = 1:2
        # The derivative of a constant is zero, and C' 1 is a boundary integral
        @test norm(matrices.Cx[d] * ones(size(matrices.M, 1))) < 1e-14
        @test matrices.Cx[d]' * ones(size(matrices.M, 1)) ≈ (matrices.G[1]*D)[:, d]
    end
    for (profile, dir) in
        ((PGSE(2000.0, 20000.0), [1.0, 0.0]), (CosOGSE(5000.0, 6000.0, 3), [1.0, 2.0]))
        gradient = ScalarGradient(dir, profile; bvalue = 1000, gamma)
        u = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 10.0))
        u0 = solve_btpde(
            model,
            assemble_matrices(model),
            gradient,
            ThetaSolver(; timestep = 10.0),
        )
        signal(u) = abs(compute_signal(matrices.M, u)) / (Lx * Ly)
        @test signal(u) ≈ exp(-1000 * gradient.dir' * D * gradient.dir) rtol = 1e-5
        @test signal(u0) > 1.5 * signal(u)
    end
    @test only(solve_hadc_tensors(model, matrices, PGSE(2000.0, 20000.0))) ≈ D

    # Meshes that are not periodic
    setup = SphereSetup(; ncell = 2, ecs = BoxECS())
    mesh, = create_geometry(setup; rng = Xoshiro(1))
    model = Model(; mesh, get_coeffs(setup)...)
    @test_throws ErrorException assemble_matrices(model; periodic = true)

    # A sample of periodic tissue: the signals and the HADC of the compartments do not
    # depend on the position of the cells in the box (unlike with impermeable box faces),
    # and the cells that do not touch the box are the same as without periodicity
    setup = DiskSetup(;
        ncell = 8,
        rmin = 1.0,
        rmax = 2.0,
        dmin = 0.2,
        nsidewall = 16,
        ecs = PeriodicBoxECS(; fraction = 0.5),
        refinement = 0.1,
    )
    coeffs = get_coeffs(setup)
    coeffs = (; coeffs..., kappa = zero(coeffs.kappa))
    cells = SpinDoctor.create_cells(setup; rng = Xoshiro(1))
    L = cells.box[2][1]
    shifted = (; cells..., centers = mod.(cells.centers .+ [L / 3, L / 5], L))
    gradient = ScalarGradient([1.0, 0.3], PGSE(5000.0, 50000.0); bvalue = 2000, gamma)
    results = map((cells, shifted)) do cells
        mesh, = create_geometry(setup; cells)
        model = Model(; mesh, coeffs...)
        map((true, false)) do periodic
            matrices = assemble_matrices(model; periodic)
            u = solve_btpde(model, matrices, gradient, ThetaSolver(; timestep = 200.0))
            signals = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
            adcs = solve_hadc(model, matrices, gradient; timestep = 200.0)
            (; signals = abs.(signals) ./ matrices.volumes, adcs)
        end
    end
    (periodic, walls), (periodic_shifted, walls_shifted) = results
    @test periodic.signals ≈ periodic_shifted.signals rtol = 5e-3
    @test periodic.adcs ≈ periodic_shifted.adcs rtol = 5e-3
    @test !isapprox(walls.signals, walls_shifted.signals; rtol = 1e-2)
    @test !isapprox(walls.adcs, walls_shifted.adcs; rtol = 1e-2)
    r = vec(cells.radii)
    inside = [all(r[j] .< cells.centers[:, j] .< L - r[j]) for j = 1:8]
    @test periodic.signals[1:8][inside] ≈ walls.signals[1:8][inside] rtol = 5e-3
    @test periodic.adcs[end] > 2 * walls.adcs[end]

    # Periodic meshes of cylinders (also along the axis) and polyhedral cells. Along the
    # cylinders, the diffusion is free with periodicity, and restricted by the top and
    # bottom faces without.
    setup = CylinderSetup(;
        ncell = 5,
        rmin = 1.0,
        rmax = 2.0,
        nsidewall = 12,
        ecs = PeriodicBoxECS(; fraction = 0.5),
        height = 2.0,
    )
    mesh, = create_geometry(setup; rng = Xoshiro(2))
    model = Model(; mesh, get_coeffs(setup)...)
    profile = PGSE(2000.0, 20000.0)
    hadc = solve_hadc_tensors(model, assemble_matrices(model; periodic = true), profile)
    @test [D[3, 3] for D in hadc] ≈ [D[3, 3] for D in model.D]
    hadc = solve_hadc_tensors(model, assemble_matrices(model), profile)
    @test all(i -> hadc[i][3, 3] < model.D[i][3, 3] / 2, eachindex(hadc))

    setup = VoronoiSetup(; ncell = 6, boxsize = 8.0, fraction = 0.8)
    mesh, = create_geometry(setup; rng = Xoshiro(2))
    model = Model(; mesh, get_coeffs(setup)...)
    hadc = solve_hadc_tensors(model, assemble_matrices(model; periodic = true), profile)
    hadc_walls = solve_hadc_tensors(model, assemble_matrices(model), profile)
    @test all(D -> all(>(0), eigvals(Symmetric(D))), hadc)
    @test tr(hadc[end]) > 2 * tr(hadc_walls[end])
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

        # Text files of periodic samples: cells with a gap and a box, and polygon facets
        setup = VoronoiSetup(; ncell = 4, boxsize = 8.0, fraction = 0.8)
        savedir = joinpath(dir, "voronoi")
        mesh, _, cells = create_geometry(setup; savedir, rng = Xoshiro(1))
        mesh_text, _, cells_text = create_geometry(setup; savedir)
        @test same(mesh_text, mesh)
        @test cells_text == cells
        cells = periodic_packing(1 .+ rand(Xoshiro(1), 5); fraction = 0.5)
        SpinDoctor.save_cells(joinpath(dir, "cells"), cells)
        @test SpinDoctor.read_cells(joinpath(dir, "cells")) == cells
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
            solver = ThetaSolver(; timestep = 5.0)
            u = solve_btpde(model, matrices, gradient, solver; callbacks = [writer])
            @test isfile(joinpath(dir, "nested", "output", "solution.pvd"))
            save_field(joinpath(dir, "field"), mesh, u)
            @test isfile(joinpath(dir, "field.vtm"))
        end
    end
end
