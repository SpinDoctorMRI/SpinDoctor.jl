# Record regression outputs of SpinDoctor for a set of small setups.
#
# Usage: julia --project=<env> record.jl <outputs.toml> [generate]
#
# With `generate`, the geometries are created from fixed seeds and saved to
# `meshes/<case>` next to this script. Without it, the saved geometries are read
# back, so the old and the new code start from identical meshes.

using SpinDoctor
using LinearAlgebra
using Random
using TOML
using OrdinaryDiffEqBDF: QNDF

outfile = ARGS[1]
generate = length(ARGS) >= 2 && ARGS[2] == "generate"
meshdir = joinpath(@__DIR__, "meshes")

const gamma = 2.67513e-4
const D0 = 0.002

function cellcoeffs(setup, nlayer, dim)
    coefficients(
        setup;
        D = (; cell = [D0 * I(dim) for _ = 1:nlayer], ecs = D0 * I(dim)),
        T2 = (; cell = fill(Inf, nlayer), ecs = Inf),
        rho = (; cell = fill(1.0, nlayer), ecs = 1.0),
        kappa = (;
            cell_interfaces = fill(1e-4, nlayer - 1),
            cell_boundaries = fill(0, nlayer),
            cell_ecs = 1e-4,
            ecs = 0,
        ),
        gamma,
    )
end

function platecoeffs(setup, n, dim)
    coefficients(
        setup;
        D = [D0 * I(dim) for _ = 1:n],
        T2 = fill(Inf, n),
        rho = fill(1.0, n),
        kappa = (; interfaces = fill(1e-4, n - 1), boundaries = fill(0.0, n)),
        gamma,
    )
end

neuroncoeffs(setup) = coefficients(
    setup;
    D = (; neuron = D0 * I(3), ecs = D0 * I(3)),
    T2 = (; neuron = Inf, ecs = Inf),
    rho = (; neuron = 1.0, ecs = 1.0),
    kappa = (; neuron_ecs = 1e-4, neuron = 0.0, ecs = 0.0),
    gamma,
)

# name => (seed, setup, coeffs, analytical)
cases = [
    "slabs" => begin
        s = SlabSetup(; depth = 20.0, widths = fill(5.0, 4), height = 20.0, refinement = 10.0)
        (1, s, platecoeffs(s, 4, 3), true)
    end,
    "plates_2d" => begin
        s = PlateSetup(; depth = 20.0, widths = fill(5.0, 4), refinement = 2.0)
        (1, s, platecoeffs(s, 4, 2), true)
    end,
    "disks_box_2d" => begin
        s = DiskSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            ecs = BoxECS(; margin = 1.0),
            refinement = 0.5,
        )
        (2, s, cellcoeffs(s, 1, 2), false)
    end,
    "cylinders_box" => begin
        s = CylinderSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            dmin = 0.2,
            dmax = 0.3,
            ecs = BoxECS(; margin = 1.0),
            height = 10.0,
            refinement = 1.0,
        )
        (3, s, cellcoeffs(s, 1, 3), false)
    end,
    "cylinders_bent_hull" => begin
        s = CylinderSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = pi / 4,
        )
        (4, s, cellcoeffs(s, 1, 3), false)
    end,
    "cylinders_layers_hull" => begin
        s = CylinderSetup(;
            ncell = 2,
            layersizes = [0.6, 1.0],
            rmin = 2.0,
            rmax = 4.0,
            ecs = ConvexHullECS(; margin = 1.0),
            height = 10.0,
            bend = 0.02,
            twist = pi / 4,
        )
        (4, s, cellcoeffs(s, 2, 3), false)
    end,
    "spheres_box" => begin
        s = SphereSetup(;
            ncell = 3,
            rmin = 2.0,
            rmax = 4.0,
            nsidewall = 60,
            ecs = BoxECS(; margin = 1.0),
        )
        (5, s, cellcoeffs(s, 1, 3), false)
    end,
    "neuron_box" => begin
        s = NeuronSetup(; ecs = BoxECS(; margin = 1.0))
        (6, s, neuroncoeffs(s), false)
    end,
]

# Deterministic test vector
testvec(n) = sin.(1:n)

re(z) = real.(z)
im_(z) = imag.(z)

function gradient_for(profile, dir, b)
    ScalarGradient(dir, profile; bvalue = b, gamma)
end

directions(dim) =
    dim == 2 ? [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]] :
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]

function tryrecord!(out, key, f)
    t = @elapsed try
        out[key] = f()
    catch e
        out[key] = "ERROR: " * sprint(showerror, e)[1:min(end, 300)]
    end
    println("  ", rpad(key, 24), round(t; digits = 2), " s")
end

results = Dict{String,Any}()
results["julia_version"] = string(VERSION)

for (name, (seed, setup, coeffs, analytical)) in cases
    println(name)
    savedir = joinpath(meshdir, name)
    if setup isa NeuronSetup
        mkpath(savedir)
        cp(joinpath(@__DIR__, "neuron.msh"), joinpath(savedir, "mesh.msh"); force = true)
    end
    if generate
        Random.seed!(seed)
        create_geometry(setup; savedir, recreate = true)
    end
    mesh, = create_geometry(setup; savedir, recreate = false)
    model = Model(; mesh, coeffs...)
    dim = size(mesh.points[1], 1)
    ncmpt = length(mesh.points)
    out = Dict{String,Any}()
    results[name] = out

    out["npoint"] = size.(mesh.points, 2)
    out["nelement"] = size.(mesh.elements, 2)
    out["volumes"] = compartment_volumes(mesh)

    matrices = nothing
    t = @elapsed try
        matrices = assemble_matrices(model)
    catch e
        out["assemble_matrices"] = "ERROR: " * sprint(showerror, e)[1:min(end, 300)]
    end
    println("  assemble_matrices       ", round(t; digits = 2), " s")
    isnothing(matrices) && continue
    v = testvec(size(matrices.M, 1))
    for key in (:M, :S, :R, :Q)
        A = getfield(matrices, key)
        out["matrix_$(key)_sum"] = sum(A)
        out["matrix_$(key)_vAv"] = v' * A * v
    end
    for d = 1:dim
        out["matrix_Mx$(d)_sum"] = sum(matrices.Mx[d])
        out["matrix_Mx$(d)_vAv"] = v' * matrices.Mx[d] * v
    end
    out["G_sum"] = [sum(matrices.G[i][:, d]) for i = 1:ncmpt for d = 1:dim]
    out["G_xn"] =
        [mesh.points[i][d, :]' * matrices.G[i][:, d] for i = 1:ncmpt for d = 1:dim]

    profile = PGSE(2000.0, 6000.0)
    dirs = directions(dim)
    grad = gradient_for(profile, dirs[end], 1000.0)

    tryrecord!(
        out,
        "adc_sta",
        () -> reduce(
            vcat,
            [compute_adc_sta(model, gradient_for(profile, d, 1000.0)) for d in dirs],
        ),
    )
    tryrecord!(
        out,
        "hadc",
        () -> reduce(
            vcat,
            [solve_hadc(model, matrices, gradient_for(profile, d, 1000.0)) for d in dirs],
        ),
    )

    tryrecord!(
        out,
        "btpde_qndf",
        () -> begin
            u = solve_btpde(model, matrices, grad, QNDF())
            s = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
            Dict("re" => re(s), "im" => im_(s))
        end,
    )
    tryrecord!(
        out,
        "btpde_cn",
        () -> begin
            solver = ThetaSolver(; timestep = 20.0)
            u = solve_btpde(model, matrices, grad, solver)
            s = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
            s0 = compute_signal.(
                matrices.M_cmpts,
                split_field(mesh, initial_conditions(model)),
            )
            Dict(
                "re" => re(s),
                "im" => im_(s),
                "fit_adc" => [fit_adc([0.0, 1000.0], [s0[i], s[i]]) for i = 1:ncmpt],
            )
        end,
    )
    if name == "slabs"
        TE = 5000.0
        phi = -pi / 6
        Rot = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]
        gvec(t) = 1.0 * Rot * [sin(2pi * t / TE), sin(20pi * t / TE) / 5, cos(2pi * t / TE)]
        general = GeneralGradient(; gvec, TE)
        tryrecord!(
            out,
            "btpde_general",
            () -> begin
                u = solve_btpde(model, matrices, general, QNDF())
                s = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
                Dict("re" => re(s), "im" => im_(s))
            end,
        )
        ogse = gradient_for(CosOGSE(5000.0, 5000.0, 2), dirs[end], 1000.0)
        tryrecord!(out, "hadc_ogse", () -> solve_hadc(model, matrices, ogse))
    end

    tryrecord!(
        out,
        "laplace_mf",
        () -> begin
            lap_eig = solve_laplace(model, matrices; neig_max = 30)
            u = solve_mf(model, matrices, lap_eig, grad)
            s = compute_signal.(matrices.M_cmpts, split_field(mesh, u))
            Dict("eigenvalues" => lap_eig.values, "re" => re(s), "im" => im_(s))
        end,
    )

    tryrecord!(
        out,
        "karger",
        () -> begin
            kdirs = dim == 2 ? unitcircle(8)[1:2, :] : unitsphere(10)
            adcs = [
                solve_hadc(model, matrices, gradient_for(profile, collect(d), 1000.0)) for d in eachcol(kdirs)
            ]
            difftensors = fit_tensors(kdirs, adcs)
            s = solve_karger(model, difftensors, grad; timestep = 5.0)
            Dict(
                "difftensors" => reduce(vcat, vec.(difftensors)),
                "re" => re(s),
                "im" => im_(s),
            )
        end,
    )

    if analytical
        tryrecord!(
            out,
            "analytical",
            () -> begin
                volumes = compartment_volumes(mesh)
                D_avg = 1 / dim * tr.(model.D)' * volumes / sum(volumes)
                eiglim = length2eig(1.0, D_avg)
                acoeffs = analytical_coefficients(setup, coeffs)
                alap = solve_analytical_laplace(; acoeffs..., eiglim, eigstep = 1e-8)
                s = solve_analytical_mf(alap, volumes, grad)
                Dict("re" => real(s), "im" => imag(s))
            end,
        )
    end
end

open(outfile, "w") do io
    TOML.print(io, results; sorted = true)
end
println("Wrote ", outfile)
