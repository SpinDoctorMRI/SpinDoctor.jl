
using SpinDoctor
using DifferentialEquations: Trapezoid

## Choose setup script

include("setups/cylinders.jl")


## Prepare experiments

cells = create_cells(cellsetup)
domain = prepare_pde(cellsetup, domainsetup)

tri = create_surface_triangulation(cellsetup, cells, domain)

tetgenmesh = create_mesh(cellsetup, domain, tri)
if any(cellsetup.deformation .> 1e-16)
    deform_domain!(tetgenmesh.pointlist, cellsetup.deformation)
end
mesh = split_mesh(domain, tetgenmesh)

directions = create_directions(experiment)
volumes = get_cmpt_volumes(mesh)

σ_avg = domain.diffusivity' * volumes / sum(volumes)


## Solve

if !isnothing(experiment.btpde)
    btpde = @time solve_btpde(mesh, domain, experiment, directions)

    if experiment.btpde.nsave == 1
        savefield(mesh, btpde.magnetization[:, 1, 1, 1], "output/$(cellsetup.name)/magnetization_btpde")
    else
        save_btpde_results(mesh, btpde, experiment, directions, "output/$(cellsetup.name)/magnetization_btpde")
    end
end

if !isnothing(experiment.mf)
    λ_max = length2eig(experiment.mf.length_scale, σ_avg)
    lap_eig = compute_laplace_eig(mesh, domain, λ_max, experiment.mf.neig_max)
    length_scales = eig2length.(lap_eig.values, σ_avg)
    mf = solve_mf(mesh, domain, experiment, lap_eig, directions)

    savefield(mesh, mf.magnetization[:, 1, 1, 1], "output/$(cellsetup.name)/magnetization_mf")

    npoint_cmpts = size.(mesh.points, 2)
    bounds = cumsum([0; npoint_cmpts])
    ϕ_cmpts = [lap_eig.funcs[bounds[i]+1:bounds[i+1], :] for i = 1:mesh.ncompartment]
    savefield(mesh, ϕ_cmpts, "output/$(cellsetup.name)/laplace_eig", "Laplace eigenfunction")
end
