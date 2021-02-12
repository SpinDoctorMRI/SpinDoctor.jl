
include("src/SpinDoctor.jl")


## Choose setup script

include("setups/cylinders.jl")


## Prepare experiments

cells = SpinDoctor.create_cells(cellsetup)
domain = SpinDoctor.prepare_pde(cellsetup, domainsetup)

tri = SpinDoctor.create_surface_triangulation(cellsetup, cells, domain)

tetgenmesh = SpinDoctor.create_mesh(cellsetup, domain, tri)
if any(cellsetup.deformation .> 1e-16)
    SpinDoctor.deform_domain!(tetgenmesh.pointlist, cellsetup.deformation)
end
mesh = SpinDoctor.split_mesh(domain, tetgenmesh)

directions = SpinDoctor.create_directions(experiment)
volumes = SpinDoctor.get_cmpt_volumes(mesh)

σ_avg = domain.diffusivity' * volumes / sum(volumes)


## Solve

if !isnothing(experiment.btpde)
    btpde = @time SpinDoctor.solve_btpde(mesh, domain, experiment, directions)

    if experiment.btpde.nsave == 1
        SpinDoctor.savefield(mesh, btpde.magnetization[:, 1, 1, 1], "output/$(cellsetup.name)/magnetization_btpde")
    else
        SpinDoctor.save_btpde_results(mesh, btpde, experiment, directions, "output/$(cellsetup.name)/magnetization_btpde")
    end
end

if !isnothing(experiment.mf)
    λ_max = SpinDoctor.length2eig(experiment.mf.length_scale, σ_avg)
    lap_eig = SpinDoctor.compute_laplace_eig(mesh, domain, λ_max, experiment.mf.neig_max)
    length_scales = SpinDoctor.eig2length.(lap_eig.values, σ_avg)
    mf = SpinDoctor.solve_mf(mesh, domain, experiment, lap_eig, directions)

    SpinDoctor.savefield(mesh, mf.magnetization[:, 1, 1, 1], "output/$(cellsetup.name)/magnetization_mf")

    npoint_cmpts = size.(mesh.points, 2)
    bounds = cumsum([0; npoint_cmpts])
    ϕ_cmpts = [lap_eig.funcs[bounds[i]+1:bounds[i+1], :] for i = 1:mesh.ncmpt]
    SpinDoctor.savefield(mesh, ϕ_cmpts, "output/$(cellsetup.name)/laplace_eig", "Laplace eigenfunction")
end
