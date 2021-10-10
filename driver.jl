using SpinDoctor
using LinearAlgebra


## Choose setup script
# include("setups/1axon_analytical.jl")
# include("setups/1sphere_analytical.jl")
include("setups/cylinders.jl")
# include("setups/spheres.jl")
# include("setups/neuron.jl")


## Prepare experiments
model = create_model(setup)
volumes = get_cmpt_volumes(model.mesh)
D_avg = 1 / 3 * tr.(model.D)' * volumes / sum(volumes)

# Sizes
ncompartment = length(model.mesh.points)
namplitude = length(experiment.gradient.values)
nsequence = length(experiment.gradient.sequences)
ndirection = size(experiment.gradient.directions, 2)

qvalues, bvalues = get_values(experiment.gradient)

## Assemble finite element matrices
matrices = assemble_matrices(model)

##
btpde = @time solve_btpde(model, matrices, experiment)

## Use manual time stepping scheme (theta rule)
if !isnothing(experiment.btpde_midpoint)
    btpde = @time solve_btpde_midpoint(model, matrices, experiment)
end

## Solve BTPDE
if !isnothing(experiment.btpde)
    # btpde = @time solve_btpde(model, experiment)

    refinement_str = ""
    if !isnothing(setup.refinement)
        refinement_str = "_refinement$(setup.refinement)"
    end
    output_dir = "output/" * setup.name * refinement_str
    isdir(output_dir) || mkpath(output_dir)

    if experiment.btpde.nsave == 1
        savefield(
            model.mesh,
            btpde.magnetization[:, end, 1, 1],
            "$output_dir/magnetization_btpde",
        )
    else
        save_btpde_results(model.mesh, experiment, btpde, "$output_dir/magnetization_btpde")
    end

    adc = [
        fit_adc(
            bvalues[:, iseq],
            real(btpde.signal[icmpt, :, iseq, idir]) / (model.ρ' * volumes),
        ) for idir = 1:ndirection, iseq = 1:nsequence, icmpt = 1:ncompartment
    ]
end


## Solve MF
if !isnothing(experiment.mf)
    λ_max = length2eig(experiment.mf.length_scale, D_avg)
    lap_eig = compute_laplace_eig(model, matrices, λ_max, experiment.mf.neig_max)
    length_scales = eig2length.(lap_eig.values, D_avg)
    mf = solve_mf(model, matrices, lap_eig, experiment)

    output_dir = "output/$(setup.name)"
    isdir(output_dir) || mkpath(output_dir)
    savefield(model.mesh, mf.magnetization[:, 1, 1, 1], "$output_dir/magnetization_mf")

    npoint_cmpts = size.(model.mesh.points, 2)
    bounds = cumsum([0; npoint_cmpts])
    ϕ_cmpts = [lap_eig.funcs[bounds[i]+1:bounds[i+1], :] for i = 1:ncompartment]
    savefield(model.mesh, ϕ_cmpts, "$output_dir/laplace_eig", "Laplace eigenfunction")
end

## Solve HADC
if !isnothing(experiment.hadc)
    hadc = solve_hadc(model, experiment)
end

## Solve Karger model
if !isnothing(experiment.karger)
    # Fit difftensors (assuming HADC has been solved)
    difftensors = fit_tensors(experiment.gradient.directions, hadc.adc)

    # Solv Karger
    karger = solve_karger(model, experiment, difftensors)
end

## Solve analytical model
if !isnothing(experiment.analytical)
    analytical_signal = solve_analytical(setup, experiment, volumes)
end
