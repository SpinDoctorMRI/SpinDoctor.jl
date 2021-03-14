using SpinDoctor

## Choose setup script
# include("setups/cylinders.jl")
# include("setups/spheres.jl")
include("setups/neuron.jl")

## Prepare experiments
prepare_pde!(setup)
mesh, surfaces, cells = create_geometry(setup)
volumes = get_cmpt_volumes(mesh)
σ_avg = setup.pde[:σ]' * volumes / sum(volumes)

# Sizes
ncompartment = length(mesh.points)
namplitude = length(setup.gradient[:values])
nsequence = length(setup.gradient[:sequences])
ndirection = size(setup.gradient[:directions], 2)

# Q-values and b-values
if setup.gradient[:values_type] == 'q'
    qvalues = repeat(setup.gradient[:values], 1, nsequence)
    bvalues = setup.gradient[:values] .^ 2 .* bvalue_no_q.(setup.gradient[:sequences])'
else
    bvalues = repeat(setup.gradient[:values], 1, nsequence)
    qvalues = .√(setup.gradient[:values] ./ bvalue_no_q.(setup.gradient[:sequences])')
end

##
btpde = @time solve_btpde(mesh, setup)

## Solve BTPDE
if !isnothing(setup.btpde)
    btpde = @time solve_btpde(mesh, setup)

    name = split(setup.name, "/")[end]
    if setup.btpde[:nsave] == 1
        savefield(
            mesh,
            btpde.magnetization[:, end, 1, 1],
            "output/$name/magnetization_btpde",
        )
    else
        save_btpde_results(mesh, btpde, setup, "output/$name/magnetization_btpde")
    end

    adc = [
        fit_adc(
            bvalues[:, iseq],
            real(btpde.signal[icmpt, :, iseq, idir]) / (setup.pde[:ρ]' * volumes),
        ) for idir = 1:ndirection, iseq = 1:nsequence, icmpt = 1:ncompartment
    ]
end

## Solve MF
if !isnothing(setup.mf)
    λ_max = length2eig(setup.mf[:length_scale], σ_avg)
    lap_eig = compute_laplace_eig(mesh, setup.pde, λ_max, setup.mf[:neig_max])
    length_scales = eig2length.(lap_eig.values, σ_avg)
    mf = solve_mf(mesh, setup, lap_eig)

    name = split(setup.name, "/")[end]
    savefield(mesh, mf.magnetization[:, 1, 1, 1], "output/$name/magnetization_mf")

    npoint_cmpts = size.(mesh.points, 2)
    bounds = cumsum([0; npoint_cmpts])
    ϕ_cmpts = [lap_eig.funcs[bounds[i]+1:bounds[i+1], :] for i = 1:mesh.ncompartment]
    savefield(mesh, ϕ_cmpts, "output/$name/laplace_eig", "Laplace eigenfunction")
end
