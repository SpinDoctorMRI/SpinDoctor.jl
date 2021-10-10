# # Custom workaround OrdinaryDiffEq: https://github.com/SciML/OrdinaryDiffEq.jl/blob/54fb35870fa402fc95d665cd5f9502e2759ea436/src/derivative_utils.jl#L606-L608
# # as proposed in: https://github.com/CliMA/ClimaCore.jl/blob/af6b01bc8f945ada0794147d8f2c7ca5303adf62/examples/hybrid/inertial_gravity_wave_implicit.jl#L411-L414
# struct CustomWFact
#     W
#     γref
# end
# Base.similar(W::CustomWFact) = deepcopy(W)
# function linsolve!(::Type{Val{:init}}, f, u0; kwargs...)
#     function _linsolve!(x, W, b, update_matrix = false; kwargs...)
#         W isa CustomWFact || error("Linsolve specialized for CustomWFact")
#         ldiv!(x, W.W, b)
#     end
# end

"""
    solve_btpde(model, matrices, experiment)

Solve the Bloch-Torrey partial differential equation using P1 finite elements.
"""
function solve_btpde(model::Model, matrices, experiment::Experiment)

    # Measure function evalutation time
    starttime = Base.time()

    # Extract input parameters
    @unpack mesh, D, T₂, ρ = model
    @unpack M, S, R, Mx, Q, M_cmpts = matrices
    @unpack directions, sequences, values, values_type = experiment.gradient
    @unpack reltol, abstol, nsave = experiment.btpde

    # Deduce sizes
    ncompartment = length(ρ)
    npoint_cmpts = size.(mesh.points, 2)
    inds_cmpts = cumsum([0; npoint_cmpts])
    ndirection = size(directions, 2)
    nsequence = length(sequences)
    namplitude = length(values)

    qvalues, bvalues = get_values(experiment.gradient)

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex.(ρ), npoint_cmpts)...)

    # Allocate output arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization =
        Array{Matrix{ComplexF64},4}(undef, ncompartment, namplitude, nsequence, ndirection)
    time = Array{Vector{Float64},3}(undef, namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    # Time dependent ODE function
    function Mdξ!(dξ, ξ, p, t)
        # @show t
        @unpack J, S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
        mul!(dξ, J, ξ)
        nothing
    end

    # Time independent ODE function, given Jacobian
    function Mdξ_constant!(dξ, ξ, p, t)
        # @show t
        J = p.J
        mul!(dξ, J, ξ)
        nothing
    end

    # Time dependent Jacobian of ODE function with respect to the state `ξ`
    function Jac!(J, ξ, p, t)
        # println("Jac at $t")
        @unpack S, Q, A, R, q, f = p
        @. J = -(S + Q + R + im * f(t) * q * A)
        nothing
    end

    # # TODO: propagate types
    # # TODO: look for caching packages
    # # TODO: detect change in f and recompute J and W
    # function Wfact(u, p, γ, t)
    #     println("Wfact at $t")
    #     @unpack J = p
    #     CustomWFact(lu(M .- γ .* J), Ref(0.0))
    # end
    # function Wfact!(W, u, p, γ, t)
    #     # println("Wfact! at $t")
    #     @unpack J, γ_cache, factorization_cache = p
    #     ind = findfirst(isapprox(γ), γ_cache)
    #     # if isnothing(ind)
    #     if !(γ ≈ W.γref[])
    #         println("Refact at γ = $γ, t = $t")
    #         # W = factorize(M .- γ .* J)
    #         lu!(W.W, M .- γ .* J)
    #         W.γref[] = γ
    #     else
    #         println("No fact at γ = $γ, t = $t")
    #     end
    #     #     push!(γ_cache, γ)
    #     #     push!(factorization_cache, W)
    #     # else
    #     #     W = factorization_cache[ind]
    #     # end
    #     W
    # end

    # Jacobian sparsity pattern
    jac_prototype = -(S + Q + im * sum(Mx))

    # Iterate over gradient amplitudes, time profiles and directions
    for idir = 1:ndirection, iseq = 1:nsequence, iamp = 1:namplitude

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Time profile
        f = sequences[iseq]
        interval = (0, echotime(f))
        if nsave == 1
            saveat = [interval[2]]
            # saveat = [interval[1], interval[2]]
        else
            saveat = LinRange(interval..., nsave)
        end

        # Gradient direction
        dir = directions[:, idir]

        # Display state of iterations
        @printf "Solving BTPDE with size %d\n" (sum(npoint_cmpts))
        @printf "  Direction %d of %d: g = [%.2f, %.2f, %.2f]\n" idir ndirection dir...
        @printf "  Sequence  %d of %d: f = %s\n" iseq nsequence f
        @printf "  Amplitude %d of %d: q = %g, b = %g\n" iamp namplitude q b

        # Gradient direction dependent finite element matrix
        A = dir' * Mx

        # ODE problem
        J = -(S + Q + im * q * A)
        # γ_cache = Float64[]
        # factorization_cache = LinearAlgebra.Factorization[]
        # p = (; J, S, Q, A, R, q, f, γ_cache, factorization_cache)
        p = (; J, S, Q, A, R, q, f)

        # Gather ODE function
        if all(constant_intervals(p.f)) # is_constant(p.f, t)
            func = Mdξ_constant!
            Jac!(p.J, ρ, p, 0)
            jac = (J, _, p, t) -> (J .= p.J; nothing)
            # jac = (J, _, p, t) -> (println("Jac called at $t"); J .= p.J; nothing)
        else
            func = Mdξ!
            jac = Jac!
        end
        odefunction = ODEFunction(
            func;
            mass_matrix = M,
            jac,
            jac_prototype, # = Wfact(ρ, p, 0.1, 0.0),
            # Wfact = Wfact!,
        )
        odeproblem = ODEProblem(odefunction, ρ, interval, p, progress = false)#, dtmax = 50)

        tstops = intervals(f)[2:end-1]

        # Callback for updating problem
        function affect!(integrator)
            @unpack u, p, t = integrator
            # print("    t = $t: constant time profile, f = $(p.f(t))")
            Jac!(p.J, u, p, t)
        end
        callback = PresetTimeCallback(tstops, affect!)

        # Measure solving time
        itertime = Base.time()

        # Solve ODE problem
        sol = solve(
            odeproblem,
            QNDF(); # , linsolve = linsolve!);
            saveat,
            reltol,
            abstol,
            callback,
        )

        itertimes[iamp, iseq, idir] = Base.time() - itertime

        # Extract solution
        time[iamp, iseq, idir] = sol.t
        ξ = hcat(sol.u...)

        if nsave == 1
            time[iamp, iseq, idir] = time[iamp, iseq, idir][[end]]
            ξ = ξ[:, [end]]
        end

        # Split solution into compartments
        for icmpt = 1:ncompartment
            inds = inds_cmpts[icmpt]+1:inds_cmpts[icmpt+1]

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = ξ[inds, :]

            # Integrate final magnetization over compartment
            signal[icmpt, iamp, iseq, idir] =
                sum(M_cmpts[icmpt] * ξ[inds, end], dims = 1)[1]

        end
    end

    signal_allcmpts = sum(signal, dims = 1)[1, :, :, :]

    totaltime = Base.time() - starttime

    (; magnetization, signal, signal_allcmpts, time, itertimes, totaltime)
end
