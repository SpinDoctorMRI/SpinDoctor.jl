"""
    solve_mf(model, matrices, lap_eig, experiment)

Solve for magnetization using Matrix Formalism.
"""
function solve_mf(model, matrices, lap_eig, experiment)

    # Measure time of function evaluation
    starttime = Base.time()

    # Extract parameters
    @unpack mesh, ρ = model
    @unpack M, M_cmpts = matrices
    @unpack directions, sequences, values, values_type = experiment.gradient
    @unpack ninterval = experiment.mf

    # Laplace eigenmodes
    λ = lap_eig.values
    ϕ = lap_eig.funcs
    Ax = lap_eig.moments
    T = lap_eig.massrelax

    # Sizes
    ndirection = size(directions, 2)
    namplitude = length(values)
    nsequence = length(sequences)
    npoint_cmpts = size.(mesh.points, 2)
    ncompartment = length(mesh.points)

    qvalues, bvalues = get_values(experiment.gradient)

    # Compartment index ranges
    inds_start = cumsum([1; npoint_cmpts[1:end-1]])
    inds_stop = cumsum(npoint_cmpts[1:end])
    get_inds(icmpt) = inds_start[icmpt]:inds_stop[icmpt]

    # Global mass matrix
    M = blockdiag(M_cmpts...)

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex.(ρ), npoint_cmpts)...)

    # Project initial spin density onto Laplace eigenfunction basis
    ν₀ = ϕ' * (M * ρ)

    # Allocate arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization =
        Array{Matrix{ComplexF64},4}(undef, ncompartment, namplitude, nsequence, ndirection)
    itertimes = zeros(namplitude, nsequence, ndirection)

    # Laplace operator in Laplace eigenfunction basis
    L = diagm(λ)

    # Iterate over gradient amplitudes, time profiles and directions
    for iamp = 1:namplitude, iseq = 1:nsequence, idir = 1:ndirection

        # Measure time of iteration
        itertime = time()

        # Gradient amplitude
        q = qvalues[iamp, iseq]
        b = bvalues[iamp, iseq]

        # Time profile
        f = sequences[iseq]

        # Gradient direction
        dir = directions[:, idir]

        # Display state of iterations
        @printf "Solving MF with size %d and neig = %d\n" (sum(npoint_cmpts)) length(λ)
        @printf "  Direction %d of %d: g = [%.2f, %.2f, %.2f]\n" idir ndirection dir...
        @printf "  Sequence %d of %d: f = %s\n" iseq nsequence f
        @printf "  Amplitude %d of %d: q = %g, b = %g\n" iamp namplitude q b

        # Gradient direction dependent finite element matrix
        A = sum(dir[i] * Ax[:, :, i] for i = 1:3)

        # Create array for magnetization coefficients
        ν = copy(ν₀)

        if isa(f, PGSE)
            # Constant Bloch-Torrey operator in Laplace eigenfunction basis
            K = L + T + im * q * A
            expmv!(-f.δ, K, ν)
            expmv!(-(f.Δ - f.δ), L + T, ν)
            expmv!(-f.δ, K', ν)
            # edK = exp(-f.δ * K)
            # edL = exp(-(f.Δ - f.δ) * (L + T))
            # ν = edK' * (edL * (edK * ν))
        else
            # Bloch-Torrey operator in Laplace eigenfunction basis for given
            # time profile value
            K(fᵢ) = L + T + im * q * fᵢ * A
            t = LinRange(0, echotime(f), ninterval + 1)
            for i = 1:ninterval
                δᵢ = t[i+1] - t[i]
                fᵢ = (f(t[i+1]) + f(t[i])) / 2
                # ν .= exp(-δᵢ * K(fᵢ)) * ν
                expmv!(-δᵢ, K(fᵢ), ν)
            end
        end

        # Final magnetization coefficients in finite element nodal basis
        mag = ϕ * ν

        # Store results
        for icmpt = 1:ncompartment
            inds = get_inds(icmpt)

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = mag[inds]

            # Integrate magnetization over compartment
            signal[icmpt, iamp, iseq, idir] = sum(M_cmpts[icmpt] * mag[inds])
        end

        itertimes[iamp, iseq, idir] = Base.time() - itertime
    end

    signal_allcmpts[:] = sum(signal, dims = 1)

    totaltime = Base.time() - starttime

    (; magnetization, signal, signal_allcmpts, totaltime, itertimes)
end
