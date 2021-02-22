"""
    M, S_cmpts, S_total = solve_mf(mesh, domain, experiment, lap_eig, directions)

Solve for magnetization using Matrix Formalism.
"""
function solve_mf(mesh, domain, experiment, lap_eig, directions)

    # Extract parameters
    @unpack ncompartment = mesh
    @unpack ρ, κ = domain
    @unpack ndirection, sequences, values, values_type, mf = experiment
    @unpack ninterval = mf

    # Laplace eigenmodes
    λ, ϕ, moments = lap_eig.values, lap_eig.funcs, lap_eig.moments

    namplitude = length(values)
    nsequence = length(sequences)
    npoint_cmpts = size.(mesh.points, 2)

    # Compartment index ranges
    inds_start = cumsum([1; npoint_cmpts[1:end-1]])
    inds_stop = cumsum(npoint_cmpts[1:end])
    get_inds(icmpt) = inds_start[icmpt]:inds_stop[icmpt]

    # Assemble mass matrices compartment-wise (for spatial integration)
    M_cmpts = []
    for icmpt = 1:ncompartment
        # Finite elements
        points = mesh.points[icmpt];
        elements = mesh.elements[icmpt];
        volumes = get_mesh_volumes(points, elements);

        # Assemble mass matrix
        push!(M_cmpts, assemble_mass_matrix(elements', volumes))
    end

    # Global mass matrix
    M = blockdiag(M_cmpts...)

    # Create initial conditions (enforce complex values)
    ρ = vcat(fill.(complex(ρ), npoint_cmpts)...)

    # Project initial spin density onto Laplace eigenfunction basis
    T0 = ϕ' * (M * ρ)

    # Q-values and b-values
    if values_type == 'q'
        qvalues = repeat(values, 1, nsequence)
        bvalues = values.^2 .* bvalue_no_q.(sequences')
    else
        bvalues = repeat(values, 1, nsequence)
        qvalues = .√(values ./ bvalue_no_q.(sequences'))
    end

    # Allocate arrays
    signal = zeros(ComplexF64, ncompartment, namplitude, nsequence, ndirection)
    signal_allcmpts = zeros(ComplexF64, namplitude, nsequence, ndirection)
    magnetization = fill(ComplexF64[], ncompartment, namplitude, nsequence, ndirection)

    # Laplace operator in Laplace eigenfunction basis
    L = diagm(λ)

    # Iterate over gradient amplitudes, time profiles and directions
    for (iamp, iseq, idir) ∈ Iterators.product(1:namplitude, 1:nsequence, 1:ndirection)

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
        A = sum(dir[i] * moments[:, :, i] for i = 1:3)

        # Create array for magnetization coefficients
        T = copy(T0)

        if typeof(f) == PGSE
            # Constant Bloch-Torrey operator in Laplace eigenfunction basis
            K = L + im * q * A
            T = expmv!(-f.δ, K, T)
            @. T *= exp(-(f.Δ - f.δ)λ)
            T = expmv!(-f.δ, K', T)
            # edK = exp(-f.δ * K)
            # edL = @. exp(-(f.Δ - f.δ)λ)
            # T = edK' * (edL .* (edK * T))
        else
            # Bloch-Torrey operator in Laplace eigenfunction basis for given
            # time profile value
            display("toto")
            K(fi) = L + im * q * fi * A
            t = LinRange(0, echotime(f), ninterval + 1)
            for i = 1:ninterval
                δi = t[i+1] - t[i]
                fi = (f(t[i+1]) + f(t[i])) / 2
                # T .= exp(-δi * K(fi)) * T
                expmv!(-δi, K(fi), T)
            end
        end

        # Final magnetization coefficients in finite element nodal basis
        mag = ϕ * T;

        # Store results
        for icmpt = 1:ncompartment
            inds = get_inds(icmpt)

            # Store magnetization in compartment
            magnetization[icmpt, iamp, iseq, idir] = mag[inds]

            # Integrate magnetization over compartment
            signal[icmpt, iamp, iseq, idir] = sum(M_cmpts[icmpt] * mag[inds])
        end
    end

    signal_allcmpts[:] = sum(signal; dims=1)

    (; magnetization, signal, signal_allcmpts)
end
