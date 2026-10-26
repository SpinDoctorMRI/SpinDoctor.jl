# Coefficients of the compartments and boundaries of the setup recipes

"""
    coefficients(setup; D, T2, rho, kappa, gamma)

Distribute the coefficients over the compartments and boundaries of `setup`, and return
them as a named tuple `(; D, T2, rho, kappa, gamma)`, from which a [`Model`](@ref) may be
created. The coefficients are given as follows:

- plates and slabs: `D`, `T2` and `rho` are vectors with one value per plate, and
  `kappa = (; interfaces, boundaries)` contains the permeabilities between consecutive
  plates and the relaxivities of the outer boundaries of the plates;
- disks, cylinders and spheres: `D`, `T2` and `rho` are named tuples `(; cell, ecs)`, with
  one value per layer for `cell`, and `kappa = (; cell_interfaces, cell_boundaries,
  cell_ecs, ecs)` contains the permeabilities between consecutive layers, the relaxivities
  of the outer boundaries of the layers (e.g. the ends of extruded cylinders), the
  permeability between the cells and the ECS and the relaxivity of the outer boundary of
  the ECS;
- neurons: `D`, `T2` and `rho` are named tuples `(; neuron, ecs)`, and
  `kappa = (; neuron_ecs, neuron, ecs)`.

The diffusion tensors `D` are matrices (µm²/µs), `T2` times in µs, permeabilities and
relaxivities in µm/µs and `gamma` the gyromagnetic ratio (see [`Model`](@ref)).
"""
function coefficients end

coefficients(setup::ExtrusionSetup; kwargs...) = coefficients(setup.groundsetup; kwargs...)

function coefficients(setup::PlateSetup{T}; D, T2, rho, kappa, gamma) where {T}
    n = length(setup.widths)
    coeffs = (;
        D = [Matrix{T}(D[i]) for i = 1:n],
        T2 = Vector{T}(T2),
        kappa = zeros(T, nboundary(n)),
        rho = Vector{Complex{T}}(rho),
        gamma = T(gamma),
    )
    for i = 1:(n-1)
        coeffs.kappa[interface_index(i, i + 1, n)] = kappa.interfaces[i]
    end
    for i = 1:n
        coeffs.kappa[outer_boundary_index(i, n)] = kappa.boundaries[i]
    end
    coeffs
end

function coefficients(
    setup::Union{DiskSetup{T},SphereSetup{T}};
    D,
    T2,
    rho,
    kappa,
    gamma,
) where {T}
    (; ncell, layersizes, ecs) = setup
    include_ecs = !(ecs isa NoECS)
    nlayer = length(layersizes)
    n = nlayer * ncell + include_ecs

    coeffs = (;
        D = [zeros(T, size(D.cell[1])) for _ = 1:n],
        T2 = zeros(T, n),
        kappa = zeros(T, nboundary(n)),
        rho = zeros(Complex{T}, n),
        gamma = T(gamma),
    )

    # The compartments are numbered layer by layer, then the ECS
    k = 1
    for i = 1:nlayer, j = 1:ncell
        coeffs.rho[k] = rho.cell[i]
        coeffs.D[k] .= D.cell[i]
        coeffs.T2[k] = T2.cell[i]

        if i < nlayer
            # Interface with the next layer of the same cell
            coeffs.kappa[interface_index(k, k + ncell, n)] = kappa.cell_interfaces[i]
        elseif include_ecs
            # Interface between the outer layer and the ECS
            coeffs.kappa[interface_index(k, n, n)] = kappa.cell_ecs
        end

        # Outer boundary: may be empty in 2D, but not in extrusions (top and bottom)
        coeffs.kappa[outer_boundary_index(k, n)] = kappa.cell_boundaries[i]

        k += 1
    end

    if include_ecs
        coeffs.rho[end] = rho.ecs
        coeffs.D[end] .= D.ecs
        coeffs.T2[end] = T2.ecs
        coeffs.kappa[end] = kappa.ecs
    end

    coeffs
end

function coefficients(setup::NeuronSetup{T}; D, T2, rho, kappa, gamma) where {T}
    if setup.ecs isa NoECS
        (;
            D = [Matrix{T}(D.neuron)],
            T2 = T[T2.neuron],
            kappa = T[kappa.neuron],
            rho = Complex{T}[rho.neuron],
            gamma = T(gamma),
        )
    else
        # Compartments: neuron, ECS. Boundaries: interface, outer boundary of the ECS.
        (;
            D = [Matrix{T}(D.neuron), Matrix{T}(D.ecs)],
            T2 = T[T2.neuron, T2.ecs],
            kappa = T[kappa.neuron_ecs, kappa.ecs],
            rho = Complex{T}[rho.neuron, rho.ecs],
            gamma = T(gamma),
        )
    end
end

"""
    analytical_coefficients(setup, coeffs)

Coefficients of the analytical matrix formalism (in SI units) for plates, slabs, disks,
cylinders and spheres with one cell, from the coefficients `coeffs` of
[`coefficients`](@ref). Disks, cylinders and spheres have the mean radius.
"""
analytical_coefficients(s::ExtrusionSetup, coeffs) =
    analytical_coefficients(s.groundsetup, coeffs)

function analytical_coefficients(setup::PlateSetup, coeffs)
    r = cumsum(setup.widths)
    radial_coefficients(setup, r, coeffs)
end

function analytical_coefficients(setup::Union{DiskSetup,SphereSetup}, coeffs)
    r = setup.layersizes * (setup.rmin + setup.rmax) / 2
    radial_coefficients(setup, r, coeffs)
end

# Layers with outer radii `r` (µm), converted from µm and µs to m and s
function radial_coefficients(setup, r, coeffs)
    n = length(r)
    (;
        rho = coeffs.rho,
        r = r * 1e-6,
        D = tr.(coeffs.D) ./ size(coeffs.D[1], 1) .* 1e-6,
        kappa = [
            [coeffs.kappa[interface_index(i, i + 1, n)] for i = 1:(n-1)];
            coeffs.kappa[end]
        ],
        T2 = coeffs.T2 * 1e-6,
        gamma = coeffs.gamma,
        dim = radial_dimension(setup),
    )
end
