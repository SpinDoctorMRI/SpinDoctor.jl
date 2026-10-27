# Biological model: mesh and coefficients

"""
    Model(; mesh, rho, D, T2, kappa, gamma)

Finite element discretized biological model:

- `mesh`: finite element mesh split into compartments ([`FEMesh`](@ref));
- `rho`: initial spin density of each compartment;
- `D`: diffusion tensor of each compartment (µm²/µs): a matrix, or a vector with one matrix
  per element of the compartment (e.g. from diffusion tensor imaging);
- `T2`: T2 relaxation time of each compartment (µs);
- `kappa`: permeability of each boundary (µm/µs), for interfaces between compartments and
  for the outer boundaries (as a boundary relaxivity);
- `gamma`: gyromagnetic ratio (rad / (µs µm T/m), `2.67513e-4` for water protons).

The coefficients for a setup recipe may be created with [`coefficients`](@ref).
"""
Base.@kwdef struct Model{T,dim}
    mesh::FEMesh{T,dim}
    rho::Vector{Complex{T}}
    D::Vector{Union{Matrix{T},Vector{Matrix{T}}}}
    T2::Vector{T}
    kappa::Vector{T}
    gamma::T
end

function Model(mesh::FEMesh{T,dim}, rho, D, T2, kappa, gamma) where {T,dim}
    D = [d isa AbstractVector ? Matrix{T}.(d) : Matrix{T}(d * I(dim)) for d in D]
    Model{T,dim}(mesh, rho, D, T2, kappa, gamma)
end

"""
    mean_diffusion_tensors(model)

Volume average of the diffusion tensor in each compartment.
"""
mean_diffusion_tensors(model) =
    map(model.D, model.mesh.points, model.mesh.elements) do D, points, elements
        D isa AbstractVector || return D
        volumes, _ = element_geometry(points, elements)
        sum(volumes .* D) / sum(volumes)
    end

"""
    initial_conditions(model)

Initial magnetization (the spin density) at all the points of the mesh.
"""
initial_conditions(model) =
    mapreduce((rho, p) -> fill(rho, size(p, 2)), vcat, model.rho, model.mesh.points)

"""
    compartment_weights(model; TE = 0)

Relative contributions of the compartments to the signal at `b = 0`: the initial spin
density times the volume, attenuated by T2 relaxation at the echo time `TE`, and normalized
to sum to one.

The ADC of the whole sample in the low-`b` limit is the weighted average of the compartment
ADCs without exchange, e.g. `compartment_weights(model; TE)' * adcs` for compartment ADCs
from [`solve_hadc`](@ref).
"""
function compartment_weights(model; TE = 0)
    volumes = compartment_volumes(model.mesh)
    w = abs.(model.rho) .* volumes .* exp.(-TE ./ model.T2)
    w / sum(w)
end
