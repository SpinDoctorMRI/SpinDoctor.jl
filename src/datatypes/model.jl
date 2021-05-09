@with_kw struct Model
    name::String
    mesh::FEMesh
    ρ::Vector{Float64}
    D::Vector{Matrix{Float64}}
    T₂::Vector{Float64}
    κ::Vector{Float64}
end
