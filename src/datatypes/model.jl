@with_kw struct Model{T}
    name::String
    mesh::FEMesh
    ρ::Vector{T}
    D::Vector{SMatrix{3,3,T,9}}
    T₂::Vector{T}
    κ::Vector{T}
end
