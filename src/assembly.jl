# P1 finite element matrices on simplices (segments, triangles, tetrahedra).
#
# The matrices are assembled column by column: column `j` collects the contributions of
# the elements containing node `j`. Each column is written by one thread only, and the
# contributions are summed in the order of the elements, so the result does not depend on
# the number of threads.

"""
    simplex(points, elements, e)

Vertices of element `e` as the columns of a static matrix.
"""
@inline function simplex(
    points::AbstractMatrix{T},
    elements,
    e,
    ::Val{dim},
    ::Val{n},
) where {T,dim,n}
    SMatrix{dim,n,T}(
        ntuple(i -> points[(i-1)%dim+1, elements[(i-1)÷dim+1, e]], Val(dim * n)),
    )
end

# Edge vectors from the first vertex
@inline edgevectors(x::SMatrix{dim,n,T}) where {dim,n,T} = SMatrix{dim,n - 1,T}(
    ntuple(i -> x[(i-1)%dim+1, (i-1)÷dim+2] - x[(i-1)%dim+1, 1], Val(dim * (n - 1))),
)

"""
    simplex_measure(x)

Length, area or volume of the simplex with vertices `x` (columns), possibly embedded in a
higher dimensional space (e.g. a triangle in 3D).
"""
@inline function simplex_measure(x::SMatrix{dim,n}) where {dim,n}
    E = edgevectors(x)
    k = n - 1
    if k == dim
        abs(det(E)) / factorial(k)
    else
        sqrt(abs(det(E' * E))) / factorial(k)
    end
end

"""
    local_mass(x)

P1 mass matrix of the simplex `x`: `|K| (1 + δᵢⱼ) / (n (n + 1))` for `n` vertices.
"""
@inline function local_mass(x::SMatrix{dim,n,T}) where {dim,n,T}
    simplex_measure(x) / (n * (n + 1)) * (ones(SMatrix{n,n,T}) + I)
end

"""
    local_mass(x, w)

P1 mass matrix of the simplex `x` weighted by the P1 function with nodal values `w`:
`∫ w φᵢ φⱼ = |K| (n - 1)! / (n + 2)! (2W + 4wᵢ if i = j, else W + wᵢ + wⱼ)`, where `W =
∑ wₗ`.
"""
@inline function local_mass(x::SMatrix{dim,n,T}, w::SVector{n}) where {dim,n,T}
    W = sum(w)
    c = simplex_measure(x) * factorial(n - 1) / factorial(n + 2)
    SMatrix{n,n,T}(ntuple(Val(n * n)) do k
        i, j = (k - 1) % n + 1, (k - 1) ÷ n + 1
        c * (i == j ? 2W + 4w[i] : W + w[i] + w[j])
    end)
end

"""
    local_stiffness(x, D)

P1 stiffness matrix `∫ ∇φᵢ ⋅ D ∇φⱼ` of the simplex `x` (of full dimension) with constant
diffusion tensor `D`.
"""
@inline function local_stiffness(x::SMatrix{dim,n,T}, D) where {dim,n,T}
    E = edgevectors(x)
    B = inv(E)

    # Gradients of the barycentric coordinates (columns)
    ∇λ = hcat(-sum(B'; dims = 2), B')
    K = abs(det(E)) / factorial(dim) * (∇λ' * D * ∇λ)
    (K + K') / 2
end

"""
    FEPattern(elements, nnode)

Sparsity pattern of the P1 matrices of a mesh with `elements`, together with the list of
the elements containing each node.
"""
struct FEPattern
    nnode::Int
    colptr::Vector{Int}
    rowval::Vector{Int}

    # Elements containing node j: elems[ptr[j]:ptr[j+1]-1], with the local index of j
    ptr::Vector{Int}
    elems::Vector{Int}
    locs::Vector{Int}
end

function FEPattern(elements, nnode)
    nbasis, nelement = size(elements)

    # Node-to-element incidence, in increasing element order
    ptr = zeros(Int, nnode + 1)
    for i ∈ elements
        ptr[i+1] += 1
    end
    ptr[1] = 1
    cumsum!(ptr, ptr)
    next = copy(ptr)
    elems = zeros(Int, nbasis * nelement)
    locs = zeros(Int, nbasis * nelement)
    for e = 1:nelement, a = 1:nbasis
        i = elements[a, e]
        elems[next[i]] = e
        locs[next[i]] = a
        next[i] += 1
    end

    # Rows of column j: the nodes of the elements containing j
    cols = Vector{Vector{Int}}(undef, nnode)
    Threads.@threads for j = 1:nnode
        rows = Int[]
        for k = ptr[j]:(ptr[j+1]-1), a = 1:nbasis
            push!(rows, elements[a, elems[k]])
        end
        cols[j] = unique!(sort!(rows))
    end
    colptr = cumsum([1; length.(cols)])
    rowval = zeros(Int, colptr[end] - 1)
    for j = 1:nnode
        rowval[colptr[j]:(colptr[j+1]-1)] = cols[j]
    end

    FEPattern(nnode, colptr, rowval, ptr, elems, locs)
end

"""
    assemble(localmatrix, pattern, elements, ::Val{nbasis})

Assemble the sparse matrix with the element contributions `localmatrix(e)` (static
`nbasis × nbasis` matrices).
"""
function assemble(
    localmatrix::F,
    pattern::FEPattern,
    elements,
    ::Val{nbasis},
) where {F,nbasis}
    (; nnode, colptr, rowval, ptr, elems, locs) = pattern
    size(elements, 2) > 0 || return spzeros(nnode, nnode)
    T = eltype(localmatrix(1))
    nzval = zeros(T, length(rowval))
    Threads.@threads for j = 1:nnode
        rows = view(rowval, colptr[j]:(colptr[j+1]-1))
        for k = ptr[j]:(ptr[j+1]-1)
            e, b = elems[k], locs[k]
            K = localmatrix(e)
            for a = 1:nbasis
                p = colptr[j] - 1 + searchsortedfirst(rows, elements[a, e])
                nzval[p] += K[a, b]
            end
        end
    end
    SparseMatrixCSC(nnode, nnode, copy(colptr), copy(rowval), nzval)
end

"""
    assemble_mass_matrix(elements, points; elementweights = nothing, nodeweights = nothing[, pattern])

Assemble the P1 mass matrix of the simplices `elements` (columns of point indices) with
`points` (columns of coordinates). The integrand may be weighted by a piecewise constant
function (`elementweights`, one value per element) or by a P1 function (`nodeweights`, one
value per point).

The elements may be of lower dimension than the points, e.g. surface triangles in 3D.
"""
function assemble_mass_matrix(
    elements,
    points;
    elementweights = nothing,
    nodeweights = nothing,
    pattern = FEPattern(elements, size(points, 2)),
)
    isnothing(elementweights) ||
        isnothing(nodeweights) ||
        error("Give either element weights or node weights, not both")
    dim = Val(size(points, 1))
    nbasis = Val(size(elements, 1))
    _mass(points, elements, elementweights, nodeweights, pattern, dim, nbasis)
end

function _mass(points, elements, ew, nw, pattern, dim::Val, nb::Val{n}) where {n}
    if !isnothing(ew)
        localmatrix = e -> ew[e] * local_mass(simplex(points, elements, e, dim, nb))
    elseif !isnothing(nw)
        localmatrix =
            e -> local_mass(
                simplex(points, elements, e, dim, nb),
                SVector{n}(ntuple(a -> nw[elements[a, e]], Val(n))),
            )
    else
        localmatrix = e -> local_mass(simplex(points, elements, e, dim, nb))
    end
    assemble(localmatrix, pattern, elements, nb)
end

"""
    assemble_stiffness_matrix(elements, points, D = I[, pattern])

Assemble the P1 stiffness matrix `∫ ∇φᵢ ⋅ D ∇φⱼ` of the simplices `elements` with `points`,
for a constant diffusion tensor `D`.
"""
function assemble_stiffness_matrix(
    elements,
    points,
    D = I;
    pattern = FEPattern(elements, size(points, 2)),
)
    dim = size(points, 1)
    _stiffness(points, elements, SMatrix{dim,dim}(D * I(dim)), pattern, Val(dim))
end

function _stiffness(points, elements, D, pattern, ::Val{dim}) where {dim}
    nb = Val(dim + 1)
    assemble(
        e -> local_stiffness(simplex(points, elements, e, Val(dim), nb), D),
        pattern,
        elements,
        nb,
    )
end

"""
    facet_integrals(points, facets)

Integrals `∫ φᵢ dS` of the P1 basis functions over the boundary `facets` (simplices of one
dimension lower than the points).
"""
function facet_integrals(points, facets)
    w = zeros(eltype(points), size(points, 2))
    _facet_integrals!(w, points, facets, Val(size(points, 1)), Val(size(facets, 1)))
end

function _facet_integrals!(w, points, facets, dim, ::Val{n}) where {n}
    for f ∈ axes(facets, 2)
        a = simplex_measure(simplex(points, facets, f, dim, Val(n))) / n
        for i = 1:n
            w[facets[i, f]] += a
        end
    end
    w
end

"""
    facet_mass_triplets!(I, J, V, points, facets, rows, cols, c)

Append the entries of `c` times the P1 mass matrix of the `facets` to the triplet lists
`I`, `J`, `V`, with the facet point `facets[a, f]` mapped to row `rows[a, f]` and column
`cols[a, f]`.
"""
function facet_mass_triplets!(I, J, V, points, facets, rows, cols, c)
    _facet_mass_triplets!(
        I,
        J,
        V,
        points,
        facets,
        rows,
        cols,
        c,
        Val(size(points, 1)),
        Val(size(facets, 1)),
    )
end

function _facet_mass_triplets!(
    I,
    J,
    V,
    points,
    facets,
    rows,
    cols,
    c,
    dim,
    ::Val{n},
) where {n}
    for f ∈ axes(facets, 2)
        K = local_mass(simplex(points, facets, f, dim, Val(n)))
        for b = 1:n, a = 1:n
            push!(I, rows[a, f])
            push!(J, cols[b, f])
            push!(V, c * K[a, b])
        end
    end
end
