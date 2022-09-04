"""
    compute_lengths(edges, nodes)

Compute edge lengths.
"""
function compute_lengths(edges, nodes)
    v = nodes[:, edges[2, :]] - nodes[:, edges[1, :]]
    norm.(eachcol(v))
end
