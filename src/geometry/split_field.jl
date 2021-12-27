"""
    split_field(mesh, ξ)

Split the field `ξ` into a vector containing a `view` of each compartment. 
"""
function split_field(mesh, ξ)
    ncompartment = length(mesh.points)
    npoint_cmpts = size.(mesh.points, 2)
    inds = [0; cumsum(npoint_cmpts[1:end])]
    [@view(ξ[1+inds[i]:inds[i+1]]) for i = 1:ncompartment]
end
