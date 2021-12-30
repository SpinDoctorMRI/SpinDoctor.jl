"""
    compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)

Compute effective diffusion tensors using the matrix formalism.
"""
function compute_mf_diffusion_tensor(mesh, M, lap_eig, gradient::ScalarGradient)
    volumes = get_cmpt_volumes(mesh)
    p = reduce(hcat, mesh.points)
    a = p * M * lap_eig.funcs
    j = j_integral.(gradient.profile, lap_eig.values)
    D = a * Diagonal(j) * a' / sum(volumes)
end
