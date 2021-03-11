# Deform domain by bending and twisting.
function deform_domain!(points, deformation)

    height =  maximum(points[3, :]) -  minimum(points[3, :])
    width = maximum(points[1, :]) -  minimum(points[1, :])

    bend = deformation[1]
    twist = deformation[2]

    # Ben points in x-direction as a percentage of height.
    zcenter = (maximum(points[3, :]) + minimum(points[3, :])) / 2
    @. points[1, :] += bend * 30. * width * ((points[3, :] - zcenter) / height)^2.

    # Twist angles by height
    thvec = points[[3], :] / height * twist

    # Twist around center
    center = (maximum(points[1:2, :], dims=2) .+ minimum(points[1:2, :], dims=2)) ./ 2
    points[1:2, :] .-= center

    points[1:2, :] .= [
        cos.(thvec) .* points[[1], :] .- sin.(thvec) .* points[[2], :]
        sin.(thvec) .* points[[1], :] .+ cos.(thvec) .* points[[2], :]
    ]

    # points[1, :] = points[1, :] - rmax
    points[1:2, :] .+= center

end
