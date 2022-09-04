"""
    echotime(f::AbstractTimeProfile)

Get echo time ``TE`` of the time profile `f`, which is the end of the last characteristic
interval.

    echotime(grad::AbstractGradient)

Get echo time ``TE`` of `gradient`.
"""
function echotime end

echotime(f::AbstractTimeProfile) = intervals(f)[end]
echotime(g::GeneralGradient) = g.TE
echotime(g::ScalarGradient) = echotime(g.profile)
