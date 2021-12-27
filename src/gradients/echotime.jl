"""
    echotime(f)

Get echo time `TE` of the time profile `f`, which is the end of the last characteristic
interval.
"""
echotime(f::TimeProfile) = intervals(f)[end]

"""
    echotime(gradient)

Get echo time `TE` of `gradient`.
"""
echotime(gradient::GeneralGradient) = gradient.TE
echotime(gradient::ScalarGradient) = echotime(gradient.profile)
