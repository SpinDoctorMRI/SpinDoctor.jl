"""
    isconstant(profile)

Return `true` if the time profile `profile` is interval-wise constant.
"""
function isconstant end

isconstant(::AbstractTimeProfile) = false
isconstant(::PGSE) = true
isconstant(::DoublePGSE) = true
