"""
function add_checksum(savepath::String,setup::AbstractSetup)
"""

function add_checksum(path::String,setup::AbstractSetup)
    dictsetup = Dict(key => getfield(setup, key) for key in propertynames(setup))
    io = IOBuffer();
    bson(io,dictsetup)
    checksum = crc32c(io.data)
    (; path, checksum)
end