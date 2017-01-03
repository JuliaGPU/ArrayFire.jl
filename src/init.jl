global const af_lib = is_unix() ? "libaf" : "af"

function __init__()
    Libdl.dlopen(af_lib)
    afinit()
    nothing
end
