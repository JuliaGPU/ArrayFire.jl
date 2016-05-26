let
    global af_lib
    succeeded = false
    if !isdefined(:af_lib)
        @unix ? begin lib = "libaf" end : nothing
        Libdl.dlopen(lib)
        succeeded = true
        succeeded || error("Gunrock library not found")
        @eval const af_lib = $lib
    end
end
