let
    global af_lib
    succeeded = false
    if !isdefined(:af_lib)
        lib = is_unix()  ? begin "libaf" end : is_windows() ? "af" : ""
        Libdl.dlopen(lib)
        succeeded = true
        succeeded || error("ArrayFire library not found")
        @eval const af_lib = $lib
    end
end
