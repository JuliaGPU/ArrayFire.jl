let
    global af_lib
    succeeded = false
    if !isdefined(:af_lib)
        lib = is_unix()  ? begin "libaf" end : is_windows() ? "af" : ""
        try
            Libdl.dlopen(lib)
            succeeded = true
        end
        @eval const af_lib = $lib
    end
    succeeded || error("ArrayFire library not found")
end
