let
    global af_lib
    succeeded = false
    if !isdefined(:af_lib)
        lib = @unix ? begin "libaf" end : @windows ? "af" : ""
        Libdl.dlopen(lib)
        succeeded = true
        succeeded || error("ArrayFire library not found")
        @eval const af_lib = $lib
    end
end
