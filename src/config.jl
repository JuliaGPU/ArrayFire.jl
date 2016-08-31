let
    global af_lib
    if !isdefined(:af_lib)
        lib = is_unix()  ? begin "libaf" end : is_windows() ? "af" : ""
        try
            Libdl.dlopen(lib)
        catch
            error("ArrayFire library not found")
        end
        @eval const af_lib = $lib
    end
end
