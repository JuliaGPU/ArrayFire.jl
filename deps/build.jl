global const af_lib = Sys.isunix() ? "libaf" : "af"

using Libdl
try
    dlopen(af_lib)
catch
    @warn "ArrayFire library not loaded:
    either download a binary from the official site http://arrayfire.com/download
    or you can build from source https://github.com/arrayfire/arrayfire"
    rethrow()
end
