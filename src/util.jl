function af_error(err::af_err)
    if err == 0; return; end
    str = af_err_to_string(err)
    throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))"))
end

af_type(::Type{Float32})          = f32
af_type(::Type{Complex{Float32}}) = c32
af_type(::Type{Float64})          = f64
af_type(::Type{Complex{Float64}}) = c64
af_type(::Type{Bool})             = b8
af_type(::Type{Int32})            = s32
af_type(::Type{UInt32})           = u32
af_type(::Type{UInt8})            = u8
af_type(::Type{Int64})            = s64
af_type(::Type{UInt64})           = u64
