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

function af_get_numdims!(arr::af_array)
    result = RefValue{UInt32}(0)
    af_error(ccall((:af_get_numdims,af_lib),af_err,(Ptr{UInt32},af_array),result,arr))
    Int(result[])
end

function af_get_type!(arr::af_array)
    _type = RefValue{af_dtype}(0)
    af_error(ccall((:af_get_type,af_lib),af_err,(Ptr{af_dtype},af_array),_type,arr))
    t = _type[]
    if t == f32
        return Float32
    elseif t == c32
        return Complex{Float32}
    elseif t == f64
        return Float64
    elseif t == c64
        return Complex{Float64}
    elseif t == b8
        return Bool
    elseif t == s32
        return Int32
    elseif t == u32
        return UInt32
    elseif t == u8
        return UInt8
    elseif t == s64
        return Int64
    elseif t == u64
        return UInt64
    end
end
