import Base.RefValue

global const af_lib = is_unix() ? "libaf" : "af"
global const bcast = Ref{Bool}(false)

function __init__()
    Libdl.dlopen(af_lib)
    afinit()
    afinfo()
    nothing
end

function _error(err::af_err)
    if err == 0; return; end
    str = err_to_string(err)
    str2 = get_last_error()[1]
    throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))\n$(unsafe_string(str2))"))
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

af_jltype(::Val{f32}) = Float32
af_jltype(::Val{c32}) = Complex{Float32}
af_jltype(::Val{f64}) = Float64
af_jltype(::Val{c64}) = Complex{Float64}
af_jltype(::Val{b8})  = Bool
af_jltype(::Val{s32}) = Int32
af_jltype(::Val{u32}) = UInt32
af_jltype(::Val{u8})  = UInt8
af_jltype(::Val{s64}) = Int64
af_jltype(::Val{u64}) = UInt64

function get_numdims(arr::af_array)
    result = RefValue{UInt32}(0)
    _error(ccall((:af_get_numdims,af_lib),af_err,(Ptr{UInt32},af_array),result,arr))
    Int(result[])
end

function get_type(arr::af_array)
    _type = RefValue{af_dtype}(0)
    _error(ccall((:af_get_type,af_lib),af_err,(Ptr{af_dtype},af_array),_type,arr))
    af_jltype(Val{_type[]}())
end

function convert_array{T,N}(data::AbstractArray{T,N})
    arr = RefValue{af_array}(0)
    sz = size(data)
    _error(ccall((:af_create_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),
                   arr,data,UInt32(length(sz)),[sz...],af_type(T)))
    AFArray{T,N}(arr[])
end

function convert_array{T,N}(a::AFArray{T,N})
    ret = Array{T,N}(size(a))
    get_data_ptr(ret, a)
    ret
end

function recast_array{T1,N,T2}(::Type{AFArray{T1}},_in::AFArray{T2,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cast,af_lib),af_err,(Ptr{af_array},af_array,af_dtype),out,_in.arr,af_type(T1)))
    AFArray{T1,N}(out[])
end

AFArray!(arr::af_array) = AFArray{get_type(arr), get_numdims(arr)}(arr)

function constant{T<:Real,N}(val::T,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant,af_lib),af_err,(Ptr{af_array},Cdouble,UInt32,Ptr{dim_t},af_dtype),
                 arr,Cdouble(val),UInt32(N),[sz...],af_type(T)))
    AFArray{T,N}(arr[])
end

function constant{T<:Complex,N}(val::T,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_complex,af_lib),af_err,(Ptr{af_array},Cdouble,Cdouble,UInt32,Ptr{dim_t},af_dtype),
                 arr,Cdouble(real(val)),Cdouble(imag(val)),UInt32(N),[sz...],af_type(_type)))
    AFArray{T,N}(arr[])
end

function constant{N}(val::Int,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_long,af_lib),af_err,(Ptr{af_array},intl,UInt32,Ptr{dim_t}),
                 arr,val,UInt32(N),[sz...]))
    AFArray{Int,N}(arr[])
end

function constant{N}(val::UInt,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_ulong,af_lib),af_err,(Ptr{af_array},uintl,UInt32,Ptr{dim_t}),
                 arr,val,UInt32(N),[sz...]))
    AFArray{UInt,N}(arr[])
end

export constant
