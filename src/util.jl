import Base: RefValue, @pure, display, show

export constant, select, get_last_error, err_to_string

function afgc(threshold = 6e9)
    alloc_bytes = RefValue{Csize_t}(0)
    alloc_buffers = RefValue{Csize_t}(0)
    lock_bytes = RefValue{Csize_t}(0)
    lock_buffers = RefValue{Csize_t}(0)
    err = ccall((:af_device_mem_info,af_lib),af_err,(Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t}),
                alloc_bytes,alloc_buffers,lock_bytes,lock_buffers)
    if err == 0 && alloc_bytes[] > threshold
        gc()
        err = ccall((:af_device_gc,af_lib),af_err,())
    end
    return err
end

function release_array(arr::AFArray)
    ccall((:af_release_array,af_lib),af_err,(af_array,),arr.arr)
end

export afgc

display(a::AFArray) = (print("AFArray: "); display(Array(a)))
show(c::IOContext, a::AFArray) = (print(c, "AFArray: "); show(c, Array(a)))

global const af_lib = is_unix() ? "libaf" : "af"
global const bcast = Ref{Bool}(false)

function __init__()
    Libdl.dlopen(af_lib)
    afinit()
    afinfo()
    nothing
end

function _error(err::af_err)
    if err == 0
        err = afgc()
        if err == 0
            return
        end
    end
    str = err_to_string(err)
    str2 = get_last_error()
    throw(ErrorException("ArrayFire Error ($err) : $str\n$str2"))
end

@pure batched(n1, n2) = max(n1, n2)

function typed{T1,T2}(::Type{T1},::Type{T2})
    if T1 == T2
        return T1
    elseif T1 == Complex{Float64} || T2 == Complex{Float64}
        return Complex{Float64}
    elseif T1 == Complex{Float32} || T2 == Complex{Float32}
        (T1 == Float64 || T2 == Float64) && return Complex{Float64}
        return Complex{Float32}
    elseif T1 == Float64 || T2 == Float64
        return Float64
    elseif T1 == Float32 || T2 == Float32
        return Float32
    elseif T1 == UInt64 || T2 == UInt64
        return UInt64
    elseif T1 == Int64 || T2 == Int64
        return Int64
    elseif T1 == UInt32 || T2 == UInt32
        return UInt32
    elseif T1 == Int32 || T2 == Int32
        return Int32
    elseif T1 == UInt8 || T2 == UInt8
        return UInt8
    elseif T1 == Bool || T2 == Bool
        return Bool
    else
        return Float32
    end
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

function check_type_numdims{T,N}(arr::AFArray{T,N})
    @assert get_type(arr) == T "type mismatch: $(get_type(arr)) != $T"
    @assert get_numdims(arr) == N "dims mismatch: $(get_numdims(arr)) != $N"
end

function convert_array{T,N}(data::AbstractArray{T,N})
    arr = RefValue{af_array}(0)
    sz = size(data)
    _error(ccall((:af_create_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),
                   arr,data,UInt32(length(sz)),[sz...],af_type(T)))
    AFArray{T,N}(arr[])
end

function convert_array{T,N}(a::AFArray{T,N})
    if issparse(a)
        a = full(a)
    end
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

function select{T1,N1,T2,N2}(cond::AFArray{Bool},a::AFArray{T1,N1},b::AFArray{T2,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_select,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),out,cond.arr,a.arr,b.arr))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function select{T1,N1,T2<:Real}(cond::AFArray{Bool},a::AFArray{T1,N1},b::T2)
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_r,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cdouble),out,cond.arr,a.arr,Cdouble(b)))
    AFArray{typed(T1,T2),N1}(out[])
end

function select{T1,T2,N2}(cond::AFArray{Bool},a::T1,b::AFArray{T2,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_l,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,af_array),out,cond.arr,Cdouble(a),b.arr))
    AFArray{typed(T1,T2),N2}(out[])
end

function err_to_string(err::af_err)
    unsafe_string(ccall((:af_err_to_string,af_lib),Cstring,(af_err,),err))
end

function get_last_error()
    msg = RefValue{Cstring}()
    len = RefValue{dim_t}(0)
    ccall((:af_get_last_error,af_lib),Void,(Ptr{Cstring},Ptr{dim_t}),msg,len)
    unsafe_string(msg[])
end
