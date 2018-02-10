import Base: RefValue, @pure, show, clamp, find
import Base: cumsum, cumprod, cummin, cummax, chol, abs2

export constant, get_last_error, err_to_string, sort_index, fir, iir
export mean_weighted, var_weighted, set_array_indexer, set_seq_param_indexer
export afeval, iota, sortbykey

const af_threshold = Ref(4*1024*1024*1024)
const af_gc_count = Ref(0)
const af_gc_frequency = Ref(256)

function _afgc()
    af_gc_count[] += 1
    if mod(af_gc_count[], af_gc_frequency[]) == 0
        afgc()
    end
end

function afgc(threshold = 0)
    if threshold != 0
        af_threshold[] = threshold
    end
    alloc_bytes, alloc_buffers, lock_bytes, lock_buffers =  device_mem_info()
    if alloc_bytes > af_threshold[]
        gc()
        device_gc()
    end
    nothing
end

function release_array(arr::AFArray)
    ccall((:af_release_array,af_lib),af_err,(af_array,),arr.arr)
end

export afgc

toa(a) = issparse(a) ?  SparseMatrixCSC(a) : Array(a)
show(c::IOContext, a::AFArray) = (print(c, "AFArray: "); show(c, toa(a)))
show(io::IO, a::AFArray) = print(io, "AFArray: ", toa(a))
function show(io::IO, m::MIME"text/plain", a::AFArray)
    print(io, "AFArray: "); show(io, m, toa(a))
end

global const af_lib = is_unix() ? "libaf" : "af"
global const bcast = Ref{Bool}(false)

function __init__()
    Libdl.dlopen(af_lib)

    backend_envvar_name = "JULIA_ARRAYFIRE_BACKEND"
    if haskey(ENV, backend_envvar_name)
        backend_str = lowercase(ENV[backend_envvar_name])
        backends = Dict(
            "cpu" => AF_BACKEND_CPU,
            "cuda" => AF_BACKEND_CUDA,
            "opencl" => AF_BACKEND_OPENCL
        )
        if haskey(backends, backend_str)
            set_backend(backends[backend_str])
        else
            error("Unknown arrayfire backend \"$backend_str\".")
        end
    end

    afinit()
    if !haskey(ENV, "JULIA_ARRAYFIRE_SILENT")
        afinfo()
    end
    set_seed(rand(RandomDevice(), UInt64))
    nothing
end

function _error(err::af_err)
    if err == 0
        _afgc()
    else
        if err == 101
            error("GPU is out of memory, to avoid this in the future you can:
  @afgc function f(input)   # free all temporary variables inside the function scope
  @afgc a = b + c           # free all temporary arrays inside the assignment scope
  @afgc a .= b + c          # replace with a new value, free the old and all temps
  afgc(threshold)           # garbage collect after GPU memory usage reaches threshold
  finalize(array)           # manually free GPU memory")
        else
            str = err_to_string(err)
            str2 = get_last_error()
            error("ArrayFire Error ($err) : $str\n$str2")
        end
    end
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

function af_jltype(i::af_dtype)::Type
    if i == f32
        return Float32
    elseif i == c32
        return Complex{Float32}
    elseif i == f64
        return Float64
    elseif i == c64
        return Complex{Float64}
    elseif i == b8
        return Bool
    elseif i == s32
        return Int32
    elseif i == u32
        return UInt32
    elseif i == u8
        return UInt8
    elseif i == s64
        return Int64
    elseif i == u64
        return UInt64
    else
        error("Unknown type: $i")
    end
end

function get_numdims(arr::af_array)
    result = RefValue{UInt32}(0)
    _error(ccall((:af_get_numdims,af_lib),af_err,
                 (Ptr{UInt32},af_array),result,arr))
    Int(result[])
end

function get_type(arr::af_array)
    _type = RefValue{af_dtype}(0)
    _error(ccall((:af_get_type,af_lib),af_err,
                 (Ptr{af_dtype},af_array),_type,arr))
    af_jltype(_type[])
end

function check_type_numdims{T,N}(arr::AFArray{T,N})
    @assert get_type(arr) == T "type mismatch: $(get_type(arr)) != $T"
    @assert get_numdims(arr) == N "dims mismatch: $(get_numdims(arr)) != $N"
end

function convert_array{T,N}(data::AbstractArray{T,N})
    arr = RefValue{af_array}(0)
    sz = size(data)
    _error(ccall((:af_create_array,af_lib),af_err,
                 (Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),
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

function convert_array_to_sparse(a::SparseMatrixCSC)
    sz = size(a)
    at = transpose(a)
    colptr = AFArray(Vector{Int32}(at.colptr-1))
    rowval = AFArray(Vector{Int32}(at.rowval-1))
    create_sparse_array(sz[1], sz[2], AFArray(at.nzval), colptr, rowval, AF_STORAGE_CSR)
end

function convert_array_to_sparse(a::AFArray)
    @assert issparse(a) "AFArray is not sparse"
    sz = size(a)
    @assert length(sz) == 2 "AFArray is not a matrix"
    nzval, colptr, rowval, d = sparse_get_info(a)
    if d == AF_STORAGE_CSR
        transpose(SparseMatrixCSC(sz[2], sz[1],
                                  Array(colptr+1), Array(rowval+1), Array(nzval)))
    else
        convert_array_to_sparse(sparse_convert_to(a, AF_STORAGE_CSR))
    end
end

function recast_array{T1,N,T2}(::Type{AFArray{T1}},_in::AFArray{T2,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cast,af_lib),af_err,
                 (Ptr{af_array},af_array,af_dtype),out,_in.arr,af_type(T1)))
    AFArray{T1,N}(out[])
end

AFArray!(arr::af_array) = AFArray{get_type(arr), get_numdims(arr)}(arr)

function constant{T<:Real,N}(val::T,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant,af_lib),af_err,
                 (Ptr{af_array},Cdouble,UInt32,Ptr{dim_t},af_dtype),
                 arr,Cdouble(val),UInt32(N),[sz...],af_type(T)))
    AFArray{T,N}(arr[])
end

function constant{N}(val::Complex{Bool},sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_complex,af_lib),af_err,
                 (Ptr{af_array},Cdouble,Cdouble,UInt32,Ptr{dim_t},af_dtype),
                 arr,Cdouble(real(val)),Cdouble(imag(val)),UInt32(N),[sz...],c32))
    AFArray{Complex{Float32},N}(arr[])
end

function constant{T<:Complex,N}(val::T,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_complex,af_lib),af_err,
                 (Ptr{af_array},Cdouble,Cdouble,UInt32,Ptr{dim_t},af_dtype),
                 arr,Cdouble(real(val)),Cdouble(imag(val)),UInt32(N),[sz...],af_type(T)))
    AFArray{T,N}(arr[])
end

function constant{N}(val::Int,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_long,af_lib),af_err,
                 (Ptr{af_array},intl,UInt32,Ptr{dim_t}),
                 arr,val,UInt32(N),[sz...]))
    AFArray{Int,N}(arr[])
end

function constant{N}(val::UInt,sz::NTuple{N,Int})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_ulong,af_lib),af_err,
                 (Ptr{af_array},uintl,UInt32,Ptr{dim_t}),
                 arr,val,UInt32(N),[sz...]))
    AFArray{UInt,N}(arr[])
end

function ifelse{T1,N1,T2,N2}(cond::AFArray{Bool},a::AFArray{T1,N1},b::AFArray{T2,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_select,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,af_array),
                 out,cond.arr,a.arr,b.arr))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function ifelse{T1,N1,T2<:Real}(cond::AFArray{Bool},a::AFArray{T1,N1},b::T2)
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_r,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,Cdouble),
                 out,cond.arr,a.arr,Cdouble(b)))
    AFArray{typed(T1,T2),N1}(out[])
end

function ifelse{T1,T2,N2}(cond::AFArray{Bool},a::T1,b::AFArray{T2,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_l,af_lib),af_err,
                 (Ptr{af_array},af_array,Cdouble,af_array),
                 out,cond.arr,Cdouble(a),b.arr))
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

function cat{T,N1,N2}(dim::Integer,first::AFArray{T,N1},second::AFArray{T,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_join,af_lib),af_err,
                 (Ptr{af_array},Cint,af_array,af_array),
                 out,Cint(dim - 1),first.arr,second.arr))
    AFArray{T,batched(N1,N2)}(out[])
end

hcat(first::AFArray, second::AFArray) = cat(2, first, second)
vcat(first::AFArray, second::AFArray) = cat(1, first, second)

function conv{T,N}(signal::AFArray{T,N}, filter::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve1,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),
                 out,signal.arr,filter.arr,AF_CONV_EXPAND,AF_CONV_AUTO))
    AFArray{T,N}(out[])
end

function conv_fft{T,N}(signal::AFArray{T,N}, filter::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),
                 out,signal.arr,filter.arr,AF_CONV_EXPAND))
    AFArray{T,N}(out[])
end

norm(a::AFArray{Float32,2})::Float32 = allowslow(AFArray) do; svd(a)[2][1]; end
norm(a::AFArray{Float64,2}) = allowslow(AFArray) do; svd(a)[2][1]; end
norm(a::AFArray) = norm(a, AF_NORM_EUCLID, 0, 0)
vecnorm(a::AFArray) = norm(a, AF_NORM_EUCLID, 0, 0)

function svd{T}(_in::AFArray{T,2})
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    _error(ccall((:af_svd,af_lib),af_err,
                 (Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),
                 u,s,vt,_in.arr))
    (AFArray{T,2}(u[]),AFArray{T,1}(s[]),AFArray{T,2}(vt[]))
end

function sort{T,N}(_in::AFArray{T,N},dim::Integer=1,isAscending::Bool=true)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sort,af_lib),af_err,
                 (Ptr{af_array},af_array,UInt32,Bool),
                 out,_in.arr,UInt32(dim - 1),isAscending))
    AFArray{T,N}(out[])
end

function sort_index{T,N}(_in::AFArray{T,N},dim::Integer=1,isAscending::Bool=true)
    out = RefValue{af_array}(0)
    indices = RefValue{af_array}(0)
    _error(ccall((:af_sort_index,af_lib),af_err,
                 (Ptr{af_array},Ptr{af_array},af_array,UInt32,Bool),
                 out,indices,_in.arr,UInt32(dim - 1),isAscending))
    (AFArray{T,N}(out[]),AFArray{UInt32,N}(indices[])+UInt32(1))
end

function sortperm{T,N}(a::AFArray{T,N}, dim::Integer=1,isAscending::Bool=true)
    sort_index(a,dim,isAscending)[2]
end

function mean{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean,af_lib),af_err,(Ptr{af_array},af_array,dim_t),
                 out,_in.arr,dim-1))
    AFArray{T,N}(out[])
end

function mean_weighted{T,N}(_in::AFArray{T,N},weights::AFArray,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean_weighted,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,dim_t),
                 out,_in.arr,weights.arr,dim-1))
    AFArray{T,N}(out[])
end

function var{T,N}(_in::AFArray{T,N},isbiased::Bool,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_var,af_lib),af_err,
                 (Ptr{af_array},af_array,Bool,dim_t),
                 out,_in.arr,isbiased,dim-1))
    AFArray{T,N}(out[])
end

function var_weighted{T,N}(_in::AFArray{T,N},weights::AFArray,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_var_weighted,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,dim_t),
                 out,_in.arr,weights.arr,dim-1))
    AFArray{T,N}(out[])
end

function stdev{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_stdev,af_lib),af_err,
                 (Ptr{af_array},af_array,dim_t),out,_in.arr,dim-1))
    AFArray{T,N}(out[])
end

function median{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_median,af_lib),af_err,
                 (Ptr{af_array},af_array,dim_t),
                 out,_in.arr,dim-1))
    AFArray{T,N}(out[])
end

function set_seq_param_indexer(_begin::Real,_end::Real,step::Real,dim::dim_t,is_batch::Bool)
    indexer = RefValue{af_index_t}(0)
    _error(ccall((:af_set_seq_param_indexer,af_lib),af_err,
                 (Ptr{af_index_t},Cdouble,Cdouble,Cdouble,dim_t,Bool),
                 indexer,Cdouble(_begin),Cdouble(_end),Cdouble(step),dim-1,is_batch))
    indexer[]
end

function clamp{T1,N1,T2,N2}(_in::AFArray{T1,N1},lo::AFArray{T2,N2},hi::AFArray{T2,N2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_clamp,af_lib),af_err,
                 (Ptr{af_array},af_array,af_array,af_array,Bool),
                 out,_in.arr,lo.arr,hi.arr,bcast[]))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function find{T,N}(_in::AFArray{T,N})
    idx = RefValue{af_array}(0)
    _error(ccall((:af_where,af_lib),af_err,(Ptr{af_array},af_array),idx,_in.arr))
    out = AFArray{UInt32,1}(idx[])
    if length(out) > 0
        out = out + UInt32(1)
    end
    return out
end

cumsum(a::AFArray, dim::Int=1) = scan(a, dim, AF_BINARY_ADD, true)
cumprod(a::AFArray, dim::Int=1) = scan(a, dim, AF_BINARY_MUL, true)
cummin(a::AFArray, dim::Int=1) = scan(a, dim, AF_BINARY_MIN, true)
cummax(a::AFArray, dim::Int=1) = scan(a, dim, AF_BINARY_MAX, true)

function sync(a::AFArray)
    afeval(a)
    sync(get_device_id(a))
    a
end

function afeval(a::AFArray)
    _error(ccall((:af_eval,af_lib),af_err,(af_array,),a.arr))
    a
end

function chol{T,N}(_in::AFArray{T,N},is_upper::Bool=false)
    out = RefValue{af_array}(0)
    info = RefValue{Cint}(0)
    _error(ccall((:af_cholesky,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,Bool),out,info,_in.arr,is_upper))
    (AFArray{T,N}(out[]),info[])
end

abs2{T<:Real}(a::AFArray{T}) = a.*a
abs2{T<:Complex}(a::AFArray{T}) = (r = real(a); i = imag(a); r.*r+i.*i)

function complex{T1,N1,T2,N2}(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2})
    batch = bcast[]
    out = RefValue{af_array}(0)
    _error(ccall((:af_cplx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Complex{typed(T1,T2)},batched(N1,N2)}(out[])
end

function fir{T,N}(b::AFArray,x::AFArray{T,N})
    y = RefValue{af_array}(0)
    _error(ccall((:af_fir,af_lib),af_err,(Ptr{af_array},af_array,af_array),y,b.arr,x.arr))
    AFArray{T,N}(y[])
end

function iir{T,N}(b::AFArray{T},a::AFArray{T},x::AFArray{T,N})
    y = RefValue{af_array}(0)
    _error(ccall((:af_iir,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),y,b.arr,a.arr,x.arr))
    AFArray{T,N}(y[])
end

function iota{T,N}(dims::NTuple{N,Int}, typ::Type{T} = Int32)
    out = RefValue{af_array}(0)
    _error(ccall((:af_iota,af_lib), af_err,
                 (Ptr{af_array},UInt32,Ptr{dim_t},UInt32,Ptr{dim_t},af_dtype),
                 out,UInt32(N),[dims...],UInt32(1),[1],af_type(T)))
    AFArray{T,N}(out[])+T(1)
end

function sortbykey{T1,T2,N}(keys::AFArray{T1,N},values::AFArray{T2,N},dim::Integer,isAscending::Bool=true)
    out_keys = RefValue{af_array}(0)
    out_values = RefValue{af_array}(0)
    _error(ccall((:af_sort_by_key,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,UInt32,Bool),out_keys,out_values,keys.arr,values.arr,UInt32(dim - 1),isAscending))
    (AFArray{T1,N}(out_keys[]),AFArray{T2,N}(out_values[]))
end
