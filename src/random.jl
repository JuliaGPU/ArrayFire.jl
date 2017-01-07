import Base: srand, rand, randn

function rand{T,N}(::Type{AFArray}, ::Type{T}, t::NTuple{N,Int})
    out = RefValue{af_array}(0)
    _error(ccall((:af_randu,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),
                 out,UInt32(N),[t...],af_type(T)))
    AFArray{T,N}(out[])
end
rand{T,N}(::Type{AFArray}, ::Type{T}, t::Vararg{Int,N}) = rand(AFArray, T, t)
rand{N}(::Type{AFArray}, t::NTuple{N,Int}) = rand(AFArray, Float32, t)
rand(::Type{AFArray}, t::Int...) = rand(AFArray, t)

function randn{T,N}(::Type{AFArray}, ::Type{T}, t::NTuple{N,Int})
    out = RefValue{af_array}(0)
    _error(ccall((:af_randn,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),
                 out,UInt32(N),[t...],af_type(T)))
    AFArray{T,N}(out[])
end
randn{T,N}(::Type{AFArray}, ::Type{T}, t::Vararg{Int,N}) = randn(AFArray, T, t)
randn{N}(::Type{AFArray}, t::NTuple{N,Int}) = randn(AFArray, Float32, t)
randn(::Type{AFArray}, t::Int...) = randn(AFArray, t)

srand(::Type{AFArray}, i) = set_seed(i)
