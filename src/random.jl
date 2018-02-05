import Base: srand, rand, randn

function rand(::Type{AFArray{T,N}}, t::NTuple{N,Int}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_randu,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),
                 out,UInt32(N),[t...],af_type(T)))
    AFArray{T,N}(out[])
end
rand{T,N}(::Type{AFArray{T}}, t::NTuple{N,Int}) = rand(AFArray{T,N}, t)
rand{N}(::Type{AFArray}, t::NTuple{N,Int}) = rand(AFArray{Float32,N}, t)
rand{T}(::Type{AFArray{T}}, t::Int...) = rand(AFArray{T}, t)
rand(::Type{AFArray}, t::Int...) = rand(AFArray, t)

function randn(::Type{AFArray{T,N}}, t::NTuple{N,Int}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_randn,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),
                 out,UInt32(N),[t...],af_type(T)))
    AFArray{T,N}(out[])
end
randn{T,N}(::Type{AFArray{T}}, t::NTuple{N,Int}) = randn(AFArray{T,N}, t)
randn{N}(::Type{AFArray}, t::NTuple{N,Int}) = randn(AFArray{Float32,N}, t)
randn{T}(::Type{AFArray{T}}, t::Int...) = randn(AFArray{T}, t)
randn(::Type{AFArray}, t::Int...) = randn(AFArray, t)

srand(::Type{AFArray}, i) = set_seed(UInt(abs(i)))
