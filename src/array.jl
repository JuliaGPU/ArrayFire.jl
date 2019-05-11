mutable struct AFArray{T,N} <: AbstractArray{T,N}
    arr::af_array
    function AFArray{T,N}(arr::af_array) where {T,N}
        @assert get_type(arr) == T
        a = new{T,N}(arr)
        finalizer(release_array, a)
    end
end

AFVector{T} = AFArray{T,1}
AFMatrix{T} = AFArray{T,2}
AFVolume{T} = AFArray{T,3}
AFTensor{T} = AFArray{T,4}

export AFArray, AFVector, AFMatrix, AFVolume, AFTensor

import Base: Array, copy, deepcopy_internal, complex, conj
import SparseArrays: SparseMatrixCSC, issparse, sparse

export sparse, select, ctranspose, atan2, full, gradient

sparse(a::AFArray{T,N}) where {T,N} = create_sparse_array_from_dense(a, AF_STORAGE_CSR)

(::Type{AFArray{T1}})(a::AFArray{T2,N}) where {T1,T2,N} = recast_array(AFArray{T1}, a)
(::Type{AFArray{T1,N}})(a::AFArray{T2,N}) where {T1,T2,N} = recast_array(AFArray{T1}, a)

(::Type{Array{T,N}})(a::AFArray{T,N}) where {T,N} = convert_array(a)
(::Type{AFArray{T,N}})(a::Array{T,N}) where {T,N}  = convert_array(a)

(::Type{Array{T}})(a::AFArray{T,N}) where {T,N}  = convert_array(a)
(::Type{AFArray{T}})(a::Array{T,N}) where {T,N}  = convert_array(a)

(::Type{Array})(a::AFArray{T,N}) where {T,N}  = convert_array(a)
(::Type{AFArray})(a::Array{T,N}) where {T,N}  = convert_array(a)

(::Type{AFArray{T,N}})(a::SubArray{T,N}) where {T,N} = convert_array(Array(a))
(::Type{AFArray{T}})(a::SubArray{T,N}) where {T,N} = convert_array(Array(a))
(::Type{AFArray})(a::SubArray{T,N}) where {T,N} = convert_array(Array(a))

(::Type{AFArray})(a::SparseMatrixCSC) = convert_array_to_sparse(a)

(::Type{SparseMatrixCSC{T}})(a::AFArray{T}) where T = convert_array_to_sparse(a)
(::Type{SparseMatrixCSC})(a::AFArray{T}) where T = convert_array_to_sparse(a)

deepcopy_internal(a::AFArray{T,N}, d::IdDict) where {T,N} = haskey(d, a) ? d[a]::AFArray{T,N} : d[a] = copy(a)

import Base: size, eltype, ndims, abs, acos, acosh, asin, asinh, atan, atanh, cbrt, ceil, clamp, cos, cosh
import Base: count, div, exp, expm1, floor, hypot
import Base: identity, imag, isinf, isnan, iszero, join, log, log10, log1p, log2, maximum
import Base: minimum, mod, prod, randn, range, real, rem, replace, round, show, inv
import Base: sign, signbit, sin, sinh, sort, sortperm, sqrt, sum, tan, tanh, transpose, trunc, any, all
import Base: cat, hcat, vcat, max, min, sizeof, similar, length, sizeof
import Base: isfinite, ifelse, isempty
import LinearAlgebra: lu, rank, det, norm, diag, diagm, svd, cholesky, dot, qr
import Statistics: cov, std, var, mean, median
import SpecialFunctions: factorial, lgamma

similar(a::AFArray) = zeros(a)
similar(a::AFArray, ::Type{T}) where {T} = zeros(AFArray{T}, size(a))
sizeof(a::AFArray) = length(a) * sizeof(eltype(a))
eltype(a::AFArray{T,N}) where {T,N} = T
ndims(a::AFArray{T,N}) where {T,N} = N
size(a::AFVector) = (s = get_dims(a); (s[1],))
size(a::AFMatrix) = (s = get_dims(a); (s[1],s[2]))
size(a::AFVolume) = (s = get_dims(a); (s[1],s[2],s[3]))
size(a::AFTensor) = (s = get_dims(a); (s[1],s[2],s[3],s[4]))
size(a::AFArray, dim::Int) = get_dims(a)[dim]
any(a::AFArray, dims::Colon) = any_true_all(a)[1] == 1
all(a::AFArray, dims::Colon) = all_true_all(a)[1] == 1
any(f, a::AFArray) = any(f(a))
all(f, a::AFArray) = all(f(a))
maximum(a::AFArray{T}) where {T<:Real} = T(max_all(a)[1])
minimum(a::AFArray{T}) where {T<:Real} = T(min_all(a)[1])
mean(a::AFArray{T}) where {T<:Real} = T(mean_all(a)[1])
std(a::AFArray{T}) where {T<:Real} = T(sqrt(var_all(a, false)[1]))
var(a::AFArray{T}) where {T<:Real} = T(var_all(a, false)[1])
median(a::AFArray{T}) where {T<:Real} = T(median_all(a)[1])
prod(a::AFArray{T}) where {T<:Real} = T(product_all(a)[1])
sum(a::AFArray{UInt8,N}) where N = UInt32(sum_all(a)[1])
sum(a::AFArray{Bool,N}) where N = Int64(sum_all(a)[1])
sum(a::AFArray{T,N}) where {T<:Real,N} = T(sum_all(a)[1])
sum(a::AFArray{T,N}) where {T<:Complex,N} = (s = sum_all(a); T(s[1] + s[2]im))
real(a::AFArray{T}) where {T<:Real} = a
imag(a::AFArray{T}) where {T<:Real} = zeros(a)
length(a::AFArray) = prod(size(a))
inv(X::AFArray{T,N}) where {T,N} = inverse(X,AF_MAT_NONE)
range(::Type{AFArray{T}}, a::Integer, b::Integer) where T = range(1, [b], 0, T) + T(a)
function range(::Type{AFArray{T1}}, a::T2, b::T2, c::Integer) where {T1, T2}
    x = b .* ones(AFArray{T1}, c)
    x[1] = a
    cumsum(x)
end
isfinite(a::AFArray) = !isinf(a) & !isnan(a)
isempty(a::AFArray) = (length(a) == 0)

import Base: /, *, +, -, ^, ==, <, >, <=, >=, !, !=, &, |, <<, >>, xor

-(a::AFArray{T}) where T  = T(0) - a
!(a::AFArray) = not(a)

+(a::Number, b::AFArray)  = add(constant(a, size(b)), b, false)
-(a::Number, b::AFArray)  = sub(constant(a, size(b)), b, false)
*(a::Number, b::AFArray)  = mul(constant(a, size(b)), b, false)
/(a::Number, b::AFArray)  = div(constant(a, size(b)), b, false)
^(a::Number, b::AFArray)  = pow(constant(a, size(b)), b, false)
==(a::Number, b::AFArray) = eq(constant(a, size(b)),  b, false)
<(a::Number, b::AFArray)  = lt(constant(a, size(b)),  b, false)
>(a::Number, b::AFArray)  = gt(constant(a, size(b)),  b, false)
!=(a::Number, b::AFArray) = neq(constant(a, size(b)), b, false)
<=(a::Number, b::AFArray) = le(constant(a, size(b)),  b, false)
>=(a::Number, b::AFArray) = ge(constant(a, size(b)),  b, false)
(&)(a::Bool, b::AFArray)     = and(constant(a, size(b)),       b, false)
|(a::Bool, b::AFArray)     = or(constant(a, size(b)),        b, false)
(&)(a::Integer, b::AFArray)   = bitand(constant(a, size(b)),    b, false)
|(a::Integer, b::AFArray)   = bitor(constant(a, size(b)),     b, false)
<<(a::Integer, b::AFArray)  = bitshiftl(constant(a, size(b)), b, false)
>>(a::Integer, b::AFArray)  = bitshiftr(constant(a, size(b)), b, false)
xor(a::Integer, b::AFArray) = bitxor(constant(a, size(b)),    b, false)
max(a::Real, b::AFArray) = maxof(constant(a, size(b)),    b, false)
min(a::Real, b::AFArray) = minof(constant(a, size(b)),    b, false)

+(a::AFArray, b::Number)  = add(a, constant(b, size(a)), false)
-(a::AFArray, b::Number)  = sub(a, constant(b, size(a)), false)
*(a::AFArray, b::Number)  = mul(a, constant(b, size(a)), false)
/(a::AFArray, b::Number)  = div(a, constant(b, size(a)), false)
^(a::AFArray, b::Number)  = pow(a, constant(b, size(a)), false)
^(a::AFArray, b::Real)    = pow(a, constant(b, size(a)), false)
^(a::AFArray, b::Integer) = pow(a, constant(b, size(a)), false)
==(a::AFArray, b::Number) = eq(a,  constant(b, size(a)), false)
<(a::AFArray, b::Number)  = lt(a,  constant(b, size(a)), false)
>(a::AFArray, b::Number)  = gt(a,  constant(b, size(a)), false)
!=(a::AFArray, b::Number) = neq(a, constant(b, size(a)), false)
<=(a::AFArray, b::Number) = le(a,  constant(b, size(a)), false)
>=(a::AFArray, b::Number) = ge(a,  constant(b, size(a)), false)
(&)(a::AFArray, b::Bool)     = and(a,       constant(b, size(a)), false)
|(a::AFArray, b::Bool)     = or(a,        constant(b, size(a)), false)
(&)(a::AFArray, b::Integer)   = bitand(a,    constant(b, size(a)), false)
|(a::AFArray, b::Integer)   = bitor(a,     constant(b, size(a)), false)
<<(a::AFArray, b::Integer)  = bitshiftl(a, constant(b, size(a)), false)
>>(a::AFArray, b::Integer)  = bitshiftr(a, constant(b, size(a)), false)
xor(a::AFArray, b::Integer) = bitxor(a,    constant(b, size(a)), false)
max(a::AFArray, b::Real) = maxof(a,    constant(b, size(a)), false)
min(a::AFArray, b::Real) = minof(a,    constant(b, size(a)), false)

+(a::AFArray, b::AFArray) = add(a, b, bcast[])
-(a::AFArray, b::AFArray) = sub(a, b, bcast[])
*(a::AFArray, b::AFArray) = bcast[] ? mul(a, b, true) : A_mul_B(a,b)
/(a::AFArray, b::AFArray) = div(a, b, bcast[])
^(a::AFArray, b::AFArray) = pow(a, b, bcast[])
==(a::AFArray, b::AFArray) = bcast[] ? eq(a, b, true) : size(a) == size(b) && all(eq(a, b, false))
<(a::AFArray, b::AFArray) = lt(a, b, bcast[])
>(a::AFArray, b::AFArray) = gt(a, b, bcast[])
!=(a::AFArray, b::AFArray) = neq(a, b, bcast[])
<=(a::AFArray, b::AFArray) = le(a, b, bcast[])
>=(a::AFArray, b::AFArray) = ge(a, b, bcast[])
(&)(a::AFArray{Bool}, b::AFArray{Bool})  = and(a, b, bcast[])
|(a::AFArray{Bool}, b::AFArray{Bool})  = or(a, b, bcast[])
(&)(a::AFArray, b::AFArray)   = bitand(a, b, bcast[])
|(a::AFArray, b::AFArray)   = bitor(a, b, bcast[])
<<(a::AFArray, b::AFArray)  = bitshiftl(a, b, bcast[])
>>(a::AFArray, b::AFArray)  = bitshiftr(a, b, bcast[])
xor(a::AFArray, b::AFArray) = bitxor(a, b, bcast[])
max(a::AFArray, b::AFArray) = maxof(a, b, bcast[])
min(a::AFArray, b::AFArray) = minof(a, b, bcast[])

import Base: transpose, vec, reshape, adjoint
export A_mul_B, Ac_mul_B, At_mul_B, A_mul_Bc, Ac_mul_Bc, A_mul_Bt, At_mul_Bt, ctranspose

A_mul_B(a::AFArray,   b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_NONE)
Ac_mul_B(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_CTRANS, AF_MAT_NONE)
At_mul_B(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_TRANS,  AF_MAT_NONE)
A_mul_Bc(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_CTRANS)
Ac_mul_Bc(a::AFArray, b::AFArray) = matmul(a, b, AF_MAT_CTRANS, AF_MAT_CTRANS)
A_mul_Bt(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_TRANS)
At_mul_Bt(a::AFArray, b::AFArray) = matmul(a, b, AF_MAT_TRANS,  AF_MAT_TRANS)

sign(a::AFArray{T,N}) where {T,N} = (AFArray{T,N}(0<a) - AFArray{T,N}(a<0))

function At_mul_B(a::AFVector{T1},  b::AFArray{T2}) where {T1, T2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),
                 out,a.arr,b.arr,AF_MAT_TRANS, AF_MAT_NONE))
    AFArray{typed(T1,T2),1}(out[])
end

function A_mul_B(a::AFMatrix{T1}, b::AFVector{T2}) where {T1,T2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),
                 out,a.arr,b.arr,AF_MAT_NONE, AF_MAT_NONE))
    AFArray{typed(T1,T2),1}(out[])
end

function transpose(_in::AFArray{T,N},conjugate::Bool=false) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_transpose,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,conjugate))
    AFArray{T,2}(out[])
end
ctranspose(in::AFArray) = transpose(in, true)
adjoint(in::AFArray) = transpose(in, true)

function vec(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_flat,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,1}(out[])
end

function reshape(_in::AFArray{T},dims::NTuple{N,Int}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_moddims,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{dim_t}),out,_in.arr,UInt32(length(dims)),[dims...]))
    AFArray{T,N}(out[])
end
reshape(a::AFArray, t::Int...) = reshape(a, t)

using SpecialFunctions
import SpecialFunctions: erf, erfc

function Base.Broadcast.broadcasted(f, A::AFArray, Bs...)
    bcast[] =  true
    try
        return f(A, Bs...)
    finally
        bcast[] = false
    end
end

function Base.Broadcast.broadcasted(f, A::Number, B::AFArray, Cs...)
    bcast[] =  true
    try
        return f(A, B, Cs...)
    finally
        bcast[] = false
    end
end

# copy!(a::AFArray, b::AFArray) = (a.=b; b)

# function broadcast!(::typeof(identity), a::AFArray, b::AFArray)
#     write_array(a, get_device_ptr(b), UInt(sizeof(b)), afDevice)
#     unlock_device_ptr(b)
#     b
# end

# function broadcast!(::typeof(identity), a::Array, b::AFArray)
#     get_data_ptr(a, b)
#     b
# end

# function broadcast!(::typeof(identity), a::AFArray, b::Array)
#     write_array(a, b, UInt(sizeof(b)), afHost)
#     b
# end

# function broadcast!(f, C::AFArray, A::Array, Bs...)
#     bcast[] =  true
#     try
#         r = f(A, Bs...)
#         write_array(C, r, UInt(sizeof(r)), afHost)
#         return r
#     finally
#         bcast[] = false
#     end
# end

# function broadcast!(f, C::AFArray, A::AFArray, Bs...)
#     bcast[] =  true
#     try
#         swap!(C, f(A, Bs...))
#         return C
#     finally
#         bcast[] = false
#     end
# end

# function broadcast!(f, C::Array, A::AFArray, Bs...)
#     bcast[] =  true
#     try
#         r = f(A, Bs...)
#         C .= r
#         return r
#     finally
#         bcast[] = false
#     end
# end

import Base: fill, zeros, ones, zero, one

fill(::Type{AFArray}, a, dims::Int...) = constant(a, dims)
fill(::Type{AFArray{T}}, a, dims::Int...) where {T} = constant(T(a), dims)
fill(::Type{AFArray{T,N}}, a, dims::Int...) where {T,N} = constant(T(a), dims)
fill(::Type{AFArray}, a, dims::NTuple{N,Int}) where {N} = constant(a, dims)
fill(::Type{AFArray{T}}, a, dims::NTuple{N,Int}) where {T,N} = constant(T(a), dims)
fill(::Type{AFArray{T,N}}, a, dims::NTuple{N,Int}) where {T,N} = constant(T(a), dims)

zeros(::Type{AFArray}, dims::Int...) = constant(0., dims)
zeros(::Type{AFArray{T}}, dims::Int...) where {T} = constant(T(0), dims)
zeros(::Type{AFArray{T,N}}, dims::Int...) where {T,N} = constant(T(0), dims)
zeros(::Type{AFArray{T}}, dims::NTuple{N,Int}) where {T,N} = constant(T(0), dims)
zeros(::Type{AFArray{T,N}}, dims::NTuple{N,Int}) where {T,N} = constant(T(0), dims)
zeros(a::AFArray{T,N}) where {T,N} = constant(T(0), size(a))

ones(::Type{AFArray}, dims::Int...) = constant(1., dims)
ones(::Type{AFArray{T}}, dims::Int...) where {T} = constant(T(1), dims)
ones(::Type{AFArray{T,N}}, dims::Int...) where {T,N} = constant(T(1), dims)
ones(::Type{AFArray{T}}, dims::NTuple{N,Int}) where {T,N} = constant(T(1), dims)
ones(::Type{AFArray{T,N}}, dims::NTuple{N,Int}) where {T,N} = constant(T(1), dims)
ones(a::AFArray{T,N}) where {T,N} = constant(T(1), size(a))

zero(a::AFArray{T,N}) where {T,N} = constant(T(0), size(a))
function one(a::AFArray{T,N}) where {T,N}
    if bcast[]
        ones(a)
    else
        out = RefValue{af_array}(0)
        _error(ccall((:af_identity,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),
                     out,UInt32(N),[size(a)...],af_type(T)))
        AFArray{T,N}(out[])
    end
end

import Base.\

\(a::AFArray, b::AFArray) = solve(a, b, AF_MAT_NONE)

export swap!
function swap!(a::AFArray{T,N}, b::AFArray{T,N}) where {T,N}
    a.arr, b.arr = b.arr, a.arr
    nothing
end

function abs(_in::AFArray{Complex{T},N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_abs,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

eye(a::AFArray{T,2}) where T = identity(2, [size(a)...], T)
eye(::Type{AFArray{T}}, n::Int) where T = identity(2, [n, n], T)

import Base: cumsum, cumprod

for op in (:minimum, :maximum, :sum, :any, :all)
    @eval $op(a::AFArray{T,N}; dims=:) where {T,N} = $op(a, dims)
end

for op in (:sort, :cumsum, :cumprod, :cummin, :cummax)
    @eval $op(a::AFArray{T,N}; dims=1) where {T,N} = $op(a, dims)
end

function clamp(_in::AFArray{T1,N1},lo::N2,hi::N2) where {T1,N1,N2}
    return clamp(_in, constant(lo, size(_in)), constant(hi, size(_in)))
end
