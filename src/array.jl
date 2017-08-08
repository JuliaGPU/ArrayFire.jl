if VERSION < v"0.6-"
    @eval begin
        type AFArray{T,N} <: AbstractArray{T,N}
            arr::af_array
            function AFArray(arr::af_array)
                a = new(arr)
                finalizer(a, release_array)
                if !isempty(scopes)
                    push!(scopes[end], a)
                end
                a
            end
        end
    end
else
    include_string("
        mutable struct AFArray{T,N} <: AbstractArray{T,N}
            arr::af_array
            function AFArray{T,N}(arr::af_array) where {T,N}
                # @assert get_type(arr) == T
                a = new{T,N}(arr)
                finalizer(a, release_array)
                if !isempty(scopes)
                    push!(scopes[end], a)
                end
                a
            end
        end")
end

@compat AFVector{T} = AFArray{T,1}
@compat AFMatrix{T} = AFArray{T,2}
@compat AFVolume{T} = AFArray{T,3}
@compat AFTensor{T} = AFArray{T,4}

export AFArray, AFVector, AFMatrix, AFVolume, AFTensor

import Base: Array, SparseMatrixCSC, copy, deepcopy_internal, issparse, sparse, full, complex, conj

sparse{T,N}(a::AFArray{T,N}) = create_sparse_array_from_dense(a, AF_STORAGE_CSR)

(::Type{AFArray{T1}}){T1,T2,N}(a::AFArray{T2,N}) = recast_array(AFArray{T1}, a)
(::Type{AFArray{T1,N}}){T1,T2,N}(a::AFArray{T2,N}) = recast_array(AFArray{T1}, a)

(::Type{Array{T,N}}){T,N}(a::AFArray{T,N}) = convert_array(a)
(::Type{AFArray{T,N}}){T,N}(a::Array{T,N}) = convert_array(a)

(::Type{Array{T}}){T,N}(a::AFArray{T,N}) = convert_array(a)
(::Type{AFArray{T}}){T,N}(a::Array{T,N}) = convert_array(a)

(::Type{Array}){T,N}(a::AFArray{T,N}) = convert_array(a)
(::Type{AFArray}){T,N}(a::Array{T,N}) = convert_array(a)

(::Type{AFArray})(a::SparseMatrixCSC) = convert_array_to_sparse(a)

(::Type{SparseMatrixCSC{T}}){T}(a::AFArray{T}) = convert_array_to_sparse(a)
(::Type{SparseMatrixCSC}){T}(a::AFArray{T}) = convert_array_to_sparse(a)

deepcopy_internal{T,N}(a::AFArray{T,N}, d::ObjectIdDict) = haskey(d, a) ? d[a]::AFArray{T,N} : copy(a)

import Base: size, eltype, ndims, abs, acos, acosh, asin, asinh, atan, atan2, atanh, cbrt, ceil, clamp, cos, cosh
import Base: count, cov, det, div, dot, exp, expm1, factorial, fft, floor, gradient, hypot
import Base: identity, ifft, imag, isinf, isnan, iszero, join, lgamma, log, log10, log1p, log2, lu, maximum, mean, median
import Base: minimum, mod, norm, prod, qr, randn, range, rank, real, rem, replace, round, select, show
import Base: sign, signbit, sin, sinh, sort, sortperm, std, sqrt, sum, svd, tan, tanh, transpose, trunc, var, any, all
import Base: cat, hcat, vcat, conv, max, min, sizeof

sizeof(a::AFArray) = length(a) * sizeof(eltype(a))
eltype{T,N}(a::AFArray{T,N}) = T
ndims{T,N}(a::AFArray{T,N}) = N
size(a::AFVector) = (s = get_dims(a); (s[1],))
size(a::AFMatrix) = (s = get_dims(a); (s[1],s[2]))
size(a::AFVolume) = (s = get_dims(a); (s[1],s[2],s[3]))
size(a::AFTensor) = (s = get_dims(a); (s[1],s[2],s[3],s[4]))
size(a::AFArray, dim::Int) = get_dims(a)[dim]
any(a::AFArray) = any_true_all(a)[1] == 1
all(a::AFArray) = all_true_all(a)[1] == 1
maximum{T<:Real}(a::AFArray{T})::T = max_all(a)[1]
minimum{T<:Real}(a::AFArray{T})::T = min_all(a)[1]
mean{T<:Real}(a::AFArray{T})::T = mean_all(a)[1]
std{T<:Real}(a::AFArray{T})::T = sqrt(var_all(a, false)[1])
var{T<:Real}(a::AFArray{T})::T = var_all(a, false)[1]
median{T<:Real}(a::AFArray{T})::T = median_all(a)[1]
prod{T<:Real}(a::AFArray{T})::T = product_all(a)[1]
sum{N}(a::AFArray{UInt8,N}) = UInt32(sum_all(a)[1])
sum{N}(a::AFArray{Bool,N}) = Int64(sum_all(a)[1])
sum{T<:Real,N}(a::AFArray{T,N})::T = sum_all(a)[1]
sum{T<:Complex,N}(a::AFArray{T,N})::T = (s = sum_all(a); s[1] + s[2]im)

import Base: /, *, +, -, ^, ==, <, >, <=, >=, !, !=, &, |, <<, >>, xor

-{T}(a::AFArray{T})       = T(0) - a
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

import Base: Ac_mul_B, At_mul_B, A_mul_Bc, Ac_mul_Bc, A_mul_Bt, At_mul_Bt, transpose, ctranspose, vec, reshape
export A_mul_B

A_mul_B(a::AFArray,   b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_NONE)
Ac_mul_B(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_CTRANS, AF_MAT_NONE)
At_mul_B(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_TRANS,  AF_MAT_NONE)
A_mul_Bc(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_CTRANS)
Ac_mul_Bc(a::AFArray, b::AFArray) = matmul(a, b, AF_MAT_CTRANS, AF_MAT_CTRANS)
A_mul_Bt(a::AFArray,  b::AFArray) = matmul(a, b, AF_MAT_NONE,   AF_MAT_TRANS)
At_mul_Bt(a::AFArray, b::AFArray) = matmul(a, b, AF_MAT_TRANS,  AF_MAT_TRANS)

sign{T,N}(a::AFArray{T,N}) = (AFArray{T,N}(0<a) - AFArray{T,N}(a<0))

function At_mul_B{T1, T2}(a::AFVector{T1},  b::AFArray{T2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),
                 out,a.arr,b.arr,AF_MAT_TRANS, AF_MAT_NONE))
    AFArray{typed(T1,T2),1}(out[])
end

function A_mul_B{T1,T2}(a::AFMatrix{T1}, b::AFVector{T2})
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),
                 out,a.arr,b.arr,AF_MAT_NONE, AF_MAT_NONE))
    AFArray{typed(T1,T2),1}(out[])
end

function transpose{T,N}(_in::AFArray{T,N},conjugate::Bool=false)
    out = RefValue{af_array}(0)
    _error(ccall((:af_transpose,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,conjugate))
    AFArray{T,2}(out[])
end
ctranspose(in::AFArray) = transpose(in, true)

function vec{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_flat,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,1}(out[])
end

function reshape{T,N}(_in::AFArray{T},dims::NTuple{N,Int})
    out = RefValue{af_array}(0)
    _error(ccall((:af_moddims,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{dim_t}),out,_in.arr,UInt32(length(dims)),[dims...]))
    AFArray{T,N}(out[])
end
reshape(a::AFArray, t::Int...) = reshape(a, t)

if VERSION < v"0.6-"
    @eval begin
        import Base: ./, .*, .+, .-, .^, .==, .<, .>, .<=, .>=, .!=, .<<, .>>
            export xor

        .-{T}(a::AFArray{T})       = T(0) - a

        .+(a::Number, b::AFArray)  = add(constant(a, size(b)), b, false)
        .-(a::Number, b::AFArray)  = sub(constant(a, size(b)), b, false)
        .*(a::Number, b::AFArray)  = mul(constant(a, size(b)), b, false)
        ./(a::Number, b::AFArray)  = div(constant(a, size(b)), b, false)
        .^(a::Number, b::AFArray)  = pow(constant(a, size(b)), b, false)
        .==(a::Number, b::AFArray) = eq(constant(a, size(b)),  b, false)
        .<(a::Number, b::AFArray)  = lt(constant(a, size(b)),  b, false)
        .>(a::Number, b::AFArray)  = gt(constant(a, size(b)),  b, false)
        .!=(a::Number, b::AFArray) = neq(constant(a, size(b)), b, false)
        .<=(a::Number, b::AFArray) = le(constant(a, size(b)),  b, false)
        .>=(a::Number, b::AFArray) = ge(constant(a, size(b)),  b, false)
        .<<(a::Integer, b::AFArray)  = bitshiftl(constant(a, size(b)), b, false)
        .>>(a::Integer, b::AFArray)  = bitshiftr(constant(a, size(b)), b, false)

        .+(a::AFArray, b::Number)  = add(a, constant(b, size(a)), false)
        .-(a::AFArray, b::Number)  = sub(a, constant(b, size(a)), false)
        .*(a::AFArray, b::Number)  = mul(a, constant(b, size(a)), false)
        ./(a::AFArray, b::Number)  = div(a, constant(b, size(a)), false)
        .^(a::AFArray, b::Number)  = pow(a, constant(b, size(a)), false)
        .==(a::AFArray, b::Number) = eq(a,  constant(b, size(a)), false)
        .<(a::AFArray, b::Number)  = lt(a,  constant(b, size(a)), false)
        .>(a::AFArray, b::Number)  = gt(a,  constant(b, size(a)), false)
        .!=(a::AFArray, b::Number) = neq(a, constant(b, size(a)), false)
        .<=(a::AFArray, b::Number) = le(a,  constant(b, size(a)), false)
        .>=(a::AFArray, b::Number) = ge(a,  constant(b, size(a)), false)
        .<<(a::AFArray, b::Integer)  = bitshiftl(a, constant(b, size(a)), false)
        .>>(a::AFArray, b::Integer)  = bitshiftr(a, constant(b, size(a)), false)

        .+(a::AFArray, b::AFArray) = add(a, b, true)
        .-(a::AFArray, b::AFArray) = sub(a, b, true)
        .*(a::AFArray, b::AFArray) = mul(a, b, true)
        ./(a::AFArray, b::AFArray) = div(a, b, true)
        .^(a::AFArray, b::AFArray) = pow(a, b, true)
        .==(a::AFArray, b::AFArray) = eq(a, b, true)
        .<(a::AFArray, b::AFArray) = lt(a, b, true)
        .>(a::AFArray, b::AFArray) = gt(a, b, true)
        .!=(a::AFArray, b::AFArray) = neq(a, b, true)
        .<=(a::AFArray, b::AFArray) = le(a, b, true)
        .>=(a::AFArray, b::AFArray) = ge(a, b, true)
        .<<(a::AFArray, b::AFArray)  = bitshiftl(a, b, true)
        .>>(a::AFArray, b::AFArray)  = bitshiftr(a, b, true)
        xor(a, b) = a $ b
    end
    import Base: erf, erfc, broadcast

    function broadcast(f, A::AFArray, Bs...)
        bcast[] =  true
        try
            return f(A, Bs...)
        finally
            bcast[] = false
        end
    end

    function broadcast(f, A0::Number, A::AFArray, Bs...)
        bcast[] =  true
        try
            return f(A, Bs...)
        finally
            bcast[] = false
        end
    end

else
    using SpecialFunctions
    import SpecialFunctions: erf, erfc

    import Base.Broadcast: promote_containertype, broadcast_c, _containertype, broadcast_c!

    _containertype(::Type{<:AFArray}) = AFArray

    promote_containertype(::Type{AFArray}, ::Type{AFArray}) = AFArray
    promote_containertype(::Type{AFArray}, ct) = AFArray
    promote_containertype(ct, ::Type{AFArray}) = AFArray

    function broadcast_c(f, ::Type{AFArray}, A, Bs...)
        bcast[] =  true
        try
            return f(A, Bs...)
        finally
            bcast[] = false
        end
    end

    function broadcast_c!(f, ::Type{Array}, ::Type{AFArray}, C, A, Bs...)
        bcast[] =  true
        try
            r = f(A, Bs...)
            get_data_ptr(C, r)
            return r
        finally
            bcast[] = false
        end
    end

    function broadcast_c!(f, ::Type{AFArray}, ::Type{Array}, C, A, Bs...)
        bcast[] =  true
        try
            r = f(A, Bs...)
            write_array(C, r, UInt(sizeof(r)), afHost)
            return r
        finally
            bcast[] = false
        end
    end

    function broadcast_c!(f, ::Type{AFArray}, ::Type{AFArray}, C, A, Bs...)
        bcast[] =  true
        try
            swap!(C, f(A, Bs...))
            return C
        finally
            bcast[] = false
        end
    end
end

import Base: fill, zeros, ones

fill(::Type{AFArray}, a, dims::Int...) = constant(a, dims)
fill{T}(::Type{AFArray{T}}, a, dims::Int...) = constant(T(a), dims)
fill{T,N}(::Type{AFArray{T,N}}, a, dims::Int...) = constant(T(a), dims)
fill{N}(::Type{AFArray}, a, dims::NTuple{N,Int}) = constant(a, dims)
fill{T,N}(::Type{AFArray{T}}, a, dims::NTuple{N,Int}) = constant(T(a), dims)
fill{T,N}(::Type{AFArray{T,N}}, a, dims::NTuple{N,Int}) = constant(T(a), dims)

zeros(::Type{AFArray}, dims::Int...) = constant(0., dims)
zeros{T}(::Type{AFArray{T}}, dims::Int...) = constant(T(0), dims)
zeros{T,N}(::Type{AFArray{T,N}}, dims::Int...) = constant(T(0), dims)
zeros{T,N}(::Type{AFArray{T}}, dims::NTuple{N,Int}) = constant(T(0), dims)
zeros{T,N}(::Type{AFArray{T,N}}, dims::NTuple{N,Int}) = constant(T(0), dims)
zeros{T,N}(a::AFArray{T,N}) = constant(T(0), size(a))

ones(::Type{AFArray}, dims::Int...) = constant(1., dims)
ones{T}(::Type{AFArray{T}}, dims::Int...) = constant(T(1), dims)
ones{T,N}(::Type{AFArray{T,N}}, dims::Int...) = constant(T(1), dims)
ones{T,N}(::Type{AFArray{T}}, dims::NTuple{N,Int}) = constant(T(1), dims)
ones{T,N}(::Type{AFArray{T,N}}, dims::NTuple{N,Int}) = constant(T(1), dims)
ones{T,N}(a::AFArray{T,N}) = constant(T(1), size(a))

export swap!
function swap!{T,N}(a::AFArray{T,N}, b::AFArray{T,N})
    a.arr, b.arr = b.arr, a.arr
    nothing
end

function abs{T,N}(_in::AFArray{Complex{T},N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_abs,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end
