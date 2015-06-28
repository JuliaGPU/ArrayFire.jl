module ArrayFire

using Cxx
using Base.Meta

import Base: rand, show, randn, ones, diag, eltype, size, elsize,
    sizeof, length, showarray, convert, ndims
import Cxx: CppEnum
export AFArray

init_library() = Libdl.dlopen(Pkg.dir("ArrayFire","deps","src","arrayfire","src","backend","opencl","libafopencl.dylib"),Libdl.RTLD_GLOBAL)
init_library()

__init__() = init_library()

cxx"""
#include <arrayfire.h>
"""

setDevice(i) = icxx"af::setDevice($i);"

const f32 = CppEnum{:af_dtype}(0) # 32-bit floating point values
const c32 = CppEnum{:af_dtype}(1) # 32-bit complex floating point values
const f64 = CppEnum{:af_dtype}(2) # 64-bit complex floating point values
const c64 = CppEnum{:af_dtype}(3) # 64-bit complex floating point values
const b8  = CppEnum{:af_dtype}(4) #  8-bit boolean values
const s32 = CppEnum{:af_dtype}(5) # 32-bit signed integral values
const u32 = CppEnum{:af_dtype}(6) # 32-bit unsigned integral values
const u8  = CppEnum{:af_dtype}(7) #  8-bit unsigned integral values
const s64 = CppEnum{:af_dtype}(8) # 64-bit signed integral values
const u64 = CppEnum{:af_dtype}(9) # 64-bit unsigned integral values

aftype(::Type{Float32})          = f32
aftype(::Type{Complex{Float32}}) = c32
aftype(::Type{Float64})          = f64
aftype(::Type{Complex{Float64}}) = c64
aftype(::Type{Bool})             = b8
aftype(::Type{Int32})            = s32
aftype(::Type{UInt32})           = u32
aftype(::Type{UInt8})            = u8
aftype(::Type{Int64})            = s64
aftype(::Type{UInt64})           = u64

function jltype(af_dtype)
    if af_dtype == f32
        return Float32
    elseif af_dtype == c32
        return Complex{Float32}
    elseif af_dtype == f64
        return Float64
    elseif af_dtype == c64
        return Complex{Float64}
    elseif af_dtype == b8
        return Bool
    elseif af_dtype == s32
        return Int32
    elseif af_dtype == u32
        return UInt32
    elseif af_dtype == u8
        return UInt8
    elseif af_dtype == s64
        return Int64
    elseif af_dtype == u64
        return UInt64
    end
end

function af_promote{T,S}(::Type{T},::Type{S})
    if T == S
        return T
    elseif T == Complex{Float64} || S == Complex{Float64}
        return Complex{Float64}
    elseif T == Complex{Float32} || S == Complex{Float32}
        (T == Float64 || S == Float64) && return Complex{Float64}
        return Complex{Float32}
    elseif T == Float64 || S == Float64
        return Float64
    elseif T == Float32 || S == Float32
        return Float32
    elseif T == UInt64 || S == UInt64
        return UInt64
    elseif T == Int64 || S == Int64
        return Int64
    elseif T == UInt32 || S == UInt32
        return UInt32
    elseif T == Int32 || S == Int32
        return Int32
    elseif T == UInt8 || S == UInt8
        return UInt8
    elseif T == Bool || S == Bool
        return Bool
    else
        return Float32
    end
end

abstract AFAbstractArray{T} <: AbstractArray{T,4}

immutable AFArray{T} <: AFAbstractArray{T}
    array::vcpp"af::array"
end
immutable AFSubArray{T} <: AFAbstractArray{T}
    array::vcpp"af::array::array_proxy"
end

call{T}(::Type{AFAbstractArray{T}}, proxy::vcpp"af::array::array_proxy") =
    AFSubArray{T}(proxy)
call{T}(::Type{AFAbstractArray{T}}, arr::vcpp"af::array") =
    AFArray{T}(arr)

Cxx.cppconvert{T}(x::AFAbstractArray{T}) = x.array

convert{T}(::Type{AFArray{T}}, arr::AFSubArray{T}) =
    AFArray{T}(icxx"(af::array)$arr")

eltype{T}(x::AFAbstractArray{T}) = T
backend_eltype(x::AFAbstractArray) = jltype(icxx"$x.type();")
sizeof{T}(a::AFAbstractArray{T}) = elsize(a) * length(a)

function convert{T}(::Type{Array{T}},x::AFAbstractArray{T})
    ret = Array(Uint8,sizeof(x))
    icxx"$x.host($(pointer(ret)));"
    ret = reinterpret(T, ret)
    ret = reshape(ret, size(x)...)
    ret
end

# Show this array, by getting a data pointer to it and using julia's printing
# mechanism. ArrayFire copies internally anyway, so there is nothing to be gained
# by using it's printing
function showarray{T}(io::IO, X::AFAbstractArray{T};
                   header::Bool=true, kwargs...)
    header && print(io, summary(X))
    !isempty(X) && println(io,":")
    showarray(io, convert(Array{T},X); header = false, kwargs...)
end

af_print(X::AFAbstractArray) = icxx"""af::print("",$X);"""

# Translate dimensions. Note that we need to translate between 0 and 1 based
# indexing
ndims(dim4) = icxx"$dim4.ndims();"
ndims(A::AFAbstractArray) = ndims(icxx"$A.dims();")
function dims_to_dim4(dims)
    if length(dims) == 1
        icxx"af::dim4($(dims[1]));"
    elseif length(dims) == 2
        icxx"af::dim4($(dims[1]),$(dims[2]));"
    elseif length(dims) == 3
        icxx"af::dim4($(dims[1]),$(dims[2]),$(dims[3]));"
    elseif length(dims) == 4
        icxx"af::dim4($(dims[1]),$(dims[2]),$(dims[3]),$(dims[4]));"
    else
        error("Too many dimensions")
    end
end
function dim4_to_dims(dim4)
    d = ndims(dim4)
    t = (Int(icxx"(int)$dim4[0];"),)
    d == 1 && return t
    t = tuple(t...,Int(icxx"(int)$dim4[1];"))
    d == 2 && return t
    t = tuple(t...,Int(icxx"(int)$dim4[2];"))
    d == 3 && return t
    t = tuple(t...,Int(icxx"(int)$dim4[3];"))
    d == 4 && return t
    error("Bad dim4 object")
end

size(x::AFAbstractArray) = dim4_to_dims(icxx"($x).dims();")

rand{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::randu($(dims_to_dim4(dims)),$(aftype(T)));")
randn{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::randn($(dims_to_dim4(dims)),$(aftype(T)));")
eye{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::identity($(dims_to_dim4(dims)),$(aftype(T)));")
diag{T}(x::AFArray{T}, k = 0) = AFArray{T}(icxx"af::diag($(dims_to_dim4(dims)),$k);");

for op in (:sin, :cos, :tan, :asin, :acos, :log, :log1p, :log10, :sqrt, :transpose,
    :exp, :expm1, :erf, :erfc, :cbrt, :lgamma, :transpose)
    @eval Base.($(quot(op))){T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::($op)(x.array))
end

Base.gamma{T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::tgamma(x))

#
import Base: +, -

# Resolve conflicts
+(x::AFAbstractArray{Bool},y::Bool) = AFArray{Bool}(@cxx +(x.array,y))
+(y::Bool,x::AFAbstractArray{Bool}) = AFArray{Bool}(@cxx +(y,x.array))
-(x::AFAbstractArray{Bool},y::Bool) = AFArray{Bool}(@cxx -(x.array,y))
-(y::Bool,x::AFAbstractArray{Bool}) = AFArray{Bool}(@cxx -(y,x.array))


for (op,cppop) in ((:+,:+),(:(.+),:+),(:-,:-),(:(.-),:-),(:.*,:*),(:./,:/),(:.>>,:>>),(:.<<,:<<))
    @eval function Base.($(quot(op))){T,S}(x::AFAbstractArray{T}, y::AFAbstractArray{S})
        a1 = x.array
        a2 = y.array
        AFArray{af_promote(T,S)}(@cxx ($(cppop))(a1,a2))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(x::AFAbstractArray{T}, y::S)
        a = x.array
        AFArray{af_promote(T,S)}(@cxx ($(cppop))(a, y))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(y::S, x::AFAbstractArray{T})
        a = x.array
        AFArray{af_promote(T,S)}(@cxx ($(cppop))(y, a))
    end
end
# TODO: add! using +=, etc.

import Base: getindex

# Note that we need to translate between 0 and 1 based indexing
to_af_idx(x::Real) = convert(Cint, x-1)
to_af_idx(x::Range) = icxx"af::seq($(first(x))-1,$(last(x))-1,$(step(x)));"
to_af_idx(x::Colon) = icxx"af::span;"
const tis = to_af_idx

IS = Union(Real,Range,Colon)

_getindex{T}(x::AFAbstractArray{T},y::AFAbstractArray)   = icxx"$x($y-1);"
# Avoid an ambiguity
_getindex{T}(x::AFAbstractArray{T},s0::Real)             = icxx"$x($(tis(s0)));"

_getindex{T}(x::AFAbstractArray{T},s0::IS)               = icxx"$x($(tis(s0)));"

_getindex{T}(x::AFAbstractArray{T},s0::IS,s1::IS)        = icxx"$x($(tis(s0)),
                                                                   $(tis(s1)));"

_getindex{T}(x::AFAbstractArray{T},s0::IS,s1::IS,s2::IS) = icxx"$x($(tis(s0)),
                                                                   $(tis(s1)),
                                                                   $(tis(s2)));"
_getindex{T}(x::AFAbstractArray{T},
                    s0::IS,s1::IS,s2::IS, s3::IS)        = icxx"$x($(tis(s0)),
                                                                   $(tis(s1)),
                                                                   $(tis(s2)),
                                                                   $(tis(s3)));"

function getindex{T}(x::AFAbstractArray{T}, idxs...)
    proxy = _getindex(x,idxs...)
    AFSubArray{T}(proxy)
end

function setindex!{T}(x::AFAbstractArray{T}, val, idxs...)
    proxy = _getindex(x,idxs...)
    icxx"$proxy = $val;"
end

# BLAS operations

import Base: dot, A_mul_Bt, At_mul_B, At_mul_Bt, A_mul_Bc,
    Ac_mul_B, Ac_mul_Bc, transpose, ctranspose, transpose!, ctranspose!

dot{T,S}(lhs::AFAbstractArray{T}, rhs::AFAbstractArray{S}) =
    icxx"dot($lhs, $rhs);"

# Matmul
*(a::AFAbstractArray, b::AFAbstractArray) =
    icxx"matmult($a, $b);"
*(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray) =
    icxx"matmult($a, $b, $c);"
*(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::AFAbstractArray) =
    icxx"matmult($a, $b, $c, $d);"

function _matmult(a::AFAbstractArray, b::AFAbstractArray;
    lhsProp = AF_MAT_NONE, rhsProp = AF_MAT_NONE)
    icxx"matmul($a,$b,$lhsProp,$rhsProp);"
end

# with transpose
A_mul_Bt(a::AFAbstractArray, b::AFAbstractArray) =
    icxx"matmulNT($a,$b);"

At_mul_B(a::AFAbstractArray, b::AFAbstractArray) =
    icxx"matmulTN($a,$b);"

At_mul_Bt(a::AFAbstractArray, b::AFAbstractArray) =
    icxx"matmulTT($a,$b);"

# with complex conjugate
A_mul_Bc(a::AFAbstractArray, b::AFAbstractArray) =
    _matmult(a,b,rhsProp=AF_MAT_CTRANS)
Ac_mul_B(a::AFAbstractArray, b::AFAbstractArray) =
    _matmult(a,b,lhsProp=AF_MAT_CTRANS)
Ac_mul_Bc(a::AFAbstractArray, b::AFAbstractArray) =
    _matmult(a,b,lhsProp=AF_MAT_CTRANS,rhsProp=AF_MAT_CTRANS)

# transpose

transpose(x::AFAbstractArray) = icxx"af::transpose($x);"
ctranspose(x::AFAbstractArray) = icxx"af::transpose($x,true);"

transpose!(x::AFAbstractArray) = icxx"af::transposeInPlace($x);"
ctranspose!(x::AFAbstractArray) = icxx"af::transposeInPlace($x,true);"

# solve

# TODO: The documentation says only AF_MAT_LOWER/AF_MAT_UPPER are supported
# once AF_MAT_(C)TRANS is supported this could be useful for A_rdiv, etc
\(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::solve($a,$b);"

# Factorizations
import Base.LinAlg: chol, chol!, PosDefException, UpperTriangular,
    LowerTriangular, lufact!, lufact

function _chol{T}(a::AFAbstractArray{T}, is_upper::Bool)
    out = similar(a)
    info = icxx"af::cholesky($out,$a,$is_upper);"
    info > 0 && throw(PosDefException(info))
    AFArray{T}(out)
end

function _chol!{T}(a::AFAbstractArray{T}, is_upper::Bool)
    info = icxx"af::choleskyInPlace($a,$is_upper);"
    info > 0 && throw(PosDefException(info))
    a
end

chol(a::AFAbstractArray, ::Type{Val{:U}}) = UpperTriangular(_chol(a), true)
chol(a::AFAbstractArray, ::Type{Val{:L}}) = LowerTriangular(_chol(a), false)
chol(a::AFAbstractArray) = chol(A,Val{:U})

chol!(a::AFAbstractArray, ::Type{Val{:U}}) = UpperTriangular(_chol!(a), true)
chol!(a::AFAbstractArray, ::Type{Val{:L}}) = LowerTriangular(_chol!(a), false)
chol!(a::AFAbstractArray) = chol!(A,Val{:U})

# Fourier Transforms
# TODO: Multidimensional
import Base: fft
fft{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"fft($a);")

end # module
