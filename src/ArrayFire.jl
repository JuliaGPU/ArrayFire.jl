module ArrayFire

using Cxx
using Base.Meta

import Base: rand, show, randn, ones, diag, eltype, size, elsize,
    sizeof, length, showarray, convert, ndims, lu, qr, svd
import Cxx: CppEnum
export AFArray, chol!, constant

# If you have a crash, enable this
const AF_DEBUG = true

function init_library()
	Libdl.dlopen(Pkg.dir("ArrayFire","deps","build","arrayfire","src","backend","opencl","libafopencl"),Libdl.RTLD_GLOBAL)
	addHeaderDir(Pkg.dir("ArrayFire","deps","src","arrayfire","include");kind = C_System)
end
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
    function AFArray(array::vcpp"af::array")
        ret = new(array)
        if AF_DEBUG && backend_eltype(ret) != T
            error("Tried to create AFArray{$T} with array of eltype $(backend_eltype(ret))")
        end
        ret
    end
end
immutable AFSubArray{T} <: AFAbstractArray{T}
    array::vcpp"af::array::array_proxy"
    function AFSubArray(array::vcpp"af::array::array_proxy")
        ret = new(array)
        if AF_DEBUG && backend_eltype(ret) != T
            error("Tried to create AFArray{$T} with array of eltype $(backend_eltype(ret))")
        end
        ret
    end
end

AFArray() = icxx"af::array();"

call{T}(::Type{AFAbstractArray{T}}, proxy::vcpp"af::array::array_proxy") =
    AFSubArray{T}(proxy)
call{T}(::Type{AFAbstractArray{T}}, arr::vcpp"af::array") =
    AFArray{T}(arr)

Cxx.cppconvert{T}(x::AFAbstractArray{T}) = x.array

convert{T}(::Type{AFArray{T}}, arr::AFSubArray{T}) =
    AFArray{T}(icxx"(af::array)$arr")

eltype{T}(x::AFAbstractArray{T}) = T
backend_eltype(x) = jltype(icxx"$x.type();")
sizeof{T}(a::AFAbstractArray{T}) = elsize(a) * length(a)

# GPU to Host
function convert{T}(::Type{Array{T}},x::AFAbstractArray{T})
    ret = Array(UInt8,sizeof(x))
    icxx"$x.host($(pointer(ret)));"
    ret = reinterpret(T, ret)
    ret = reshape(ret, size(x)...)
    ret
end
convert{T}(::Type{Array},x::AFAbstractArray{T}) = convert(Array{T},x)

# Host to GPU
function convert{T}(::Type{AFArray{T}}, x::Array{T})
    AFArray{T}(icxx"return af::array{$(dims_to_dim4(size(x))),$(pointer(x)),afHost};")
end
convert{T}(::Type{AFArray}, x::Array{T}) = convert(AFArray{T}, x)



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
        icxx"af::dim4{$(dims[1])};"
    elseif length(dims) == 2
        icxx"af::dim4{$(dims[1]),$(dims[2])};"
    elseif length(dims) == 3
        icxx"af::dim4{$(dims[1]),$(dims[2]),$(dims[3])};"
    elseif length(dims) == 4
        icxx"af::dim4{$(dims[1]),$(dims[2]),$(dims[3]),$(dims[4])};"
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

#Functions to generate arrays
constant(val::Real, dims::Integer...) = AFArray{typeof(val)}(icxx"af::constant($val, $(dims_to_dim4(dims)), $(aftype(typeof(val))));")
rand{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::randu($(dims_to_dim4(dims)),$(aftype(T)));")
randn{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::randn($(dims_to_dim4(dims)),$(aftype(T)));")
eye{T}(::Type{AFArray{T}}, dims...) = AFArray{T}(icxx"af::identity($(dims_to_dim4(dims)),$(aftype(T)));")
diag{T}(x::AFArray{T}, k = 0) = AFArray{T}(icxx"af::diag($(dims_to_dim4(dims)),$k);");

import Base: getindex

# Note that we need to translate between 0 and 1 based indexing
to_af_idx(x::Real) = convert(Cint, x-1)
to_af_idx(x::Range) = icxx"af::seq($(first(x))-1,$(last(x))-1,$(step(x)));"
to_af_idx(x::Colon) = icxx"af::span;"
const tis = to_af_idx

IS = Union{Real,Range,Colon}

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
    num_elems = icxx"$proxy.elements();"
    if num_elems == 1
        return Array(AFSubArray{T}(proxy))[1]
    else
        return AFSubArray{T}(proxy)
    end
end

#TODO : return 0-element array when b is all falses
function getindex(x::AFAbstractArray, b::AFAbstractArray{Bool})
    out = AFArray()
    if any(b)
        icxx"$out = $x($b);"
        AFArray{backend_eltype(out)}(out)
    else
        return 0
    end
end

function setindex!{T}(x::AFAbstractArray{T}, val, idxs...)
    proxy = _getindex(x,idxs...)
    icxx"$proxy = $val;"
end

# Avoid crash when referencing with boolean arrays with all falses 
function setindex!(x::AFAbstractArray, val, b::AFAbstractArray{Bool})
    if any(b) 
        icxx"$x($b) = $val;"
    else
        return val
    end
end

#Helper functions
import Base: isnan

function isnan(a::AFAbstractArray)
    AFArray{Bool}(icxx"af::isNaN($a);")
end

# Fourier Transforms
# TODO: Multidimensional
import Base: fft
function fft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{Complex{T}}(icxx"fft2($a);")
	end
end
		
function fft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{T}(icxx"fft2($a);")
	end
end

# Exception handling
import Base: showerror
@exception function showerror(io::IO, e::rcpp"af::exception")
    println(io, bytestring(icxx"$e.what();"))
end

#import other files
include("math.jl")
include("image.jl")
include("stats.jl")
include("vector.jl")
include("linalg.jl")

end # module
