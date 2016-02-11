# MATH OPERATIONS

for op in (:sin, :cos, :tan, :asin, :acos, :log, :log1p, :log10, :sqrt, :transpose,
    :exp, :expm1, :erf, :erfc, :cbrt, :lgamma, :transpose)
    @eval Base.($(quot(op))){T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::($op)(x.array))
end

Base.gamma{T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::tgamma(x))

#
import Base: +, -, abs

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
        # This is special behavior hardcoded in arrayfire for real floats
        ST = S <: FloatingPoint ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(a, y))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(y::S, x::AFAbstractArray{T})
        a = x.array
        # This is special behavior hardcoded in arrayfire for real floats
        ST = S <: FloatingPoint ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(y, a))
    end
end
# TODO: add! using +=, etc.

import Base: abs, min, max

function abs(a::AFAbstractArray)
    b = AFArray()
    icxx"$b = af::abs($a);"
    AFArray{backend_eltype(b)}(b)
end

#Max
function max(a::AFAbstractArray, val::Real)
    out = AFArray();
    icxx"""
        float val = $val;
        $out = af::max($a, val);
    """
    AFArray{backend_eltype(out)}(out)
end
max(val::Real, a::AFAbstractArray) = max(a, val)

#Min
function min(a::AFAbstractArray, val::Real)
    out = AFArray();
    icxx"""
        float val = $val;
        $out = af::min($a, val);
    """
    AFArray{backend_eltype(out)}(out)
end
min(val::Real, a::AFAbstractArray) = min(a, val)

#Negation
-{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"-$a;")

#Logical ops
import Base: ==, .==, .>, .<, .>=, .<=

#Equals
==(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a == $b;"
.==(a::AFAbstractArray, val::Real) = AFArray{Bool}(icxx"$a == $val;")

#Greater
.>(a::AFAbstractArray, b::AFAbstractArray) = AFArray{Bool}(icxx"$a > $b;")
.>(a::AFAbstractArray, b::Real) = AFArray{Bool}(icxx"$a > $b;")
.>=(a::AFAbstractArray, b::AFAbstractArray) = AFArray{Bool}(icxx"$a >= $b;")
.>=(a::AFAbstractArray, b::Real) = AFArray{Bool}(icxx"$a >= $b;")

#Lesser
.<(a::AFAbstractArray, b::AFAbstractArray) = AFArray{Bool}(icxx"$a < $b;")
.<(a::AFAbstractArray, b::Real) = AFArray{Bool}(icxx"$a < $b;")
.<=(a::AFAbstractArray, b::AFAbstractArray) = AFArray{Bool}(icxx"$a <= $b;")
.<=(a::AFAbstractArray, b::Real) = AFArray{Bool}(icxx"$a <= $b")
