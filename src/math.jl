# MATH OPERATIONS

for op in (:sin, :cos, :tan, :atan, :asin, :acos, :log, :log1p, :log10, :sqrt, 
    :exp, :expm1, :erf, :erfc, :cbrt, :lgamma, :transpose, :acosh, :cosh, :asinh, :sinh, :atanh, :tanh, :factorial)
    @eval Base.($(quot(op))){T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::($op)(x.array))
end

Base.gamma{T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::tgamma(x))

import Base: +, -, abs

# TODO: add! using +=, etc.

import Base: abs, min, max

#Abs
abs{T}(a::AFAbstractArray{T}) = AFArray{T}(af_abs(a))

#Max
function max(a::AFAbstractArray, val::Real)
    out = af_max(a, val)
    AFArray{backend_eltype(out)}(out)
end
max(val::Real, a::AFAbstractArray) = max(a, val)

#Min
function min(a::AFAbstractArray, val::Real)
    out = af_min(a, val)
    AFArray{backend_eltype(out)}(out)
end
min(val::Real, a::AFAbstractArray) = min(a, val)

#Negation
-{T}(a::AFSubArray{T}) = AFArray{T}(icxx"0-$a;")
-{T}(a::AFArray{T}) = AFArray{T}(af_neg(a))

#Logical ops
import Base: ==, .==, .>, .<, .>=, .<=

#Equals
.==(a::AFAbstractArray, b::AFAbstractArray) = af_equals(a,b)
.==(a::AFAbstractArray, b::Real) = af_equals(a,b)
.==(a::Real, b::AFAbstractArray) = af_equals(a,b)

#Greater
.>(a::AFAbstractArray, b::AFAbstractArray) = af_gt(a,b)
.>(a::Real, b::AFAbstractArray) = af_gt(a,b)
.>(a::AFAbstractArray, b::Real) = af_gt(a,b)
.>=(a::AFAbstractArray, b::AFAbstractArray) = af_ge(a,b)
.>=(a::Real, b::AFAbstractArray) =  af_ge(a,b)
.>=(a::AFAbstractArray, b::Real) =  af_ge(a,b)

#Lesser
.<(a::AFAbstractArray, b::AFAbstractArray) = af_lt(a,b)
.<(a::Real, b::AFAbstractArray) = af_lt(a,b)
.<(a::AFAbstractArray, b::Real) = af_lt(a,b)
.<=(a::AFAbstractArray, b::AFAbstractArray) = af_le(a,b)
.<=(a::Real, b::AFAbstractArray) = af_le(a,b)
.<=(a::AFAbstractArray, b::Real) = af_le(a,b)

#Or
|(a::AFAbstractArray, b::AFAbstractArray) = af_bitor(a,b)
