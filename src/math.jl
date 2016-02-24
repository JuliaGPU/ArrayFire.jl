# MATH OPERATIONS

for op in (:sin, :cos, :tan, :atan, :asin, :acos, :log, :log1p, :log10, :sqrt, :transpose,
    :exp, :expm1, :erf, :erfc, :cbrt, :lgamma, :transpose, :acosh, :cosh, :asinh, :sinh, :atanh, :tanh, :factorial)
    @eval Base.($(quot(op))){T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::($op)(x.array))
end

Base.gamma{T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::tgamma(x))

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
        ST = S <: AbstractFloat ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(a, y))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(y::S, x::AFAbstractArray{T})
        a = x.array
        # This is special behavior hardcoded in arrayfire for real floats
        ST = S <: AbstractFloat ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(y, a))
    end
end
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
.==(a::AFAbstractArray, b::Union{AFAbstractArray, Real}) = af_equals(a,b)
.==(a::Union{AFAbstractArray,Real}, b::AFAbstractArray) = af_equals(a,b)

#Greater
.>(a::AFAbstractArray, b::Union{AFAbstractArray,Real}) = af_gt(a,b)
.>(a::Union{Real, AFAbstractArray}, b::AFAbstractArray) = af_gt(a,b)
.>=(a::AFAbstractArray, b::Union{AFAbstractArray, Real}) = af_ge(a,b)
.>=(a::Union{Real,AFAbstractArray}, b::AFAbstractArray) =  af_ge(a,b)

#Lesser
.<(a::AFAbstractArray, b::Union{AFAbstractArray, Real}) = af_lt(a,b)
.<(a::Union{AFAbstractArray,Real}, b::AFAbstractArray) = af_lt(a,b)
.<=(a::AFAbstractArray, b::Union{AFAbstractArray, Real}) = af_le(a,b)
.<=(a::Union{AFAbstractArray, Real}, b::AFAbstractArray) = af_le(a,b)

#Or
|(a::AFAbstractArray, b::AFAbstractArray) = af_bitor(a,b)
