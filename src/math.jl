# MATH OPERATIONS

for op in (:sin, :cos, :tan, :atan, :asin, :acos, :log, :log1p, :log10, :sqrt, 
    :exp, :expm1, :erf, :erfc, :cbrt, :lgamma, :acosh, :cosh, :asinh, :sinh, :atanh, :tanh, :factorial)
    @eval Base.($(quot(op))){T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::($op)(x.array))
end

Base.gamma{T}(x::AFAbstractArray{T}) = AFArray{T}(@cxx af::tgamma(x))

import Base: +, -, abs, .^ 

# TODO: add! using +=, etc.

.^{T,S}(a::AFAbstractArray{T}, b::AFAbstractArray{S}) = AFArray{af_promote(T,S)}(af_pow(a,b))
.^{T}(a::Real, b::AFAbstractArray{T}) = AFArray{T}(af_pow(a,b))
.^{T}(a::AFAbstractArray{T}, b::Real) = AFArray{T}(af_pow(a,b))

export root
root{T,S}(a::AFAbstractArray{T}, b::AFAbstractArray{S}) = AFArray{af_promote(T,S)}(af_root(a,b))
root{T}(a::Real, b::AFAbstractArray{T}) = AFArray{T}(af_root(a,b))
root{T}(a::AFAbstractArray{T}, b::Real) = AFArray{T}(af_root(a,b))

import Base: abs, min, max

#Abs
abs{T}(a::AFAbstractArray{T}) = AFArray{T}(af_abs(a))

#Max
function max(a::AFAbstractArray, b::AFAbstractArray)
    out = af_max(a,b)   
    AFArray{backend_eltype(out)}(out)
end
function max(a::AFAbstractArray, b::Real)
    out = af_max(a,b)
    AFArray{backend_eltype(out)}(out)
end
max(a::Real, b::AFAbstractArray) = max(b, a)

#Min
function min(a::AFAbstractArray, b::AFAbstractArray)
    out = af_min(a, b)
    AFArray{backend_eltype(out)}(out)
end
function min(a::AFAbstractArray, b::Real)
    out = af_min(a,b)
    AFArray{backend_eltype(out)}(out)
end
min(a::Real, b::AFAbstractArray) = max(b, a)

#Negation
-{T}(a::AFSubArray{T}) = AFArray{T}(icxx"0-$a;")
-{T}(a::AFArray{T}) = AFArray{T}(af_neg(a))

#Logical ops
import Base: ==, .==, .>, .<, .>=, .<=, |

#Equals
.==(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = AFArray{Bool}(af_equals(a,b))

#Greater
.>(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = AFArray{Bool}(af_gt(a,b))
.>=(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = AFArray{Bool}(af_ge(a,b))

#Lesser
.<(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = AFArray{Bool}(af_lt(a,b))
.<=(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = AFArray{Bool}(af_le(a,b))

#Or
|(a::AFAbstractArray{Bool}, b::AFAbstractArray{Bool}) = AFArray{Bool}(af_bitor(a,b))

#Complex functions
import Base:complex, real, imag

complex{T}(a::AFAbstractArray{T}) = AFArray{Complex{T}}(af_complex(a))
complex{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{Complex{T}}(af_complex(a,b))
real{T}(a::AFAbstractArray{Complex{T}}) = AFArray{T}(af_real(a))
imag{T}(a::AFAbstractArray{Complex{T}}) = AFArray{T}(af_imag(a))
