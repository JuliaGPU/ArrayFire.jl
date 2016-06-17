using Base.Meta

import Base: complex, conj, real, imag, max, min, abs, round, sign, floor, hypot
import Base: &, |, $, .>, .>=, .<, .<=, !, .==, .!=, ^, .^, /, ./

export sigmoid

for (op,fn) in ((:+,:af_add), (:.+,:af_add), (:-,:af_sub), (:.-,:af_sub), (:.*,:af_mul), 
                (:./,:af_div), (:%, :af_mod), (:.%, :af_mod), (:<<, :af_bitshiftl),
                (:.<<, :af_bitshiftl),(:>>, :af_bitshiftr), (:.>>, :af_bitshiftr),
                (:^, :af_pow), (:.^, :af_pow ))

    @eval function Base.($(quot(op))){T,S}(a::AFArray{T}, b::AFArray{S}; batched = true)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, batched)
        AFArray{af_promote(T,S)}(ptr[]) 
    end

end

Base.(:-)(a::AFArray) = 0 - a
Base.(:+)(a::AFArray{Bool}, v::Bool) = +(a,Int(v))
Base.(:+)(v::Bool, a::AFArray{Bool}) = +(a,v)
Base.(:-)(a::AFArray{Bool}, v::Bool) = -(a,Int(v))
Base.(:-)(v::Bool, a::AFArray{Bool}) = -(a,v)

^{T<:Real}(a::AFArray{T}, v::Integer) = ^(a, Real(v))
^{T<:Real}(a::AFArray{Complex{T}}, v::Integer) = ^(a, Real(v))
.^{T<:Real}(v::Irrational{:e}, a::AFArray{T}) = ^(Real(e), a) 
.^{T<:Real}(v::Irrational{:e}, a::AFArray{Complex{T}}) = ^(Real(e), a) 

for (op,fn) in ((:+, :af_add), (:.+, :af_add), (:-, :af_sub), (:.-, :af_sub), (:*, :af_mul), 
                (:.*, :af_mul), (:./, :af_div), (:%, :af_mod), (:.%, :af_mod),
                (:^, :af_pow), (:.^, :af_pow))

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(a::AFArray{T}, v::S)
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{af_promote(T,S)}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(a::AFArray{T}, v::Complex{S})
        b = constant(Complex{af_promote(T,S)}(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(v::S, a::AFArray{T})
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, b, a, true)
        AFArray{af_promote(T,S)}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(v::Complex{S}, a::AFArray{T})
        b = constant(Complex{af_promote(T,S)}(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, b, a, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(a::AFArray{Complex{T}}, v::S)
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(a::AFArray{Complex{T}}, v::Complex{S})
        b = constant(Complex{af_promote(T,S)}(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(v::S, a::AFArray{Complex{T}})
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, b, a, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(v::Complex{S}, a::AFArray{Complex{T}})
        b = constant(Complex{af_promote(T,S)}(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, b, a, true)
        AFArray{Complex{af_promote(T,S)}}(ptr[]) 
    end
end
/{T,S<:Real}(a::AFArray{T}, v::S) = ./(a, v)
/{T,S<:Real}(a::AFArray{T}, v::Complex{S}) = ./(a, v)

for (op,fn) in ((:sin, :af_sin), (:cos, :af_cos), (:tan, :af_tan), (:asin, :af_asin), 
                (:acos, :af_acos), (:atan, :af_atan), (:sinh, :af_sinh), (:cosh, :af_cosh),
                (:tanh, :af_tanh), (:asinh, :af_asinh), (:acosh, :af_acosh), (:atanh, :af_atanh),
                (:cbrt, :af_cbrt), (:erf, :af_erf), (:erfc, :af_erfc), (:exp, :af_exp),
                (:expm1, :af_expm1), (:factorial, :af_factorial), (:lgamma, :af_lgamma),
                (:log, :af_log), (:log10, :af_log10), (:log1p, :af_log1p), (:sqrt, :af_sqrt),
                (:gamma, :af_tgamma), (:log2, :af_log2))
    @eval function Base.($(quot(op)))(a::AFArray) 
        out = new_ptr()
        eval($(quot(fn)))(out, a)
        AFArray{backend_eltype(out[])}(out[])
    end

end

for (op,fn) in ((:atan2, :af_atan2),)
    @eval function Base.($(quot(op))){T,S}(a::AFArray{T}, b::AFArray{S}; batched = true)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, batched)
        AFArray{af_promote(T,S)}(ptr[])
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(a::AFArray{T}, v::S)
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{af_promote(T,S)}(ptr[])
    end

    @eval function Base.($(quot(op))){T<:Real,S<:Real}(v::S, a::AFArray{T})
        b = constant((af_promote(T,S))(v), size(a)...)
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, b, a, true)
        AFArray{af_promote(T,S)}(ptr[])
    end
end

function sigmoid(a::AFArray)
    out = new_ptr()
    af_sigmoid(out, a)
    AFArray{backend_eltype(out[])}(out[])
end

# Complex functions

function complex{T<:Real}(a::AFArray{T})
    out = new_ptr()
    af_cplx(out, a)
    AFArray{Complex{T}}(out[])
end

function complex{T<:Real}(a::AFArray{T}, b::AFArray{T}; batch = true)
    out = new_ptr()
    af_cplx2(out, a, b, batch)
    AFArray{Complex{T}}(out[])
end

function conj{T<:Complex}(a::AFArray{T})
    out = new_ptr()
    af_conjg(out, a)
    AFArray{T}(out[])
end

function real{T<:Real}(a::AFArray{Complex{T}})
    out = new_ptr()
    af_real(out, a)
    AFArray{T}(out[])
end 

function real{T<:Real}(a::AFArray{T})
    out = new_ptr()
    af_real(out, a)
    AFArray{T}(out[])
end 

function imag{T}(a::AFArray{T})
    out = new_ptr()
    af_imag(out, a)
    AFArray{T}(out[])
end

# Numeric functions

for (op,fn) in ((:max, :af_maxof), (:min, :af_minof))

    @eval function ($op){T,S}(a::AFArray{T}, b::AFArray{S}; batched = true)
        out = new_ptr()
        eval($fn)(out, a, b, true)
        AFArray{af_promote(T,S)}(out[])
    end

    @eval function ($op){T,S<:Real}(a::AFArray{T}, b::S)
        out = new_ptr()
        tmp = constant(b, size(a))
        eval($fn)(out, a, tmp, true)
        AFArray{af_promote(T,S)}(out[])
    end
    @eval ($op){T,S<:Real}(b::S, a::AFArray{T}) = max(a, b)

    @eval function ($op){T<:Complex,S<:Complex}(a::AFArray{T}, b::S)
        out = new_ptr()
        tmp = constant(b, size(a))
        eval($fn)(out, a, tmp, true)
        AFArray{af_promote(T,S)}(out[])
    end
    @eval ($op){T<:Complex,S<:Complex}(b::S, a::AFArray{T}) = max(a, b)

end

for (op,fn) in ((:abs, :af_abs), (:arg, :af_arg),
                (:ceil, :af_ceil), (:sign, :af_sign),
                (:floor, :af_floor), (:round, :af_round), 
                (:trunc, :af_trunc))

    @eval function ($op){T}(a::AFArray{T})
        out = new_ptr()
        eval($fn)(out, a)
        AFArray{T}(out[])
    end

end

function hypot{T}(a::AFArray{T}, b::AFArray{T}; batched = true)
    out = new_ptr()
    af_hypot(out, a, b, batched)
    AFArray{T}(out[])
end

# Logical Operations

for (op,fn) in ((:&, :af_bitand),(:|, :af_bitor), (:$, :af_bitxor),
                (:.==, :af_eq), (:.!=, :af_neq), (:.>, :af_gt), 
                (:.>=, :af_ge), (:.<, :af_lt), (:.<=, :af_le),
                (:.!=, :af_neq))

    @eval function ($op)(a::AFArray, b::AFArray; batched = true)
        out = new_ptr()
        eval($fn)(out, a, b, batched)
        AFArray{backend_eltype(out[])}(out[])
    end
    
    @eval function ($op)(a::AFArray, b::Real; batched = true)
        out = new_ptr()
        tmp = constant(b, size(a))
        eval($fn)(out, a, tmp, batched)
        AFArray{backend_eltype(out[])}(out[])
    end
    @eval ($op)(b::Real, a::AFArray) = ($op)(a, b)

end

function !{T}(a::AFArray{T})
    out = new_ptr()
    af_not(out, a)
    AFArray{T}(out[]) 
end
