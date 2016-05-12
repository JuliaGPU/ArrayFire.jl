using Base.Meta

for (op,fn) in ((:+,:af_add), (:.+,:af_add), (:-,:af_sub), (:.-,:af_sub), (:.*,:af_mul), (:./,:af_div))
    @eval function Base.($(quot(op))){T,S}(a::AFArray{T}, b::AFArray{S})
        ptr = new_ptr()
        eval($(quot(fn)))(ptr, a, b, true)
        AFArray{af_promote(T,S)}(ptr[]) 
    end
end

Base.(:+)(a::AFArray{Bool}, v::Bool) = +(a,Int(v))
Base.(:+)(v::Bool, a::AFArray{Bool}) = +(a,v)
Base.(:-)(a::AFArray{Bool}, v::Bool) = -(a,Int(v))
Base.(:-)(v::Bool, a::AFArray{Bool}) = -(a,v)

for op in (:+, :.+, :-, :.-, :*, :.*, :/, :./)
    @eval function Base.($(quot(op))){T<:Number,S<:Number}(a::AFArray{T}, v::S)
        b = constant(v, size(a)...)
        eval($(quot(op)))(a, b)
    end
    @eval Base.($(quot(op))){S<:Number,T<:Number}(v::S, a::AFArray{T}) = Base.($(quot(op)))(a, v)
end

for (op,fn) in ((:sin, :af_sin), (:cos, :af_cos), (:tan, :af_tan), (:asin, :af_asin), 
                (:acos, :af_acos), (:atan, :af_atan), (:sinh, :af_sinh), (:cosh, :af_cosh),
                (:tanh, :af_tanh), (:asinh, :af_asinh), (:acosh, :af_acosh), (:atanh, :af_atanh))
    @eval function Base.($(quot(op)))(a::AFArray) 
        out = new_ptr()
        eval($(quot(fn)))(out, a)
        AFArray{backend_eltype(out[])}(out[])
    end
end
