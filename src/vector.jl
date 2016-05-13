### Vector Algorithms

# Reduction 

import Base: sum, min, max, minimum, maximum

for (op,fn) in ((:sum, :af_sum_all), (:product, :af_product_all),
                (:maximum, :af_max_all), (:minimum, af_min_all))

    @eval function ($op){T<:Real}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($(quot(fn)))(real, imag, a)
        real[]
    end

    @eval function ($op){T<:Complex}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($(quot(fn)))(real, imag, a)
        complex(real[], imag[])
    end

end

for (op, fn) in ((:sum, :af_sum), (:product, :af_product), 
                (:max, :af_max), (:min, :af_min))

    @eval function ($op){T}(a::AFArray{T}, dim::Integer)
        dim = dim - 1
        out = new_ptr()
        eval($(quot(fn)))(out, a, dim)
        AFArray{T}(out[])
    end

end
