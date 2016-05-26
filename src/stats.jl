# Statistics functions

import Base: mean, median, std, var, cov

export meanWeighted, varWeighted, corrcoef

for (op, fn) in ((:mean, :af_mean_all), (:median, :af_median_all),
                    (:std, :af_stdev_all))

    @eval function ($op){T<:Real}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($fn)(real, imag, a)
        T(real[])
    end 

    @eval function ($op){T<:Complex}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($fn)(real, imag, a)
        complex(real[], imag[])
    end

end

for (op, fn) in ((:meanWeighted, :af_mean_all_weighted), 
                    (:varWeighted, :af_var_all_weighted))

    @eval function ($op){T<:Real,S}(a::AFArray{T}, w::AFArray{S})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($fn)(real, imag, a, w)
        T(real[])
    end 

    @eval function ($op){T<:Complex,S}(a::AFArray{T}, w::AFArray{S})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($fn)(real, imag, a, w)
        complex(real[], imag[])
    end

end

for (op, fn) in ((:mean, :af_mean), (:median, :af_median), 
                    (:std, :af_stdev))

    @eval function ($op){T}(a::AFArray{T}, dim::Integer)
        out = new_ptr()
        eval($fn)(out, a, Cuint(dim-1))
        AFArray{T}(out[])
    end

end

for (op, fn) in ((:meanWeighted, :af_mean_weighted), (:varWeighted, :af_var_weighted))

    @eval function ($op){T}(a::AFArray{T}, w::AFArray{T}, dim::Integer)
        out = new_ptr()
        eval($fn)(out, a, w, Cuint(dim-1))
        AFArray{T}(out[])
    end

end

function var{T}(a::AFArray{T}; is_biased = false)
    real = Base.Ref{Cdouble}(0)
    imag = Base.Ref{Cdouble}(0)
    af_var_all(real, imag, a, is_biased)
    Float64(real[])
end

function var{T}(a::AFArray{T}, dim::Integer; is_biased = false)
    out = new_ptr()
    af_var(out, a, is_biased, Cuint(dim-1))
    AFArray{backend_eltype(out[])}(out[])
end

function cov(a::AFArray, b::AFArray; is_biased = false)
    out = new_ptr()
    af_cov(out, a, b, is_biased)
    AFArray{backend_eltype(out[])}(out[])
end

function corrcoef(a::AFArray, b::AFArray)
    real = Base.Ref{Cdouble}(0)
    imag = Base.Ref{Cdouble}(0)
    af_corrcoef(real, imag, a, b)
    if imag[] == 0
        real[]
    else
        complex(real[], imag[])
    end
end
    
