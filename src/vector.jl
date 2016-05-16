### Vector Algorithms

import Base: sum, min, max, minimum, maximum, countnz, any, all, sort

export sortIndex

# Reduction 

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

function countnz{T}(a::AFArray{T})
    real = Base.Ref{Cdouble}(0)
    imag = Base.Ref{Cdouble}(0)
    af_count_all(real, imag, a)
    Int(real[])
end

function countnz{T}(a::AFArray{T}, dim::Integer)
    dim = dim - 1
    out = new_ptr()
    af_count(out, a, dim)
    AFArray{Cuint}(out[])
end

for (op,fn) in ((:any, :af_any_true_all),(:all, :af_all_true_all))

    @eval function ($op){T}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        eval($fn)(real, imag, a)
        Bool(real[])
    end

end

for (op,fn) in ((:any, :af_any_true), (:all, :af_all_true))
 
    @eval function ($op){T}(a::AFArray{T}, dim::Integer)
        dim = dim - 1
        out = new_ptr()
        eval($fn)(out, a, dim)
        AFArray{Bool}(out[])
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

# Sorting 

function sort{T}(a::AFArray{T}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out = new_ptr()
    af_sort(out, a, Cuint(dim-1), !rev)
    AFArray{T}(out[])
end

function sortIndex{T}(a::AFArray{T}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out = new_ptr()
    indices = new_ptr()
    af_sort_index(out, indices, a, Cuint(dim-1), !rev)
    AFArray{T}(out[]), AFArray{Int32}(indices[]) + 1
end

function sortByKey{S,T}(keys::AFArray{S}, values::AFArray{T}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out_keys = new_ptr()
    out_values = new_ptr()
    af_sort_by_key(out_keys, out_values, keys, values, Cuint(dim - 1), !rev)
    AFArray{S}(out_keys[]), AFArray{T}(out_values[])
end
