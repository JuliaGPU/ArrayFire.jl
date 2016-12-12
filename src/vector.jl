### Vector Algorithms

import Base: sum, min, max, minimum, maximum, countnz, any, all, sort,
                union, find, cumsum, diff, findmax, findmin, sortperm

export sortIndex, sortByKey, diff2, minidx, maxidx

# Reduction 

for (op,fn) in ((:sum, :af_sum_all), (:product, :af_product_all),
                (:maximum, :af_max_all), (:minimum, af_min_all))

    @eval function ($op){T<:Real}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        $(fn)(real, imag, a)
        real[]
    end

    @eval function ($op){T<:Complex}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        $(fn)(real, imag, a)
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
        $(fn)(real, imag, a)
        Bool(real[])
    end

end

for (op,fn) in ((:any, :af_any_true), (:all, :af_all_true))
 
    @eval function ($op){T}(a::AFArray{T}, dim::Integer)
        dim = dim - 1
        out = new_ptr()
        $(fn)(out, a, dim)
        AFArray{Bool}(out[])
    end

end

for (op, fn) in ((:sum, :af_sum), (:product, :af_product), 
                (:maximum, :af_max), (:minimum, :af_min))

    @eval function ($op){T}(a::AFArray{T}, dim::Integer)
        dim = dim - 1
        out = new_ptr()
        $(fn)(out, a, dim)
        AFArray{T}(out[])
    end

end

# Sorting 

function sort{T,N}(a::AFArray{T,N}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out = new_ptr()
    af_sort(out, a, Cuint(dim-1), !rev)
    AFArray{T,N}(out[])
end

function sortperm{T,N}(a::AFArray{T,N}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out = new_ptr()
    indices = new_ptr()
    af_sort_index(out, indices, a, Cuint(dim-1), !rev)
    AFArray{Int32,N}(indices[]) + 1
end

function sortByKey{S,T,N}(keys::AFArray{S,N}, values::AFArray{T,N}, dim::Integer = 1; rev = false)
    if dim == 2
        error("ArrayFire doesn't support sorting along this dimension")
    end
    out_keys = new_ptr()
    out_values = new_ptr()
    af_sort_by_key(out_keys, out_values, keys, values, Cuint(dim - 1), !rev)
    AFArray{S,N}(out_keys[]), AFArray{T,N}(out_values[])
end

# Set Operations

function union{T,S}(a::AFArray{T}, b::AFArray{S}; is_unique = false)
    out = new_ptr()
    af_set_union(out, a, b, is_unique)
    AFArray{af_promote(T,S)}(out[])
end

function unique{T}(a::AFArray{T}; is_sorted = false)
    out = new_ptr()
    af_set_unique(out, a, is_sorted)
    AFArray{T}(out[])
end

function setdiff{T,S}(a::AFArray{T}, b::AFArray{S}; is_unique = false)
    out = new_ptr()
    af_set_intersect(out, a, b, is_unique)
    AFArray{af_promote(T,S)}(out[])
end

# Inclusive Scan Operations

function cumsum{T}(a::AFArray{T}, dim::Integer = 1)
    out = new_ptr()
    af_accum(out, a, Cint(dim-1))
    AFArray{T}(out[])
end

function find(a::AFArray)
    out = new_ptr()
    af_where(out, a)
    AFArray{backend_eltype(out[])}(out[]) + 1
end

# Numerical Differentiation

function diff{T}(a::AFArray{T}, dim::Integer = 1)
    out = new_ptr()
    af_diff1(out, a, Cint(dim-1))
    AFArray{T}(out[])
end

function diff2{T}(a::AFArray{T}, dim::Integer = 1)
    out = new_ptr()
    af_diff2(out, a, Cint(dim-1))
    AFArray{T}(out[])
end

# Max, Min with Indices

for (op, fn) in ((:findmax, :af_imax_all), (:findmin, :af_imin_all))

    @eval function ($op){T<:Real}(a::AFArray{T})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        idx = Base.Ref{Cuint}(0)
        $(fn)(real, imag, idx, a)
        real[], Int(idx[]) + 1
    end

    @eval function ($op){T<:Real}(a::AFArray{Complex{T}})
        real = Base.Ref{Cdouble}(0)
        imag = Base.Ref{Cdouble}(0)
        idx = Base.Ref{Cuint}(0)
        $(fn)(real, imag, idx, a)
        complex(real[], imag[]), Int(idx[]) + 1
    end

end

for (op, fn) in ((:minidx, :af_imin), (:maxidx, :af_imax))
    
    @eval function ($op){T}(a::AFArray{T}, dim::Int)
        out = new_ptr()
        idx = new_ptr()
        $(fn)(out, idx, a, dim-1)
        AFArray{T}(out[]), 
        AFArray{backend_eltype(idx[])}(idx[]) + 1
    end

end
