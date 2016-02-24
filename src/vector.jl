#####################
# Vector operations #
#####################


#Reduction operations

import Base: sum, maximum, minimum, any, max, min

#Sum 
sum(a::AFAbstractArray) = af_sum(a)
sum{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_sum(a,dim-1))

#Product
product(a::AFAbstractArray) = af_product(a)
product{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_product(a,dim-1))

#Maximum
maximum(a::AFAbstractArray) = af_max(a)
maximum{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_max_dim(a,dim-1))

#Minimum
minimum(a::AFAbstractArray) = af_min(a)
minimum{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_min_dim(a,dim-1))

#Any
any(a::AFAbstractArray) = af_anyTrue(a)
any(a::AFAbstractArray, dim::Integer) = AFArray{Bool}(af_anyTrue(a, dim-1))

#Alltrue
alltrue(a::AFAbstractArray) = af_allTrue(a)
alltrue(a::AFAbstractArray, dim::Integer) = AFArray{Bool}(icxx"af::allTrue($a, $(dim - 1));")

#Count 
countnz(a::AFAbstractArray) = af_count(a)
countnz(a::AFAbstractArray, dim::Integer) = AFArray{UInt32}(af_count(a, dim-1))


#Inclusive Scan Operations

import Base: cumsum, find

#Cumsum
cumsum{T}(a::AFAbstractArray{T}, dim::Integer = 1) = AFArray{T}(af_accum(a,dim-1))

#Find
find(a::AFAbstractArray) = AFArray{UInt32}(af_where(a))


#Numerical differentiation

import Base: diff, gradient

#1D Diff
diff{T}(a::AFAbstractArray{T}, dim::Integer = 1) = AFArray{T}(af_diff(a,dim-1))

#2D Siff
diff2{T}(a::AFAbstractArray{T}, dim::Integer = 1) = AFArray{T}(af_diff2(a,dim-1))

#Gradient
function gradient{T}(a::AFAbstractArray{T})
    dx = AFArray()
    dy = AFArray()
    af_grad(dx, dy, a)
    AFArray{T}(dx), AFArray{T}(dy)
end


#Sorting

import Base: sort
export sortIndex

#Sort
#NOTE: ArrayFire currently restricts dim to the value 0, which means you can sort only columns
function sort{T}(A::AFArray{T}; rev = false)
    ndims(A) == 1 ||
        error("Must explicitly specify dimension when sorting mutlidimensional array")
    AFArray{T}(af_sort(A,0,!rev))
end
function sort{T}(A::AFArray{T}, dim; rev = false) 
    if dim > 1
        error("ArrayFire currently lets you sort only by columns")
    else
        AFArray{T}(af_sort(A,dim-1,!rev))
    end
end

#sortIndex
function sortIndex(a::AFArray, dim = 1; rev = false)
    (ndims(a) == 2 && dim == 2) && error("ArrayFire currently supports sorting along columns")
    out  = AFArray()
    ind = AFArray()
    af_sort(out, ind, a, dim-1, !rev)
    AFArray{backend_eltype(out)}(out), AFArray{backend_eltype(ind)}(ind)
end

#sortByKey
function sortByKey(keys::AFArray, vals::AFArray, dim = 1; rev = false)
    (ndims(a) == 2 && dim == 2) && error("ArrayFire currently supports sorting along columns")
    out_keys = AFArray()
    out_vals = AFArray()
    af_sort(out_keys, out_vals, keys, vals, dim-1, !rev)
    AFArray{backend_eltype(out_keys)}(out_keys), AFArray{backend_eltype(out_vals)}(out_vals)
end

#Set Operations

import Base: setdiff, union, unique

#Setdiff
setdiff{T}(a::AFAbstractArray{T} , b::AFAbstractArray{T}) = AFArray{T}(af_setIntersect(a,b))

#Union
union{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_setUnion(a,b))

#Unique
unique{T}(a::AFAbstractArray{T}) = AFArray{T}(af_setUnique(a))


#Reordering functions

import Base:flipdim, vec

flipdim{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_flip(a,dim-1))
vec{T}(a::AFAbstractArray{T}) = AFArray{T}(af_flat(a))
