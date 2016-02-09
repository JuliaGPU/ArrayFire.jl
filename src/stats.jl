import Base: mean
    
#Mean
#TODO : generate mean of all elements
mean{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::mean($a, $(dim-1));")
mean{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray(icxx"af::mean($a, $weights, $dim);")
