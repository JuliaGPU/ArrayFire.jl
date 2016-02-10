import Base: mean
export average
    
#Mean
#TODO Import use Base.mean for this instead of a different name
average{T}(a::AFAbstractArray{T}) = icxx"af::mean<float>($a);"
mean{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::mean($a, $(dim-1));")
mean{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray(icxx"af::mean($a, $weights, $dim);")
