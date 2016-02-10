import Base: mean, std, median, var, cov
export average, corrcoef
    
#Mean
#TODO Import use Base.mean for this instead of a different name
average{T}(a::AFAbstractArray{T}) = icxx"af::mean<float>($a);"
mean{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::mean($a, $(dim-1));")
mean{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray(icxx"af::mean($a, $weights, $dim);")

#Median
median{T}(a::AFAbstractArray{T}) = icxx"af::median<float>($a);"
median{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::median($a, $(dim-1));")

#Standard Deviation
std{T}(a::AFAbstractArray{T}) = icxx"af::stdev<float>($a);"
std{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::stdev($a, $(dim-1));")

#Variance
var{T}(a::AFAbstractArray{T}) = icxx"af::var<float>($a);"
var{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::var($a, false, $(dim-1));")
var{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray(icxx"af::var($a, $weights, $dim);")

#Covariance
cov{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(icxx"af::cov($a, $b);")

#Correlation coefficient
corrcoef{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = icxx"af::corrcoef($a, $b);"
