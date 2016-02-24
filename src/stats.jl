import Base: mean, std, median, var, cov
export corrcoef
    
#Mean
mean(a::AFAbstractArray) = af_average(a)
mean{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_mean(a, dim-1))
mean{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_mean(a,weights,dim-1))

#Median
median(a::AFAbstractArray) = af_median(a)
median{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_median(a, dim-1))

#Standard Deviation
std(a::AFAbstractArray) = af_stdev(a)
std{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_stdev(a,dim-1))

#Variance
var(a::AFAbstractArray) = af_var(a)
var{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_var(a,dim-1))
var{T<:Real}(a::AFAbstractArray{T}, weights::AFAbstractArray{T}, dim::Integer) = AFArray{T}(af_var(a,weights,dim-1))

#Covariance
cov{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_cov(a,b))

#Correlation coefficient
corrcoef{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = af_corrcoef(a,b)
