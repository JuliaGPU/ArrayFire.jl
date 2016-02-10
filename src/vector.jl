#####################
# Vector operations #
#####################


#Reduction operations

import Base: sum, max, min, any 
export product, alltrue

#Sum 
sum{T}(a::AFAbstractArray{T}) = icxx"af::sum<float>($a);"
sum{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::sum($a, $(dim -1);")

#Product
product{T}(a::AFAbstractArray{T}) = icxx"af::product<float>($a);"
product{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::product($a, $(dim -1);")

#Max
max{T}(a::AFAbstractArray{T}) = icxx"af::max<float>($a);"
max{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::max($a, $(dim -1);")

#Min
min{T}(a::AFAbstractArray{T}) = icxx"af::min<float>($a);"
min{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{T}(icxx"af::min($a, $(dim -1);")

#Any
any{T}(a::AFAbstractArray{T}) = icxx"af::anyTrue<boolean_t>($a);" != 0
any(a::AFAbstractArray, dim::Integer) = AFArray{Bool}(icxx"af::anyTrue($a, $(dim -1);")

#Alltrue
alltrue{T}(a::AFAbstractArray{T}) = icxx"af::allTrue<boolean_t>($a);" != 0
alltrue(a::AFAbstractArray{T}, dim::Integer) = AFArray{Bool}(icxx"af::allTrue($a, $(dim -1);")

#Count 
countnz{T}(a::AFAbstractArray{T}) = icxx"af::count<float>($a);"
countnz{T}(a::AFAbstractArray{T}, dim::Integer) = AFArray{UInt32}(icxx"af::count($a, $(dim -1);")
