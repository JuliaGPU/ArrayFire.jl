type AFArray{T,N}
    arr::af_array
    function AFArray(arr::af_array)
        a = new(arr)
        finalizer(a, x -> af_release_array(x))
        @assert af_get_type!(arr) == T
        @assert af_get_numdims!(arr) == N
        a
    end
end

typealias AFVector{T} AFArray{T,1}
typealias AFMatrix{T} AFArray{T,2}
typealias AFVolume{T} AFArray{T,3}
typealias AFTensor{T} AFArray{T,4}

export AFArray, AFVector, AFMatrix, AFVolume, AFTensor

import Base: convert, copy, deepcopy_internal

convert{T,N}(::Type{AFArray{T,N}}, a::AbstractArray{T,N}) = af_create_array(a)
copy{T,N}(a::AFArray{T,N}) = af_copy_array(a)
deepcopy_internal{T,N}(a::AFArray{T,N}, d::ObjectIdDict) = haskey(d, a) ? d[a]::AFArray{T,N} : copy(a)
