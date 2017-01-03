type AFArray{T,N}
    arr::af_array
    function AFArray(arr::af_array)
        a = new(arr)
        finalizer(a, x -> af_release_array(x))
        a
    end
end

typealias AFVector{T} AFArray{T,1}
typealias AFMatrix{T} AFArray{T,2}
typealias AFVolume{T} AFArray{T,3}
typealias AFTensor{T} AFArray{T,4}

export AFArray, AFVector, AFMatrix, AFVolume, AFTensor

import Base: convert, copy

convert{T,N}(::Type{AFArray{T,N}}, a::AbstractArray{T,N}) = af_create_array(a)
copy{T,N}(a::AFArray{T,N}) = af_copy_array(a)
