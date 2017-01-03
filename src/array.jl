type AFArray{T,N}
    arr::af_array
    function AFArray(arr::af_array)
        a = new(arr)
        finalizer(a, x -> release_array(x))
        @assert get_type!(arr) == T "type mismatch: $(af_get_type!(arr)) != $T"
        @assert get_numdims!(arr) == N "dims mismatch: $(af_get_numdims!(arr)) != $N"
        a
    end
end

typealias AFVector{T} AFArray{T,1}
typealias AFMatrix{T} AFArray{T,2}
typealias AFVolume{T} AFArray{T,3}
typealias AFTensor{T} AFArray{T,4}

export AFArray, AFVector, AFMatrix, AFVolume, AFTensor

import Base: convert, copy, deepcopy_internal

convert{T,N}(::Type{AFArray{T,N}}, a::AbstractArray{T,N}) = create_array(a)
deepcopy_internal{T,N}(a::AFArray{T,N}, d::ObjectIdDict) = haskey(d, a) ? d[a]::AFArray{T,N} : copy(a)

import Base: abs, acos, acosh, asin, asinh, atan, atan2, atanh, cbrt, ceil, clamp, cos, cosh
import Base: count, cov, det, div, dot, erf, erfc, exp, expm1, factorial, fft, floor, gradient, hypot
import Base: identity, ifft, imag, isinf, isnan, join, lgamma, log, log10, log1p, log2, lu, maximum, mean, median
import Base: minimum, mod, norm, prod, qr, randn, range, rank, real, rem, replace, round, scale, select, show
import Base: signbit, sin, sinh, sort, sqrt, sub, sum, svd, tan, tanh, transpose, trunc, var
