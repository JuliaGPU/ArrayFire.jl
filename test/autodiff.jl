using AutoDiffSource

function checkdiff_inferred(f, δf, x0...)
    x = [x0...]
    y0 = f(x...)
    @assert length(y0) == 1 "Scalar functions only"
    y, ∇f = Test.@inferred δf(x...)
    @assert isapprox(y0, y) "Return values do not match"
    @assert typeof(y0) === typeof(y) "Return type doesn't match"
    ∂x = @inferred ∇f(1.0)
    if isa(∂x, Tuple)
        @assert typeof(x0) === typeof(∂x) "Gradient type doesn't match: $(typeof(x0)) vs $(typeof(∂x))"
    else
        @assert typeof(x0) === typeof((∂x,)) "Gradient type doesn't match: : $(typeof(x0)) vs $(typeof((∂x,)))"
    end
    checkgrad(f, x, ∂x)
end

typealias AFVector AFArray{Float64,1}
typealias AFMatrix AFArray{Float64,2}
safediv{T<:AbstractFloat}(x::T, y)::T = y == 0 ? 0::T : x / y
safediv(x::AbstractArray, y) = safediv.(x, y)
safediv(x::AFArray, y) = x ./ y
δdot_power(x::AbstractArray, y::AbstractFloat) = (t = x.^y; (t, z->(safediv(z.*y.*t, x), sum(z.*t.*log.(x)))))
δdot_power(x::AbstractMatrix, y::AbstractVector) = (t = x.^y; (t, z->(safediv(z.*y.*t, x), vec(sum(z.*t.*log(x), 2)))))
δdot_power(x::AbstractVector, y::AbstractMatrix) = (t = x.^y; (t, z->(vec(sum(safediv(z.*y.*t, x), 2)), z.*t.*log(x))))
δdot_power_const2(x, y) = (t = x.^y; (t, z-> y == 2 ? z.*2x : safediv(z.*y.*t, x)))
δdot_power{T}(x::T, y::T) = (t = x.^y; (t, z->(safediv(z.*y.*t, x), z.*t.*log.(x))))
δsum(x::AbstractArray) = (t = size(x); (sum(x), z->fill(z, t)))
δsum(x::AFArray) = (t = size(x); (sum(x), z->AFArray(fill(z, t))))


@test checkdiff_inferred(sum, δsum, AFVector(rand(3)))
@test checkdiff_inferred(sum, δsum, AFMatrix(rand(3, 2)))
