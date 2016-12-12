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
δsum(x::AFArray) = (t = size(x); (sum(x), z->constant(z, t)))

rnd(len) = rand(AFVector, len)
rnd(len1, len2) = rand(AFMatrix, (len1, len2))

@test checkdiff_inferred(sum, δsum, rnd(3))
@test checkdiff_inferred(sum, δsum, rnd(3, 2))
