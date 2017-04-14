@testset "Indexing" begin
a1 = AFArray(rand(3, 2))
b1 = AFArray([1. 2.])
a1[1, :] = b1
@test all(@inferred(a1[1, :]) == b1)
a1[:, :] = 0
@test all(a1 == 0)
@test typeof(a1[1,1]) == Float64
@test @inferred(a1[1,1]) == 0
end
