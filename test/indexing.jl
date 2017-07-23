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

@testset "Inplace assignment" begin
    a = randn(10)
    af = AFArray(a)
    c = zeros(a)
    cf = AFArray(c)

    c .= af
    cf .= c

    @test c == a
    @test af == cf
    @test typeof(a) == Array{Float64, 1}
    @test typeof(af) == AFArray{Float64, 1}
    @test typeof(c) == Array{Float64, 1}
    @test typeof(cf) == AFArray{Float64, 1}
end
