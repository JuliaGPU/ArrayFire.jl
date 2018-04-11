@testset "Array" begin
    @testset "range" begin
        x = range(AFArray{Float32}, 1, 5)
        @test Array(x) == Array{Float32}(collect(range(1, 5)))
        @test typeof(x) == AFArray{Float32, 1}
        x = range(AFArray{Float32}, 3, 5)
        @test Array(x) == Array{Float32}(collect(range(3, 5)))
        @test typeof(x) == AFArray{Float32, 1}
        x = range(AFArray{Float32}, 2, 3, 10)
        @test Array(x) == collect(range(2,3,10))
        @test typeof(x) == AFArray{Float32, 1}

        @test sum(@inferred iota((2,3))) == 21
    end
end
