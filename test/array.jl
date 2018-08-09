@testset "Array" begin
    @testset "range" begin
        x = range(AFArray{Float32}, 1, 5)
        @test Array(x) == Array{Float32}(collect(1.0:5.0))
        @test typeof(x) == AFArray{Float32, 1}
        x = range(AFArray{Float32}, 3, 5)
        @test Array(x) == Array{Float32}(collect(3.0:7.0))
        @test typeof(x) == AFArray{Float32, 1}
        x = range(AFArray{Float32}, 2, 3, 10)
        @test Array(x) == collect(2.0:3:29.0)
        @test typeof(x) == AFArray{Float32, 1}

        @test sum(@inferred iota((2,3))) == 21
    end
end
