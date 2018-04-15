using ArrayFire
using Test
using Libdl,Random,SparseArrays,LinearAlgebra
using FFTW

@testset "Main" begin
    include("scope.jl")

    @testset "Bugs" begin
        include("bugs.jl")
    end

    allowslow(AFArray, false)

    @testset "FFT" begin
        include("fft.jl")
    end

    allowslow(AFArray) do
        include("indexing.jl")
    end
    include("sparse.jl")
    include("math.jl")
    include("blackscholes.jl")
    include("array.jl")
end
