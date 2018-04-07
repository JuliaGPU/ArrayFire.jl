
using ArrayFire
using Test
using Libdl,Random,SparseArrays,LinearAlgebra
using FFTW

allowslow(AFArray, false)

@testset "Main" begin
    include("scope.jl")
    allowslow(AFArray) do
        include("indexing.jl")
    end
    include("sparse.jl")
    include("math.jl")
    include("blackscholes.jl")
    include("autodiff.jl")
    include("array.jl")
end

@testset "Bugs" begin
    include("bugs.jl")
end

@testset "FFT" begin
    include("fft.jl")
end
