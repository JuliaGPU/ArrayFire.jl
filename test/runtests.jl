
using ArrayFire
using Base.Test

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
end

@testset "Bugs" begin
    include("bugs.jl")
end

@testset "FFT" begin
    include("fft.jl")
end

gc()
device_gc()
