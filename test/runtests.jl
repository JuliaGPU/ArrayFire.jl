
using ArrayFire
using Base.Test

@testset "FFT" begin
    include("fft.jl")
end

@testset "Main" begin
    include("scope.jl")
    include("indexing.jl")
    include("sparse.jl")
    include("math.jl")
    include("blackscholes.jl")
    include("autodiff.jl")
end

gc()
device_gc()
