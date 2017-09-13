
using ArrayFire
using Base.Test

@testset "Main" begin
    include("scope.jl")
    include("indexing.jl")
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
