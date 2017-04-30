
using ArrayFire
using Base.Test

@testset "Main" begin
    include("indexing.jl")
    include("sparse.jl")
    include("math.jl")
end

gc()
device_gc()
