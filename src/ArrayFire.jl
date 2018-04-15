__precompile__()
module ArrayFire

using Libdl,Random,SparseArrays,LinearAlgebra,FFTW,DSP

include("common.jl")
include("array.jl")
include("util.jl")
include("wrap.jl")
include("random.jl")
include("graphics.jl")
include("scope.jl")
include("indexing.jl")
include("fft.jl")

end
