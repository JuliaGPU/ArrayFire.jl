module ArrayFire

using Compat

export AFArray

include("config.jl")


type AFArray{T,N} <: AbstractArray{T,N}
    ptr::Ptr{Void}
    function AFArray(ptr::Ptr{Void})
        a = new(ptr)
        finalizer(a, af_release_array)
        a
    end
end

typealias AFVector{T} AFArray{T,1}
typealias AFMatrix{T} AFArray{T,2}

include("error.jl")
include("defns.jl")
include("indexing.jl")
include("wrap.jl")
include("utils.jl")
include("create.jl")
include("math.jl")
include("vector.jl")
include("backend.jl")
include("stats.jl")
include("linalg.jl")
include("image.jl")
include("cv.jl")
include("signal.jl")
include("device.jl")

end
