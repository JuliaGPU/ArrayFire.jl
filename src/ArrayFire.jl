module ArrayFire

export AFArray

include("config.jl")

abstract AFAbstractArray{T,N} <: AbstractArray{T,4}

immutable AFArray{T,N} <: AFAbstractArray{T,N}
    ptr::Ptr{Void}
end


include("error.jl")
include("defns.jl")
include("wrap.jl")
include("utils.jl")
include("create.jl")
include("math.jl")
include("vector.jl")

end
