module ArrayFire

export AFArray

include("config.jl")

abstract AFAbstractArray{T,N} <: AbstractArray{T,4}

immutable AFArray{T,N} <: AFAbstractArray{T,N}
    ptr::Ptr{Void}
end


include("defns.jl")
include("utils.jl")
include("error.jl")
include("create.jl")

end
