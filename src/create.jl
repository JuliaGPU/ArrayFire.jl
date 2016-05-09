import Base: rand, randn

function rand{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randu!(ptr, dimensions, T)
    AFArray{T}(ptr[])
end

function randn{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randn!(ptr, dimensions, T)
    AFArray{T}(ptr[])
end

