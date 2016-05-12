import Base: rand, randn, convert
export constant

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

function convert{T,N}(::Type{AFArray{T,N}}, a::Array{T,N}) 
    n = ndims(a)
    d = get_all_dims(a) 
    ptr = new_ptr()
    err = ccall((:af_create_array, af_lib), 
                Cint, (Ptr{Void}, Ptr{T}, Cuint, Ptr{Cuint}, Cint),
                ptr, pointer(a), n, pointer(d), aftype(T))
    err == 0 || throwAFerror(err)
    AFArray{T,N}(ptr[])
end
            
call{T,N}(::Type{AFArray}, a::Array{T,N}) = convert(AFArray{T,N}, a)
convert{T,N}(::Type{AFArray}, a::Array{T,N}) = AFArray(a)

function constant{T<:Real}(val::T, dims::Integer...)
    n = length(dims)
    dims = [dims...]
    for i = n+1:4
        push!(dims, 1)
    end
    ptr = new_ptr()
    af_constant!(ptr, val, n, dims, T)
    AFArray{T}(ptr[])
end 

function constant{T}(val::Complex{T}, dims::Integer...)
    n = length(dims)
    dims = [dims...]
    for i = n+1:4
        push!(dims, 1)
    end
    ptr = new_ptr()
    if T <: Integer
        val = convert(Complex{Float32}, val)
    end
    af_constant_complex!(ptr, val, n, dims, Float32)
    if T <: Integer
        AFArray{Complex{Float32}}(ptr[])
    else
        AFArray{Complex{T}}(ptr[])
    end
end
