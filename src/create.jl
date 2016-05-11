import Base: rand, randn, convert

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
