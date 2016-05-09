import Base: rand, randn

function rand{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    err = ccall((:af_randu, af_lib), 
                Cint, (Ptr{Ptr{Void}}, Cuint, Ptr{dim}, Cuint),
                ptr, length(dimensions), pointer(dimensions), aftype(T))
    err == 0 || throwAFerror(err)
    AFArray{T}(ptr[])
end
                
function randn{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    err = ccall((:af_randn, af_lib), 
                Cint, (Ptr{Ptr{Void}}, Cuint, Ptr{dim}, Cuint),
                ptr, length(dimensions), pointer(dimensions), aftype(T))
    err == 0 || throw(ErrorException("$err"))
    AFArray{T}(ptr[])
end
    
