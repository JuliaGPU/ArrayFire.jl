import Base: elsize, size, ndims, convert, showarray

immutable Dim4
    dim1::Integer
    dim2::Integer
    dim3::Integer
    dim4::Integer
end

new_ptr() = Base.RefValue{Ptr{Void}}(C_NULL)

sizeof{T}(a::AFArray{T}) = elsize(a) * length(a)
    
function convert{T,N}(::Type{Array{T,N}}, x::AFAbstractArray{T,N})
    ret = Array(UInt8, sizeof(x))
    ccall((:af_get_data_ptr, af_lib), 
            Cint, (Ptr{T}, Ptr{Ptr{Void}}), 
            pointer(ret), x.ptr)
    ret = reinterpret(T, ret)
    ret = reshape(ret, size(x)...)
    ret
end

function size(a::AFArray)
    dim1 = Base.RefValue{Cuint}(0)
    dim2 = Base.RefValue{Cuint}(0)
    dim3 = Base.RefValue{Cuint}(0)
    dim4 = Base.RefValue{Cuint}(0)   
    ccall((:af_get_dims, af_lib),
            Cint, (Ptr{Cuint}, Ptr{Cuint}, Ptr{Cuint}, Ptr{Cuint}, Ptr{Void}),
            dim1, dim2, dim3, dim4, a.ptr)
    n = ndims(a)
    dim4_to_dims(Dim4(dim1[], dim2[], dim3[], dim4[]), n)
end

function dim4_to_dims(d::Dim4, n::Integer)
    if n == 1
        return (d.dim1,)
    elseif n == 2
        return (Int(d.dim1), Int(d.dim2))
    elseif n == 3
        return (Int(d.dim1), Int(d.dim2), Int(d.dim3))
    elseif n == 4
        return (Int(d.dim1), Int(d.dim2), Int(d.dim3), Int(d.dim4))
    else
        throw(ArgumentError("Too many dimensions for an ArrayFire Array"))
    end
end
        
function ndims(a::AFArray)
    n = Base.RefValue{Cuint}(0)
    err = ccall((:af_get_numdims, af_lib),
            Cint, (Ptr{Cuint}, Ptr{Void}),
            n, a.ptr)
    Int(n[])
end

function ndims(ptr::Ptr{Void})
    n = Base.RefValue{Cuint}(0)
    err = ccall((:af_get_numdims, af_lib),
            Cint, (Ptr{Cuint}, Ptr{Void}),
            n, ptr)
    Int(n[])
end 

call{T}(::Type{AFArray{T}}, ptr::Ptr{Void}) = AFArray{T, ndims(ptr)}(ptr)

function showarray{T,N}(io::IO, X::AFAbstractArray{T,N};
                   header::Bool=true, kwargs...)
    header && print(io, summary(X))
    !isempty(X) && println(io,":")
    showarray(io, convert(Array{T,N},X); header = false, kwargs...)
end
