import Base: elsize, size, ndims, convert, showarray, vec, flipdim, vcat, hcat, cat, reshape, permutedims, circshift, repeat

export AFInfo, replace!, mergeArrays

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
    err = ccall((:af_get_data_ptr, af_lib), 
                Cint, (Ptr{T}, Ptr{Ptr{Void}}), 
                pointer(ret), x.ptr)
    err == 0 || throwAFerror(err)
    #af_get_data_ptr!(ret, x, T)
    ret = reinterpret(T, ret)
    ret = reshape(ret, size(x)...)
    ret
end

call{T,N}(::Type{Array}, a::AFArray{T,N}) = convert(Array{T,N}, a)

function size(a::AFArray)
    dim1 = Base.RefValue{Cuint}(0)
    dim2 = Base.RefValue{Cuint}(0)
    dim3 = Base.RefValue{Cuint}(0)
    dim4 =Base.RefValue{Cuint}(0)   
    af_get_dims!(dim1, dim2, dim3, dim4, a)
    n = ndims(a)
    dim4_to_dims(Dim4(dim1[], dim2[], dim3[], dim4[]), n)
end

function dim4_to_dims(d::Dim4, n::Integer)
    if n == 1
        return (Int(d.dim1),)
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

function get_all_dims(a::Array)
    n = ndims(a)
    dims = [size(a)...]
    for i = n+1:4
        push!(dims, 1)
    end
    dims
end
        
function ndims(a::AFArray)
    n = Base.RefValue{Cuint}(0)
    af_get_numdims!(n, a.ptr) 
    Int(n[])
end

function ndims(ptr::Ptr{Void})
    n = Base.RefValue{Cuint}(0)
    af_get_numdims!(n, ptr)
    Int(n[])
end 

call{T}(::Type{AFArray{T}}, ptr::Ptr{Void}) = AFArray{T, ndims(ptr)}(ptr)

function showarray{T,N}(io::IO, X::AFAbstractArray{T,N};
                   header::Bool=true, kwargs...)
    header && print(io, summary(X))
    !isempty(X) && println(io,":")
    showarray(io, convert(Array{T,N},X); header = false, kwargs...)
end

function backend_eltype(p::Ptr{Void})
    t = Base.Ref{Cuint}(0)
    af_get_type(t, p)
    jltype(t[])
end

backend_eltype(a::AFArray) = backend_eltype(a.ptr)

constant(val, t::Tuple) = constant(val, t...)

function AFInfo() 
    af_info()
    nothing
end

# Array Utils

function vec{T}(a::AFArray{T})
    out = new_ptr()
    af_vec(out, a)
    AFArray{T}(out[])
end

function flipdim{T}(a::AFArray{T}, dim::Int)
    out = new_ptr()
    af_flip(out, a, Cuint(dim-1))
    AFArray{T}(out[])
end

function vcat{T,S}(a::AFArray{T}, b::AFArray{S})
    out = new_ptr()
    af_join(out, 0, a, b)
    AFArray{af_promote(T,S)}(out[])
end

function hcat{T,S}(a::AFArray{T}, b::AFArray{S})
    out = new_ptr()
    af_join(out, 1, a, b)
    AFArray{af_promote(T,S)}(out[])
end

function cat(dim::Int, a::AFArray...)
    out = new_ptr()
    n = length(a)
    b = [a[i].ptr for i = 1:n]
    af_join_many(out, dim - 1, Cuint(n), b)
    AFArray{backend_eltype(out[])}(out[])
end

function reshape{T}(a::AFArray{T}, dims::Int...)
    out = new_ptr()
    l = length(dims)
    dims = [dims...]
    af_moddims(out, a, Cuint(l), dims)
    AFArray{T}(out[])
end

function permutedims{T}(a::AFArray{T}, perm::Vector{Int})
    make4Dperm!(perm)
    out = new_ptr()
    perm = perm - 1
    af_reorder(out, a, Cuint(perm[1]), Cuint(perm[2]), Cuint(perm[3]), Cuint(perm[4]))
    AFArray{T}(out[])
end

function make4Dperm!(perm::Vector{Int})
    l = length(perm)
    for i = l+1:4
        push!(perm, i)
    end
end

function replace!(a::AFArray, b::AFArray, cond::AFArray{Bool})
    af_replace(a, cond, b)
    nothing
end

function replace!(a::AFArray, b::Real, cond::AFArray{Bool})
    af_replace_scalar(a, cond, b)
    nothing
end

function mergeArrays{T,S}(a::AFArray{T}, b::AFArray{S}, cond::AFArray{Bool})
    out = new_ptr()
    af_select(out, cond, a, b)
    AFArray{af_promote(T,S)}(out[])
end

function mergeArrays{T<:Number,S}(a::T, b::AFArray{S}, cond::AFArray{Bool})
    out = new_ptr()
    af_select_scalar_l(out, cond, a, b)
    AFArray{af_promote(T,S)}(out[])
end

function mergeArrays{T,S<:Number}(a::AFArray{T}, b::S, cond::AFArray{Bool})
    out = new_ptr()
    af_select_scalar_r(out, cond, a, b)
    AFArray{af_promote(T,S)}(out[])
end

function circshift{T}(a::AFArray{T}, shifts::Vector{Int})
    make4Dshift!(shifts)
    out = new_ptr()
    af_shift(out, a, shifts...)
    AFArray{T}(out[])
end 

function make4Dshift!(s::Vector{Int})
    l = length(s)
    for i = l+1:4
        push!(s, 0)
    end
end

function repeat{T}(a::AFArray{T}; outer::Vector{Int} = [1,1,1,1])
    make4Dtile!(outer)
    out = new_ptr()
    outer = map(Cuint, outer)
    af_tile(out, a, outer...)
    AFArray{T}(out[])
end

function make4Dtile!(d::Vector{Int})
    l = length(d)
    for i = l+1:4
        push!(d, 1)
    end
end
