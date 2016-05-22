import Base: rand, randn, convert, diagm, eye, range
export constant, getSeed, setSeed, iota

function rand{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randu(ptr, dimensions, T)
    AFArray{T}(ptr[])
end

function randn{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randn(ptr, dimensions, T)
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

function diagm{T}(val::AFVector{T}, k::Integer = 1)
    out = new_ptr()
    af_diag_create(out, val, k - 1)
    AFArray{T}(out[])
end

function getSeed()
    out = Base.Ref{Cuint}(0)
    af_get_seed(out)
    Int(out[])
end

function setSeed(s::Int)
    af_set_seed(Cuint(s))
    s
end

function eye{T}(::Type{AFArray{T}}, dims::Integer...)
    out = new_ptr()
    ndims = length(dims)
    dims = [dims...]
    af_identity(out, Cuint(ndims), dims, T)
    AFArray{T}(out[])
end

function iota{T}(::Type{AFArray{T}}, dims::Integer...; tile_dims = [1])
    out = new_ptr()
    t_ndims = length(tile_dims)
    dims = [dims...]
    ndims = length(dims)
    af_iota(out, Cuint(ndims), dims, Cuint(t_ndims), tdims, T)
    AFArray{T}(out[])
end

function range{T}(::Type{AFArray{T}}, dims::Integer...; seq_dim = -1)
    out = new_ptr()
    dims = [dims...]
    l = length(dims)
    af_range(out, Cuint(l), dims, seq_dim, T)
    AFArray{T}(out[])
end
