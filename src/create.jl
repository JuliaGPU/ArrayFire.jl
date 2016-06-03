import Base: rand, randn, convert, diagm, eye, range, zeros, ones, trues, falses
export constant, getSeed, setSeed, iota

function rand{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randu(ptr, dimensions, T)
    AFArray{T}(ptr[])
end
rand{T}(::Type{AFArray{T}}, t::Tuple) = rand(AFArray{T}, t...)

function randn{T}(::Type{AFArray{T}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randn(ptr, dimensions, T)
    AFArray{T}(ptr[])
end
randn{T}(::Type{AFArray{T}}, t::Tuple) = randn(AFArray{T}, t...)

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
constant{T}(val::T, t::Tuple) = constant(val, t...)

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
eye{T}(::Type{AFArray{T}}, t::Tuple) = eye(AFArray{T}, t...)

function iota{T}(::Type{AFArray{T}}, dims::Integer...; tile_dims = [1])
    out = new_ptr()
    t_ndims = length(tile_dims)
    dims = [dims...]
    ndims = length(dims)
    af_iota(out, Cuint(ndims), dims, Cuint(t_ndims), tdims, T)
    AFArray{T}(out[])
end
iota{T}(::Type{AFArray{T}}, t::Tuple) = iota(AFArray{T}, t...)
iota{T}(::Type{AFArray{T}}, t::Tuple, tile_dims) = iota(AFArray{T}, t..., tile_dims)

function range{T}(::Type{AFArray{T}}, dims::Integer...; seq_dim = -1)
    out = new_ptr()
    dims = [dims...]
    l = length(dims)
    af_range(out, Cuint(l), dims, seq_dim, T)
    AFArray{T}(out[])
end
range{T}(::Type{AFArray{T}}, t::Tuple; seq_dim = -1) = range(AFArray{T}, t..., seq_dim)

function zeros{T}(::Type{AFArray{T}}, dims::Integer...)
    constant(T(0), dims)
end
zeros{T}(::Type{AFArray{T}}, t::Tuple) = zeros(AFArray{T}, t...)

function ones{T}(::Type{AFArray{T}}, dims::Integer...)
    constant(T(1), dims)
end
ones{T}(::Type{AFArray{T}}, t::Tuple) = ones(AFArray{T}, t...)

function trues(::Type{AFArray{Bool}}, dims::Integer...)
    constant(true, dims)
end
trues(::Type{AFArray{Bool}}, t::Tuple) = trues(AFArray{T}, t...)

function falses(::Type{AFArray{Bool}}, dims::Integer...)
    constant(true, dims)
end
falses(::Type{AFArray{Bool}}, t::Tuple) = falses(AFArray{T}, t...)
