import Base: rand, randn, convert, diagm, eye, range, zeros, ones, trues, falses, deepcopy_internal
export constant, getSeed, setSeed, iota

function rand{T,N}(::Type{AFArray{T,N}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randu(ptr, dimensions, T)
    AFArray{T,N}(ptr[])
end
rand{T,N}(::Type{AFArray{T,N}}, t::Tuple{Vararg{Int64,N}}) = rand(AFArray{T,N}, t...)

function randn{T,N}(::Type{AFArray{T,N}}, dims::Integer...)
    ptr = new_ptr()
    dimensions = [dims...]
    af_randn(ptr, dimensions, T)
    AFArray{T,N}(ptr[])
end
randn{T,N}(::Type{AFArray{T,N}}, t::Tuple{Vararg{Int64,N}}) = randn(AFArray{T,N}, t...)

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

@compat (::Type{AFArray}){T,N}(a::Array{T,N}) = convert(AFArray{T,N}, a)
convert{T,N}(::Type{AFArray}, a::Array{T,N}) = AFArray(a)

function constant{T<:Real,N}(val::T, dims::Vararg{Int64,N})
    n = length(dims)
    dims = Int[dims...]
    for i = n+1:4
        push!(dims, 1)
    end
    ptr = new_ptr()
    af_constant!(ptr, val, n, dims, T)
    AFArray{T,N}(ptr[])
end

function constant{T<:Integer,N}(val::Complex{T}, dims::Vararg{Int64,N})
    n = length(dims)
    dims = [dims...]
    for i = n+1:4
        push!(dims, 1)
    end
    ptr = new_ptr()
    val2 = convert(Complex{Float32}, val)
    af_constant_complex!(ptr, val2, n, dims, Float32)
    AFArray{Complex{Float32},N}(ptr[])
end

function constant{T,N}(val::Complex{T}, dims::Vararg{Int64,N})
    n = length(dims)
    dims = [dims...]
    for i = n+1:4
        push!(dims, 1)
    end
    ptr = new_ptr()
    af_constant_complex!(ptr, val, n, dims, T)
    AFArray{Complex{T},N}(ptr[])
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
    af_iota(out, Cuint(ndims), dims, Cuint(t_ndims), tile_dims, T)
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
    constant(false, dims)
end
falses(::Type{AFArray{Bool}}, t::Tuple) = falses(AFArray{T}, t...)

function deepcopy_internal{T,N}(a::AFArray{T,N}, stackdict::ObjectIdDict)
    haskey(stackdict, a) && return stackdict[a]
    out = new_ptr()
    af_copy_array(out, a.ptr)
    c = AFArray{T,N}(out[])
    stackdict[a] = c
    return c
end
