import Base: getindex, setindex!, lastindex

if VERSION < v"0.6-"
    Base.LinearIndexing(::Type{AFArray}) = Base.LinearSlow()
else
    Base.IndexStyle(::Type{AFArray}) = Base.IndexCartesian()
end

export allowslow

const _allowslow = Ref(true)
allowslow(::Type{AFArray}, flag = true) = (_allowslow[] = flag)
function allowslow(f, ::Type{AFArray}, flag = true)
    old, _allowslow[] = _allowslow[], flag
    try
        return f()
    finally
        _allowslow = old
    end
end

assertslow(op) = _allowslow[] || error("$op is disabled")

create_seq(r::UnitRange) = af_seq(r.start-1, r.stop-1, 1)
create_seq(r::StepRange) = af_seq(r.start-1, r.stop-1, r.step)
create_seq(i::Int) = af_seq(i-1, i-1, 1)
create_seq(::Colon) = af_seq(1, 1, 0)

set_indexer!(indexers, i, s::Union{AbstractRange,Int,Colon}) = set_seq_indexer(indexers, create_seq(s), i, true)
set_indexer!(indexers, i, s::AFArray{Bool}) = set_array_indexer(indexers, find(s)-UInt32(1), i)
set_indexer!(indexers, i, s::AFArray) = set_array_indexer(indexers, s-UInt32(1), i)
set_indexer!(indexers, i, s::Vector) = set_indexer!(indexers, i, AFArray(s))

function set_seq_indexer(indexer, idx, dim::dim_t, is_batch::Bool)
    _error(ccall((:af_set_seq_indexer,af_lib),af_err,
                 (Ptr{af_index_t},Ptr{af_seq},dim_t,Bool),
                 indexer, RefValue{af_seq}(idx), dim, is_batch))
end

function set_array_indexer(indexer, idx, dim::dim_t)
    _error(ccall((:af_set_array_indexer,af_lib),af_err,
                 (Ptr{af_index_t},af_array,dim_t),
                 indexer,idx.arr,dim))
end

function release_indexers(indexers)
    _error(ccall((:af_release_indexers,af_lib),af_err,(Ptr{af_index_t},), indexers))
end

function assign_gen!(lhs::AFArray,ndims::dim_t,indices,rhs::AFArray)
    out = RefValue{af_array}(lhs.arr)
    _error(ccall((:af_assign_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t},af_array),out,lhs.arr,ndims,indices,rhs.arr))
    if lhs.arr != out[]
        release_array(lhs)
        lhs.arr = out[]
    end
end

function create_indexers(idx)
    indexers = create_indexers()
    for (i, thing) in enumerate(idx)
        set_indexer!(indexers, i-UInt32(1), thing)
    end
    indexers
end

function getindex(a::AFArray{T}, idx::Union{AbstractRange,Colon,AFArray,Vector}, idx1::Int...) where {T}
    assertslow("getindex")
    indexers = create_indexers((idx, idx1...))
    out = index_gen_1(a, length(idx1)+1, indexers)
    release_indexers(indexers)
    out
end

function getindex(a::AFArray{T}, idx0::Union{AbstractRange,Colon,AFArray,Vector,Int},
                  idx::Union{AbstractRange,Colon,AFArray,Vector}, idx1::Int...) where {T}
    assertslow("getindex")
    indexers = create_indexers((idx0, idx, idx1...))
    out = index_gen_2(a, length(idx1)+2, indexers)
    release_indexers(indexers)
    out
end

function index_gen_1(_in::AFArray{T,N},ndims::dim_t,indices) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_index_gen,af_lib),
                 af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t}),
                 out,_in.arr,ndims,indices))
    AFArray{T,1}(out[])
end

function index_gen_2(_in::AFArray{T,N},ndims::dim_t,indices) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_index_gen,af_lib),
                 af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t}),
                 out,_in.arr,ndims,indices))
    AFArray{T,2}(out[])
end

function getindex(a::AFArray{T}, idx::Int...) where {T}
    assertslow("getindex")
    @assert length(idx) <= length(size(a))
    indexers = create_indexers(idx)
    out = index_gen(a, length(idx), indexers)
    release_indexers(indexers)
    Array(out)[1]
end

function setindex!(lhs::AFArray{T}, rhs::AFArray{S}, idx::Union{AbstractRange,Int,Colon,AFArray}...) where {T,S}
    assertslow("setindex!")
    @assert length(idx) <= length(size(lhs))
    indexers = create_indexers(idx)
    if T == S
        assign_gen!(lhs, length(idx), indexers, rhs)
    else
        assign_gen!(lhs, length(idx), indexers, convert(AFArray{T}, rhs))
    end
    release_indexers(indexers)
    rhs
end

function setindex!(lhs::AFArray{T}, val::S, idx::Union{AbstractRange,Int,Colon,AFArray}...) where {T,S}
    assertslow("setindex!")
    sz = get_sizes(idx, lhs)
    rhs = constant(T(val), sz)
    setindex!(lhs, rhs, idx...)
    val
end

function get_sizes(idx::Tuple, lhs::AFArray)
    s = Array{Int}(undef, length(idx))
    for i = 1:length(idx)
        if typeof(idx[i]) <: AbstractRange
            s[i] = length(idx[i])
        elseif idx[i] == Colon()
            s[i] = size(lhs,i)
        elseif typeof(idx[i]) <: Integer
            s[i] = 1
        elseif typeof(idx[i]) <: AFArray{Bool}
            s[i] = sum(idx[i])
        elseif typeof(idx[i]) <: AFArray
            s[i] = length(idx[i])
        end
    end
    (s...,)
end

lastindex(a::AFArray) = length(a)
lastindex(a::AFArray, d) = size(a, d)
