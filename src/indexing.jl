import Base: getindex, setindex!

Base.IndexStyle(::Type{AFArray}) = Base.IndexCartesian()

create_seq(r::UnitRange) = af_seq(r.start-1, r.stop-1, 1)
create_seq(r::StepRange) = af_seq(r.start-1, r.stop-1, r.step)
create_seq(i::Int) = af_seq(i-1, i-1, 1)
create_seq(::Colon) = af_seq(1, 1, 0)

set_indexer!(indexers, i, s::Union{Range,Int,Colon}) = set_seq_indexer(indexers, create_seq(s), i, true)
set_indexer!(indexers, i, s::AFArray) = set_array_indexer(indexers, s-1, i)

function set_seq_indexer(indexer, idx, dim::dim_t, is_batch::Bool)
    _error(ccall((:af_set_seq_indexer,af_lib),af_err,(Ptr{af_index_t},Ptr{af_seq},dim_t,Bool),
                 indexer, RefValue{af_seq}(idx), dim, is_batch))
end

function release_indexers(indexers)
    _error(ccall((:af_release_indexers,af_lib),af_err,(Ptr{af_index_t},), indexers))
end

function assign_gen(lhs::AFArray,ndims::dim_t,indices,rhs::AFArray)
    out = RefValue{af_array}(lhs.arr)
    _error(ccall((:af_assign_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t},af_array),out,lhs.arr,ndims,indices,rhs.arr))
end

function create_indexers(idx)
    indexers = create_indexers()
    for (i, thing) in enumerate(idx)
        set_indexer!(indexers, i-1, thing)
    end
    indexers
end

function getindex{T}(a::AFArray{T}, idx::Union{Range,Int,Colon,AFArray}...)
    @assert length(idx) <= length(size(a))
    indexers = create_indexers(idx)
    out = index_gen(a, length(idx), indexers)
    release_indexers(indexers)
    out
end

function setindex!{T,S}(lhs::AFArray{T}, rhs::AFArray{S}, idx::Union{Range,Int,Colon,AFArray}...)
    @assert length(idx) <= length(size(lhs))
    indexers = create_indexers(idx)
    if T == S
        assign_gen(lhs, length(idx), indexers, rhs)
    else
        assign_gen(lhs, length(idx), indexers, convert(AFArray{T}, rhs))
    end
    release_indexers(indexers)
    rhs
end
