import Base: getindex, setindex!

immutable seq
    start::Cdouble
    stop::Cdouble
    step::Cdouble
end

function create_seq_object(r::Range)
    start = r.start - 1
    stop = r.stop - 1
    if :step in fieldnames(r)
        step = r.step
    else
        step = 1
    end
    s = seq(start, stop, step)
end

function create_seq_object(i::Int)
    start = i - 1
    stop = i - 1
    step = 1
    s = seq(start, stop, step)
end

function create_seq_object(::Colon)
    start = 1
    stop = 1
    step = 0
    s = seq(start, stop, step)
end

immutable index
    ptr::Ptr{Void}
end

function getindex{T}(a::AFArray{T}, idx::Union{Range,Int,Colon,AFArray}...)
    indexers = new_ptr()
    af_create_indexers(indexers)
    indexers = index(indexers[])
    set_up_index!(idx, indexers)
    out = new_ptr()
    n = ndims(a)
    l = length(idx)
    l <= n || throw(DimensionMismatch("Number of dimensions is lesser than number of indices"))
    af_index_gen(out, a, l, indexers)
    af_release_indexers(indexers)
    outarr = AFArray{T}(out[])
    if length(outarr) == 1
        return Array(outarr)[1]
    else
        return outarr
    end
end 

function set_up_index!(idx::Tuple, indexers::index)
    for (i,thing) in enumerate(idx)
        t = typeof(thing)
        if t <: Range || t <: Int || t <: Colon
            s = create_seq_object(thing)
            af_set_seq_indexer(indexers, s, i-1, true)
        elseif t <: AFArray{Bool}
            _shift = find(thing) - 1
            af_set_array_indexer(indexers, _shift, i-1)
        elseif t <: AFArray
            _shift = thing - 1
            af_set_array_indexer(indexers, _shift, i-1)
        end
    end
end

function setindex!{T,S}(lhs::AFArray{T}, rhs::AFArray{S}, idx::Union{Range,Int,Colon,AFArray}...)
    indexers = new_ptr()
    af_create_indexers(indexers)
    indexers = index(indexers[])
    set_up_index!(idx, indexers)
    n = ndims(lhs)
    l = length(idx)
    rhs = convert(AFArray{T}, rhs)
    l <= n || throw(DimensionMismatch("Number of dimensions is lesser than number of indices"))
    af_assign_gen(pointer_from_objref(lhs), lhs, l, indexers, rhs)
    af_release_indexers(indexers)
    rhs
end

function setindex!{T,S<:Real}(lhs::AFArray{T}, val::S, idx::Union{Range,Int,Colon,AFArray}...)
    sz = get_sizes(idx, lhs)
    rhs = constant(T(val), sz...)
    setindex!(lhs, rhs, idx...)
    val
end

function get_sizes(idx::Tuple, lhs::AFArray)
    s = Array(Int, length(idx))
    for i = 1:length(idx)
        if typeof(idx[i]) <: Range
            s[i] = length(idx[i])
        elseif idx[i] == Colon()
            s[i] = size(lhs,i)
        elseif typeof(idx[i]) <: Integer
            s[i] = 1
        elseif typeof(idx[i]) <: AFArray
            s[i] = length(idx[i])
        end
    end
    tuple(s...)
end

Base.LinearIndexing(::Type{AFArray}) = Base.LinearSlow()
