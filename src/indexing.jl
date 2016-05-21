import Base: getindex

immutable seq
    start::Cdouble
    stop::Cdouble
    step::Cdouble
end

function create_seq_object(r::Range, ind::Int)
    start = r.start - 1
    stop = r.stop - 1
    if :step in fieldnames(r)
        step = r.step
    else
        step = 1
    end
    s = seq(start, stop, step)
end

function create_seq_object(i::Int, ind::Int)
    start = i - 1
    stop = i - 1
    step = 1
    s = seq(start, stop, step)
end

function create_seq_object(::Colon, i::Int)
    start = 0
    stop = i - 1
    step = 1
    s = seq(start, stop, step)
end

function getindex{T}(a::AFArray{T}, idx::Union{Range,Int,Colon}...)
    out = new_ptr()
    l = length(idx)
    n = ndims(a)
    l <= n || throw(DimensionMismatch("Number of dimensions is lesser than number of indices"))
    s = [create_seq_object(idx[i], size(a,i)) for i = 1:l]
    af_index(out, a, Cuint(l), s)
    AFArray{T}(out[])
end

