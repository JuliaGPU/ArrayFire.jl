macro afgc(expr)
    @assert expr.head == :function || expr.head == :(=) "Only works on functions or assignments"
    sc = (:(ArrayFire.scope() do ; end))
    sc.args[2].args[2] = expr.args[2]
    expr.args[2] = sc
    esc(expr)
end

global const scopes = Vector{Vector{AFArray}}()

matches(arr, except::Union{Void,Number,Array{Number}}) = false
matches(arr, except::AFArray) = arr === except
matches(arr, except::Tuple) = any(ex -> matches(arr, ex), except)
matches(arr, except::Any) = error("@afgc return value can be Void, Number, Array, or Tuple")

function push_to_scope(arr)
    if !isempty(scopes)
        push!(scopes[end], arr)
    end
    return arr
end

function leave_scope(except)
    scope = pop!(scopes)
    for arr in scope
        if matches(arr, except)
            push_to_scope(arr)
        else
            finalize(arr)
        end
    end
    return except
end

function scope(f)
    push!(scopes, Vector{AFArray}())
    local except
    try
        except = f()
    catch
        leave_scope(nothing)
        rethrow()
    end
    return leave_scope(except)
end

export @afgc
