export @afscope, afscope

macro afscope(expr)
    esc(expr)
end

global const scopes = Vector{Vector{AFArray}}()

function enter_scope()
    push!(scopes, Vector{AFArray}())
end

function leave_scope(num)
    scope = pop!(scopes)
    for k = 1:length(scope) - num
        finalize(scope[k])
    end
    if !(isempty(scopes))
        append!(scopes[end], scope[end-num+1:end])
    end
end

function afscope(f, num)
    enter_scope()
    try
        f()
    finally
        leave_scope(num)
    end
end
