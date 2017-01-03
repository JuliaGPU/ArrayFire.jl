using Clang.cindex
using Clang.wrap_c
using Compat

AF_INCLUDE = "/usr/include/af"

af_header = ["/usr/include/arrayfire.h"]

# Set up include paths
clang_includes = [".", AF_INCLUDE]

# Callback to test if a header should actually be wrapped (for exclusion)
function wrap_header(top_hdr, cursor_header)
    return startswith(dirname(cursor_header), AF_INCLUDE)
end

lib_file(hdr) = "af_lib"
output_file(hdr) = "../src/wrap.jl"

function wrap_cursor(name, cursor)
    if startswith(name, "AF") || startswith(name, "DEPRECATED") ||
        name == "bool" || isempty(name) || startswith(name, "SIZE_T") || endswith(name, ".h")
        return false
    end
    true
end

function return_val(typ, arg)
    if typ.args[2] == :af_array
        Expr(:call, :AFArray!, Expr(:ref, arg))
    else
        Expr(:ref, arg)
    end
end

function return_val(typ, arg, expr)
    if typ.args[2] == :af_array
        Expr(:call, expr, Expr(:ref, arg))
    else
        Expr(:ref, arg)
    end
end

function rewrite(line::Expr)
    if line.head == :function
        hdr = line.args[1].args
        name = hdr[1]
        args = hdr[2:end]
        body = line.args[2].args
        types = body[1].args[3].args
        vals = view(body[1].args, 4:length(body[1].args))

        num_input_arrays = 0
        num_output_arrays = 0

        for k in 1:length(args)
            if isa(args[k], Expr)
                if args[k].args[2] == :af_array
                    args[k].args[2] = :AFArray
                    vals[k] = Expr(:., vals[k], QuoteNode(:arr))
                    num_input_arrays += 1
                elseif args[k].args[2] == :Cint
                    args[k].args[2] = :Integer
                    vals[k] = Expr(:call, :Cint, vals[k])
                elseif args[k].args[2] == :UInt32
                    args[k].args[2] = :Integer
                    vals[k] = Expr(:call, :UInt32, vals[k])
                elseif args[k].args[2] == :Cdouble
                    args[k].args[2] = :Real
                    vals[k] = Expr(:call, :Cdouble, vals[k])
                end
            end
        end

        ret_type = body[1].args[2]
        if ret_type == :af_err
            body[1] = Expr(:call, :af_error, body[1])
        end
        num_out = 0
        for k = 1:length(types)
            t = types[k]
            if isa(t, Expr) && t.args[1] == :Ptr && t.args[2] == :af_array
                num_output_arrays += 1
            end
            if isa(t, Expr) && t.args[1] == :Ptr && t.args[2] != :Void && t.args[2] != :Cstring
                num_out = k
            else
                break
            end
        end
        if num_out > 0
            if num_input_arrays == 1 && num_output_arrays == 1 && num_out == 1 &&
                !contains("$name", "_sparse") && !contains("$name", "_true") && !contains("$name", "_is")
                hdr[1] = Expr(:curly, hdr[1], :T, :N)
                for k = 1:length(args)
                    if isa(args[k], Expr) && args[k].args[2] == :AFArray
                        args[k].args[2] = Expr(:curly, :AFArray, :T, :N)
                    end
                end
                push!(body, return_val(types[1], args[1], Expr(:curly, :AFArray, :T, :N)))
            elseif num_out == 1
                if num_output_arrays > 0
                    @printf("No type inference for %s\n", name)
                end
                push!(body, return_val(types[1], args[1]))
            else
                if num_output_arrays > 0
                    @printf("No type inference for %s\n", name)
                end
                push!(body, Expr(:tuple, map(k->return_val(types[k], args[k]), 1:num_out)...))
            end
            for k in num_out:-1:1
                deleteat!(hdr, 2)
                t = Expr(:curly, :RefValue, types[k].args[2])
                c = Expr(:call, t, 0)
                unshift!(body, Expr(:(=), args[k], c))
            end
        end
        return [line, Expr(:export, name)]
    end
    line
end

rewrite(line) = line

function rewriter(input)
    vcat(map(rewrite, input)...)
end

const wc = wrap_c.init(;
                       headers             = af_header,
                       common_file         = "../src/common.jl",
                       clang_includes      = clang_includes,
                       header_wrapped      = wrap_header,
                       header_library      = lib_file,
                       header_outputfile   = output_file,
                       rewriter            = rewriter,
                       cursor_wrapped      = wrap_cursor)

run(wc)
