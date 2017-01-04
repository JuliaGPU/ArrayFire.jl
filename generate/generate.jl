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

const renames = Dict("sign" => "signbit", "product" => "prod", "init" => "afinit", "info" => "afinfo",
                     "copy_array" => "copy", "get_version" => "afversion", "eval" => "afeval",
                     "min" => "minimum", "max" => "maximum", "any_true" => "any", "all_true" => "all")

const ignore = Set(["example_function", "create_array", "retain_array", "get_data_ref_count", "info_string",
                    "device_info", "alloc_host", "free_host", "alloc_pinned", "free_pinned",
                    "get_type", "get_numdims", "join_many", "eval_multiple", "get_size_of", "cast",
                    "constant", "constant_complex", "constant_long", "constant_ulong"])

const booleans = Set(["lt", "gt", "le", "ge", "eq", "neq", "iszero", "isinf", "isnan"])

const maths = Set(["add", "sub", "mul", "div", "rem", "mod", "atan2", "root", "pow", "dot"])

const recast = Dict(:Cint => :Integer, :UInt32 => :Integer, :Cdouble => :Real)

function rewrite(line::Expr)
    if line.head == :function
        hdr = line.args[1].args
        name = replace("$(hdr[1])", "af_", "", 1)
        name = get(renames, name, name)
        in(name, ignore) && return []
        hdr[1] = Symbol(name)
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
                elseif args[k].args[2] == :af_dtype
                    args[k].args[2] = :Type
                    vals[k] = Expr(:call, :af_type, vals[k])
                elseif haskey(recast, args[k].args[2])
                    key = args[k].args[2]
                    val = recast[key]
                    args[k].args[2] = val
                    vals[k] = Expr(:call, key, vals[k])
                end
            end
        end

        ret_type = body[1].args[2]
        if ret_type == :af_err
            body[1] = Expr(:call, :_error, body[1])
        end
        num_out = 0
        for k = 1:length(types)
            t = types[k]
            if isa(t, Expr) && t.args[1] == :Ptr && t.args[2] == :af_array
                num_output_arrays += 1
            end
            if isa(t, Expr) && t.args[1] == :Ptr && t.args[2] != :Void
                num_out = k
            else
                break
            end
        end
        if num_out > 0
            if num_output_arrays == 1 && num_out == 1 && name in booleans
                hdr[1] = Expr(:curly, hdr[1], :T, :N)
                for k = 1:length(args)
                    if isa(args[k], Expr) && args[k].args[2] == :AFArray
                        args[k].args[2] = Expr(:curly, :AFArray, :T, :N)
                    end
                end
                push!(body, return_val(types[1], args[1], Expr(:curly, :AFArray, :Bool, :N)))
            elseif num_input_arrays == 1 && num_output_arrays == 1 && num_out == 1 &&
                !contains(name, "_sparse") && !contains(name, "_true")
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
                c = types[k].args[2] == :Cstring ? Expr(:call, t) : Expr(:call, t, 0)
                unshift!(body, Expr(:(=), args[k], c))
            end
        end
        return [line, Expr(:export, Symbol(name))]
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
