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
    if startswith(name, "AF") || startswith(name, "DEPRECATED") || name == "af_err_to_string" ||
        name == "bool" || isempty(name) || startswith(name, "SIZE_T") || endswith(name, ".h")
        return false
    end
    true
end

function rewrite(line::Expr)
    if line.head == :function
        name = line.args[1].args[1]
        body = line.args[2].args
        body[1] = Expr(:call, :af_error, body[1])
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
