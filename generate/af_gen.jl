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
output_file(hdr) = "../src/af_wrap.jl"

function wrap_cursor(name, cursor)
    !startswith(name, "AF") && !startswith(name, "DEPRECATED") &&
    name != "bool" && !isempty(name) && !startswith(name, "SIZE_T")
end

const wc = wrap_c.init(;
                       headers             = af_header,
                       common_file         = "../src/af_common.jl",
                       clang_includes      = clang_includes,
                       header_wrapped      = wrap_header,
                       header_library      = lib_file,
                       header_outputfile   = output_file,
                       cursor_wrapped      = wrap_cursor)

run(wc)
