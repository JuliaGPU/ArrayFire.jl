function af_error(err::af_err)
    if err != 0
        str = ccall((:af_err_to_string, af_lib), Cstring, (Cint, ), err)
        throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))"))
    end
end
