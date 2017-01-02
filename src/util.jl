function af_err_to_string(err::af_err)
    ccall((:af_err_to_string, af_lib), Cstring, (Cint, ), err)
end

function af_error(err::af_err)
    if err != 0
        str = af_err_to_string(err)
        throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))"))
    end
end
