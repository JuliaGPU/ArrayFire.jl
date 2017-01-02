function af_error(err::af_err)
    if err != 0
        str = af_err_to_string(err)
        throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))"))
    end
    nothing
end
