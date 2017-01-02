using ArrayFire
using Base.Test

@test af_info() == nothing
@test isa(af_err_to_string(Cuint(0)), Cstring)
@test_throws ErrorException af_set_device(Cint(-5))
