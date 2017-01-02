using ArrayFire
using Base.Test

@test af_info() == nothing
@test unsafe_string(af_err_to_string(Cuint(0))) == "Success"
@test_throws ErrorException af_set_device(Cint(-5))
