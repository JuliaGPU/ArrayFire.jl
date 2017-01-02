using ArrayFire
using Base.Test

@test af_info() == nothing
@test unsafe_string(af_err_to_string(Cuint(0))) == "Success"
@test_throws ErrorException af_set_device(Cint(-5))
@test af_get_manual_eval_flag() == false
@test af_set_manual_eval_flag(true) == nothing
@test af_get_manual_eval_flag() == true
