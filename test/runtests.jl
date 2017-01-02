using ArrayFire
using Base.Test

@test af_info() == nothing
@test_throws ErrorException af_set_device(Cint(5))
