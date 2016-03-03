using ArrayFire
using Base.Test

#Basic math
a = rand(Float32, 10, 10)
ad = AFArray(a)
@test sumabs2(Array(ad + 2) - (a + 2)) < 1e-6
@test sumabs2(Array(ad - 2) - (a - 2)) < 1e-6
@test sumabs2(Array(ad * 2) - (a * 2)) < 1e-6
@test sumabs2(Array(ad / 2) - (a / 2)) < 1e-6
