using ArrayFire
using Base.Test

#Basic math
a = rand(Float32, 10, 10)
ad = AFArray(a)
@test sumabs2(Array(ad + 2) - (a + 2)) < 1e-6
@test sumabs2(Array(ad - 2) - (a - 2)) < 1e-6
@test sumabs2(Array(ad * 2) - (a * 2)) < 1e-6
@test sumabs2(Array(ad / 2) - (a / 2)) < 1e-6

#Trig functions
@test sumabs2(Array(sin(ad)) - sin(a)) < 1e-6
@test sumabs2(Array(cos(ad)) - cos(a)) < 1e-6
@test sumabs2(Array(tan(ad)) - tan(a)) < 1e-6
@test sumabs2(Array(sinh(ad)) - sinh(a)) < 1e-6
@test sumabs2(Array(cosh(ad)) - cosh(a)) < 1e-6
@test sumabs2(Array(tanh(ad)) - tanh(a)) < 1e-6

#Measures
@test minimum(a) == minimum(ad)
@test maximum(a) == maximum(ad)
@test sumabs2(Array(maximum(ad,1)) - maximum(a,1)) < 1e-6
@test sumabs2(Array(maximum(ad,2)) - maximum(a,2)) < 1e-6
@test sumabs2(Array(minimum(ad,1)) - minimum(a,1)) < 1e-6
@test sumabs2(Array(minimum(ad,2)) - minimum(a,2)) < 1e-6
@test sumabs2(Array(max(ad,0.5)) - max(a,0.5)) < 1e-6
@test sumabs2(Array(min(ad,0.5)) - min(a,0.5)) < 1e-6
@test sumabs2(Array(mean(ad,1)) - mean(a,1)) < 1e-6
@test sumabs2(Array(mean(ad,2)) - mean(a,2)) < 1e-6
@test sumabs2(Array(median(ad,1)) - median(a,1)) < 1e-6
@test sumabs2(Array(median(ad,2)) - median(a,2)) < 1e-6
@test sumabs2(Array(std(ad,1)) - std(a,1)) < 1e-6
@test sumabs2(Array(std(ad,2)) - std(a,2)) < 1e-6
@test sumabs2(Array(var(ad,1)) - var(a,1)) < 1e-6
@test sumabs2(Array(var(ad,2)) - var(a,2)) < 1e-6
