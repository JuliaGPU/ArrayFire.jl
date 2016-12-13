using ArrayFire
using Base.Test

#Basic math
a = rand(Float32, 10, 10)
ad = AFArray(a)
@test sumabs2(Array(ad + 2) - (a + 2)) < 1e-6
@test sumabs2(Array(ad .+ 2) - (a .+ 2)) < 1e-6
@test sumabs2(Array(2 + ad) - (2 + a)) < 1e-6
@test sumabs2(Array(2 .+ ad) - (2 .+ a)) < 1e-6
@test sumabs2(Array(ad - 2) - (a - 2)) < 1e-6
@test sumabs2(Array(ad .- 2) - (a .- 2)) < 1e-6
@test sumabs2(Array(2 - ad) - (2 - a)) < 1e-6
@test sumabs2(Array(2 .- ad) - (2 .- a)) < 1e-6
@test sumabs2(Array(ad * 2) - (a * 2)) < 1e-6
@test sumabs2(Array(ad .* 2) - (a .* 2)) < 1e-6
@test sumabs2(Array(2 * ad) - (2 * a)) < 1e-6
@test sumabs2(Array(2 .* ad) - (2 .* a)) < 1e-6
@test sumabs2(Array(ad / 2) - (a / 2)) < 1e-6
@test sumabs2(Array(ad ./ 2) - (a ./ 2)) < 1e-6
@test sumabs2(Array(2 ./ ad) - (2 ./ a)) < 1e-6
@test sumabs2(Array(ad .^ 2) - (a .^ 2)) < 1e-6
@test sumabs2(Array(2 .^ ad) - (2 .^ a)) < 1e-6

#Trig functions
@test sumabs2(Array(sin(ad)) - sin(a)) < 1e-6
@test sumabs2(Array(cos(ad)) - cos(a)) < 1e-6
@test sumabs2(Array(tan(ad)) - tan(a)) < 1e-6
@test sumabs2(Array(sinh(ad)) - sinh(a)) < 1e-6
@test sumabs2(Array(cosh(ad)) - cosh(a)) < 1e-6
@test sumabs2(Array(tanh(ad)) - tanh(a)) < 1e-6
@test sumabs2(Array(atan2(ad, ad)) - atan2(a, a)) < 1e-6
@test sumabs2(Array(atan2(ad, 2)) - atan2(a, 2)) < 1e-6
@test sumabs2(Array(atan2(2, ad)) - atan2(2, a)) < 1e-6

#Measures
@test minimum(a) == minimum(ad)
@test maximum(a) == maximum(ad)
@test sumabs2(Array(maximum(ad,1)) - maximum(a,1)) < 1e-6
@test sumabs2(Array(maximum(ad,2)) - maximum(a,2)) < 1e-6
@test sumabs2(Array(minimum(ad,1)) - minimum(a,1)) < 1e-6
@test sumabs2(Array(minimum(ad,2)) - minimum(a,2)) < 1e-6
@test sumabs2(Array(max(ad,0.5f0)) - max(a,0.5f0)) < 1e-6
@test sumabs2(Array(min(ad,0.5f0)) - min(a,0.5f0)) < 1e-6
@test_approx_eq  mean(ad) mean(a)
@test sumabs2(Array(mean(ad,1)) - mean(a,1)) < 1e-6
@test sumabs2(Array(mean(ad,2)) - mean(a,2)) < 1e-6
@test median(ad) == median(a)
@test sumabs2(Array(median(ad,1)) - median(a,1)) < 1e-6
@test_approx_eq var(ad) var(a)

#Linalg 
@test sumabs2(Array(ad') - a') < 1e-6
ld, ud, pd = lu(ad)
l, u, p = lu(a)
@test sumabs(Array(ld) - l) < 1e-5
@test sumabs(Array(ud) - u) < 1e-5 
@test sumabs(Array(pd) - p) < 1e-5 
@test sumabs2(Array(chol(ad*ad')) - chol(a*a')) < 1e-5
@test sumabs2(Array(ctranspose(chol(ad*ad'))) - ctranspose(chol(a*a'))) < 1e-5
ud, sd, vtd = svd(ad)
u, s, v = svd(a)
@test sumabs(Array(ud) - u) < 1e-4
@test sumabs(Array(sd) - s) < 1e-4
@test sumabs(Array(vtd') - v) < 1e-4

#Complex numbers 
@test Array(complex(ad,ad)) == complex(a,a)
@test Array(real(complex(ad,ad))) == real(complex(a,a))
@test Array(imag(complex(a,a))) == imag(complex(a,a))

# FFT - Issue #81
@test sumabs2(fft(a) - Array(fft(ad))) < 1e-6
@test sumabs2(ifft(a) - Array(ifft(ad, norm_factor = 0.01))) < 1e-6 # Note the scaling factor. Not sure why. 

# Inference
@test_throws MethodError ad | ad

# Deepcopy
ac = deepcopy(ad)
@test sumabs(ac - ad) == 0
ac[1] = 0
@test ac[1] != ad[1]

# Indexing - issue #96
let
    b = Float32[1.,2.,3.]
    bd = AFArray(b)
    ind = AFArray([false, true, true])
    @test Array(bd[ind]) == Float32[2., 3.]
end
