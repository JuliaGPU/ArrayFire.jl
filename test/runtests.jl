using ArrayFire
using Base.Test

include("autodiff.jl")

#Basic math
a = rand(Float32, 10, 10)
ad = AFArray(a)
@test sumabs2(Array(@inferred ad + 2) - (a + 2)) < 1e-6
@test sumabs2(Array(@inferred ad .+ 2) - (a .+ 2)) < 1e-6
@test sumabs2(Array(@inferred 2 + ad) - (2 + a)) < 1e-6
@test sumabs2(Array(@inferred 2 .+ ad) - (2 .+ a)) < 1e-6
@test sumabs2(Array(@inferred ad - 2) - (a - 2)) < 1e-6
@test sumabs2(Array(@inferred ad .- 2) - (a .- 2)) < 1e-6
@test sumabs2(Array(@inferred 2 - ad) - (2 - a)) < 1e-6
@test sumabs2(Array(@inferred 2 .- ad) - (2 .- a)) < 1e-6
@test sumabs2(Array(@inferred ad * 2) - (a * 2)) < 1e-6
@test sumabs2(Array(@inferred ad .* 2) - (a .* 2)) < 1e-6
@test sumabs2(Array(@inferred 2 * ad) - (2 * a)) < 1e-6
@test sumabs2(Array(@inferred 2 .* ad) - (2 .* a)) < 1e-6
@test sumabs2(Array(@inferred ad / 2) - (a / 2)) < 1e-6
@test sumabs2(Array(@inferred ad ./ 2) - (a ./ 2)) < 1e-6
@test sumabs2(Array(@inferred 2 ./ ad) - (2 ./ a)) < 1e-6
@test sumabs2(Array(@inferred ad .^ 2) - (a .^ 2)) < 1e-6
@test sumabs2(Array(@inferred 2 .^ ad) - (2 .^ a)) < 1e-6

#Trig functions
@test sumabs2(Array(@inferred sin(ad)) - sin(a)) < 1e-6
@test sumabs2(Array(@inferred cos(ad)) - cos(a)) < 1e-6
@test sumabs2(Array(@inferred tan(ad)) - tan(a)) < 1e-6
@test sumabs2(Array(@inferred sinh(ad)) - sinh(a)) < 1e-6
@test sumabs2(Array(@inferred cosh(ad)) - cosh(a)) < 1e-6
@test sumabs2(Array(@inferred tanh(ad)) - tanh(a)) < 1e-6
@test sumabs2(Array(@inferred atan2(ad, ad)) - atan2(a, a)) < 1e-6
@test sumabs2(Array(@inferred atan2(ad, 2)) - atan2(a, 2)) < 1e-6
@test sumabs2(Array(@inferred atan2(2, ad)) - atan2(2, a)) < 1e-6

#Measures
@test minimum(a) == @inferred minimum(ad)
@test maximum(a) == @inferred maximum(ad)
@test sumabs2(Array(@inferred maximum(ad,1)) - maximum(a,1)) < 1e-6
@test sumabs2(Array(@inferred maximum(ad,2)) - maximum(a,2)) < 1e-6
@test sumabs2(Array(@inferred minimum(ad,1)) - minimum(a,1)) < 1e-6
@test sumabs2(Array(@inferred minimum(ad,2)) - minimum(a,2)) < 1e-6
@test sumabs2(Array(@inferred max(ad,0.5f0)) - max(a,0.5f0)) < 1e-6
@test sumabs2(Array(@inferred min(ad,0.5f0)) - min(a,0.5f0)) < 1e-6
@test_approx_eq  @inferred(mean(ad)) mean(a)
@test sumabs2(Array(@inferred mean(ad,1)) - mean(a,1)) < 1e-6
@test sumabs2(Array(@inferred mean(ad,2)) - mean(a,2)) < 1e-6
@test @inferred median(ad) == median(a)
@test sumabs2(Array(@inferred median(ad,1)) - median(a,1)) < 1e-6
@test_approx_eq @inferred(var(ad)) var(a)

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

ac = deepcopy(ad)
@test sumabs(ac - ad) == 0
ac[1] = 0
@test ac[1] != ad[1]

c1 = Complex(-1,2)
@test_approx_eq @inferred(constant(c1, 1))  c1
c2 = Complex(-1.,2.)
@test_approx_eq @inferred(constant(c2, 1))  c2
