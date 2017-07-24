#Basic math
a = rand(Float32, 10, 10)
ad = AFArray(a)

@testset "Math" begin
    @test sum(abs2, Array(ad + 2) - (a + 2)) < 1e-6
    @test sum(abs2, Array(ad .+ 2) - (a .+ 2)) < 1e-6
    @test sum(abs2, Array(2 + ad) - (2 + a)) < 1e-6
    @test sum(abs2, Array(2 .+ ad) - (2 .+ a)) < 1e-6
    @test sum(abs2, Array(ad - 2) - (a - 2)) < 1e-6
    @test sum(abs2, Array(ad .- 2) - (a .- 2)) < 1e-6
    @test sum(abs2, Array(2 - ad) - (2 - a)) < 1e-6
    @test sum(abs2, Array(2 .- ad) - (2 .- a)) < 1e-6
    @test sum(abs2, Array(ad * 2) - (a * 2)) < 1e-6
    @test sum(abs2, Array(ad .* 2) - (a .* 2)) < 1e-6
    @test sum(abs2, Array(2 * ad) - (2 * a)) < 1e-6
    @test sum(abs2, Array(2 .* ad) - (2 .* a)) < 1e-6
    @test sum(abs2, Array(ad / 2) - (a / 2)) < 1e-6
    @test sum(abs2, Array(ad ./ 2) - (a ./ 2)) < 1e-6
    @test sum(abs2, Array(2 ./ ad) - (2 ./ a)) < 1e-6
    @test sum(abs2, Array(ad .^ 2) - (a .^ 2)) < 1e-6
    @test sum(abs2, Array(2 .^ ad) - (2 .^ a)) < 1e-6

    #Trig functions
    @test sum(abs2, Array(sin.(ad)) - sin.(a)) < 1e-6
    @test sum(abs2, Array(cos.(ad)) - cos.(a)) < 1e-6
    @test sum(abs2, Array(tan.(ad)) - tan.(a)) < 1e-6
    @test sum(abs2, Array(sinh.(ad)) - sinh.(a)) < 1e-6
    @test sum(abs2, Array(cosh.(ad)) - cosh.(a)) < 1e-6
    @test sum(abs2, Array(tanh.(ad)) - tanh.(a)) < 1e-6

    @test sum(fill(AFArray, 1, (1, 2))) == 2
    @test sum(fill(AFArray{Float32}, 1, (1, 2))) == 2
    @test sum(fill(AFArray{Float32,2}, 1, (1, 2))) == 2
    @test sum(ones(AFArray{Float32}, (1, 2))) == 2f0
    @test sum(ones(AFArray{Float32,2}, (1, 2))) == 2f0
    @test eltype(zeros(AFArray{Float32}, (1, 2))) == Float32
    @test eltype(zeros(AFArray{Float32,2}, (1, 2))) == Float32
    @test typeof(zeros(ad)) == AFArray{Float32, 2}
    @test typeof(ones(ad)) == AFArray{Float32, 2}
end

@testset "Arrays" begin
    @test err_to_string(Cuint(0)) == "Success"
    @test_throws ErrorException set_device(-5)
    arr1 = @inferred AFArray{Int,1}([1, 2])
    @test eltype(arr1) == Int
    @test ndims(arr1) == 1
    @test (@inferred size(arr1)) == (2,)
    arr2 = @inferred copy(arr1)
    arr3 = @inferred deepcopy(arr2)
    @test typeof(arr3) == typeof(arr2)
    arr4 = @inferred Array{Int,1}(arr1)
    @test arr4 == [1, 2]
    arr5 = AFArray{Int,2}([1 2 3; 4 5 6])
    @test (@inferred size(arr5)) == (2,3)
    @test size(@inferred cat(2, arr1, arr5)) == (2, 4)
    @test size(@inferred hcat(arr1, arr5)) == (2, 4)
    @test size(@inferred vcat(arr1', arr5')) == (4, 2)
    arr6 = @inferred AFArray([1., 2.])
    arr7 = @inferred Array(arr6)
    @test @inferred eltype(arr6) == eltype(arr7)
    @test @inferred ndims(arr6) == ndims(arr7)
    @test @inferred size(arr6) == size(arr7)
    @test @inferred any(@inferred isnan(arr6)) == false
    @test @inferred any(@inferred isinf(arr6)) == false
    @test @inferred any(@inferred any(isnan(arr6), 1)) == false
    @test @inferred any(@inferred any(isinf(arr6), 1)) == false
    @test @inferred all(@inferred isnan(arr6)) == false
    @test @inferred all(@inferred isinf(arr6)) == false
    @test @inferred all(@inferred all(isnan(arr6), 1)) == false
    @test @inferred all(@inferred all(isinf(arr6), 1)) == false
    @test sum(arr5) == 21
    @test sum(arr6) == 3.
    @test sum(complex(arr6)) == 3 + 0im
    @test sum(complex(arr6, arr6)) == 3 + 3im
    @test sum(abs2(arr5)) == 91
    @test sum(abs2(arr6)) == 5.0
    arr9 = AFArray{Complex{Float64},2}([1.0+3im 2. 3.; 4 5 6])
    @test sum(arr9) == 21+3im
    @test sum(abs2(arr9)) == 100.0
    @test sum(conj(arr9)) == 21-3im
    @test sum(real(arr9)) == 21
    @test sum(imag(arr9)) == 3
    @test typeof(@inferred AFArray{Complex{Float32},2}(arr9)) == AFArray{Complex{Float32},2}
    @test typeof(@inferred AFArray{UInt32}(arr9)) == AFArray{UInt32,2}
    b = @inferred(1 + arr1)
    @test sum(b) == 5
end

a1 = rand(2, 3) + rand(2, 3)im
a2 = rand(2, 3) + rand(2, 3)im
af1 = @inferred AFArray(a1)
af2 = @inferred AFArray(a2)

for op in [:+, :-, :.+, :.-, :.*, :./, :.^]
    @test @eval sum($op(a1, a2)) ≈ sum(@inferred $op(af1, af2))
end

a3 = randn(Float32, 2, 3)
af3 = AFArray(a3)
a4 = randn(Float32, 2, 3)
af4 = AFArray(a4)

for op in [:.>, :.>=, :.!=, :.<, :.<=]
    @test @eval all($op(a3, a4) .== Array(@inferred $op(af3, af4)))
    @test @eval all($op(0, a4) .== Array(@inferred $op(0, af4)))
    @test @eval all($op(a3, 0) .== Array(@inferred $op(af3, 0)))
end

for op in [:.==]
    @test @eval all($op(a3, a4) .== Array($op(af3, af4)))
    @test @eval all($op(0, a4) .== Array($op(0, af4)))
    @test @eval all($op(a3, 0) .== Array($op(af3, 0)))
end

c1 = rand()

for op in [:+, :-, :*, :.+, :.-, :.*, :./, :.^]
    @test @eval sum($op(c1, a2)) ≈ sum(@inferred $op(c1, af2))
end

for op in [:+, :-, :*, :/, :.^, :.+, :.-, :.*, :./]
    @test @eval sum($op(a1, c1)) ≈ sum(@inferred $op(af1, c1))
end

as = Int32[-2 -1 0 1 2]
@test all(AFArray(sign.(as)) .== sign(AFArray(as)))
@test !any(!(AFArray(sign.(as)) .== sign(AFArray(as))))

@test (@inferred size(AFArray(rand(1,2,3)))) == (1,2,3)
@test (@inferred size(AFArray(rand(1,2,3,4)))) == (1,2,3,4)

am = rand(Float32, 3, 3)
amf = AFArray(am)

@test sum(amf * amf) ≈ sum(am * am)
@test sum(amf * amf') ≈ sum(am * am')
@test sum(A_mul_Bt(amf, amf)) ≈ sum(am * am')
@test sum(amf' * amf) ≈ sum(am' * am)
@test sum(At_mul_B(amf, amf)) ≈ sum(am' * am)
@test sum(amf' * amf') ≈ sum(am' * am')
@test sum(At_mul_Bt(amf, amf)) ≈ sum(am' * am')

@test Array(sum(amf, 1)) ≈ sum(am, 1)
@test Array(sum(amf, 2)) ≈ sum(am, 2)

@test Array(prod(amf, 1)) ≈ prod(am, 1)
@test Array(prod(amf, 2)) ≈ prod(am, 2)

@test all(Array(minimum(amf, 1)) .== minimum(am, 1))
@test all(Array(minimum(amf, 2)) .== minimum(am, 2))

@test all(Array(maximum(amf, 1)) .== maximum(am, 1))
@test all(Array(maximum(amf, 2)) .== maximum(am, 2))
@test maximum(amf) == maximum(am)
@test minimum(amf) == minimum(am)
@test mean(amf) == mean(am)
@test std(amf) ≈ std(am)
@test var(amf) ≈ var(am)
@test median(amf) == median(am)

@test Array(vec(amf)) == vec(am)
@test typeof(norm(amf)) == Float32
@test @inferred(norm(amf)) ≈ norm(am)
u,s,v = @inferred(svd(amf))

signal = rand(Float32, 100)
signalf = AFArray(signal)
filt = rand(Float32, 3)
filtf = AFArray(filt)
@test sum(conv(signal, filt)) ≈ sum(@inferred conv(signalf, filtf))

srand(AFArray, rand(Int))

@testset "Sizes" begin
    @test size(rand(AFArray, 1)) == (1,)
    @test size(rand(AFArray, 1, 2)) == (1,2)
    @test size(rand(AFArray{Float64}, 1)) == (1,)
    @test size(rand(AFArray{Float64}, 1, 2)) == (1,2)
    @test size(rand(AFArray{Float64}, 1)) == (1,)
    @test size(rand(AFArray{Float64}, 1, 2)) == (1,2)

    @test size(randn(AFArray, 1)) == (1,)
    @test size(randn(AFArray, 1, 2)) == (1,2)
    @test size(randn(AFArray{Float64}, 1)) == (1,)
    @test size(randn(AFArray{Float64}, 1, 2)) == (1,2)
    @test size(randn(AFArray{Float64}, 1)) == (1,)
    @test size(randn(AFArray{Float64}, 1, 2)) == (1,2)
end

a3 = rand(UInt8, 2, 3)
af3 = AFArray(a3)
a4 = rand(UInt8, 2, 3)
af4 = AFArray(a4)
c1 = rand(UInt8)

for op in [:&, :|, :xor, :min, :max]
    @test @eval all($op.(a3, a4) .== Array(@inferred $op(af3, af4)))
    @test @eval all($op.(c1, a4) .== Array(@inferred $op(c1, af4)))
    @test @eval all($op.(a3, c1) .== Array(@inferred $op(af3, c1)))
end

a3 = rand(UInt64, 2, 3)
af3 = AFArray(a3)
a4 = rand(UInt64, 2, 3)
af4 = AFArray(a4)
c1 = rand(UInt64)

for op in [:&, :|, :xor, :min, :max]
    @test @eval all($op.(a3, a4) .== Array(@inferred $op(af3, af4)))
    @test @eval all($op.(c1, a4) .== Array(@inferred $op(c1, af4)))
    @test @eval all($op.(a3, c1) .== Array(@inferred $op(af3, c1)))
end

a3 = abs.(rand(Int64, 2, 3))
af3 = AFArray(a3)
a4 = mod.(rand(Int64, 2, 3), 10)
af4 = AFArray(a4)
c1 = mod(rand(Int64), 10)

for op in [:<<, :>>, :&, :|, :xor, :min, :max]
    @test @eval all($op.(a3, a4) .== Array($op.(af3, af4)))
    @test @eval all($op.(a3, a4) .== Array(@inferred $op(af3, af4)))
    @test @eval all($op.(c1, a4) .== Array(@inferred $op(c1, af4)))
    @test @eval all($op.(a3, c1) .== Array(@inferred $op(af3, c1)))
end

a3 = rand(Bool, 2, 3)
af3 = AFArray(a3)
a4 = rand(Bool, 2, 3)
af4 = AFArray(a4)
c1 = rand(Bool)

for op in [:&, :|, :xor, :min, :max]
    @test @eval all($op.(a3, a4) .== Array(@inferred $op(af3, af4)))
    @test @eval all($op.(c1, a4) .== Array(@inferred $op(c1, af4)))
    @test @eval all($op.(a3, c1) .== Array(@inferred $op(af3, c1)))
end

a5 = rand(Float32, 5)
af5 = AFArray(a5)

for op in [:abs, :acos, :asin, :atan, :cbrt, :ceil, :cos, :cosh,
           :exp, :expm1, :factorial, :floor,
           :imag, :isinf, :isnan, :lgamma, :log, :log10,
           :log1p, :log2, :real, :round, :sign, :signbit, :sin, :sinh, :sqrt, :tan,
           :tanh, :trunc]
    @testset "$op" begin
        @test @eval $op.(a5) ≈ Array(@inferred $op(af5))
    end
end

a6 = rand(Float32, 3, 3)
af6 = AFArray(a6)

for op in [:identity, :cov, :vec, :zeros, :ones, :transpose]
    @testset "$op" begin
        @test @eval $op(a6) ≈ Array(@inferred $op(af6))
    end
end

for op in [:maximum, :minimum, :sum, :prod, :size, :sizeof]
    @testset "$op" begin
        @test @eval $op(a6) ≈ @inferred $op(af6)
    end
end

a = AFArray([true false true])
b = AFArray([1 2 3])
c = AFArray([4 5 6; 7 8 9])

@test Array(select(a, b, c)) == [1 5 3; 1 8 3]
@test Array(select(a, 0, b)) == [0 2 0]
@test Array(select(a, b, 0)) == [1 0 3]

@test size(reshape(c, 6)) == (6, )
@test size(reshape(c, (3,2))) == (3, 2)

#plot(rand(AFArray, 10), rand(AFArray, 10))

s = rand(AFArray{Float32}, 10)
sh = Array(s)
@testset "Sort" begin
    @test Array(@inferred sort(s)) == sort(sh)
    @test Array(@inferred sortperm(s)) == sortperm(sh)
end
