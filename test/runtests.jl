using ArrayFire
using Base.Test

include("blackscholes.jl")

@test unsafe_string(err_to_string(Cuint(0))) == "Success"
@test_throws ErrorException set_device(-5)
@test get_manual_eval_flag() == false
@test set_manual_eval_flag(true) == nothing
@test get_manual_eval_flag() == true
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
arr9 = AFArray{Complex{Float64},2}([1.0+3im 2. 3.; 4 5 6])
@test sum(arr9) == 21+3im
@test typeof(@inferred AFArray{Complex{Float32},2}(arr9)) == AFArray{Complex{Float32},2}
@test typeof(@inferred AFArray{UInt32}(arr9)) == AFArray{UInt32,2}
b = @inferred(1 + arr1)
@test sum(b) == 5

a1 = rand(2, 3) + rand(2, 3)im
a2 = rand(2, 3) + rand(2, 3)im
af1 = @inferred AFArray(a1)
af2 = @inferred AFArray(a2)

for op in [:+, :-, :.+, :.-, :.*, :./, :.^]
    @assert @eval sum($op(a1, a2)) ≈ sum(@inferred $op(af1, af2))
end

a3 = randn(Float32, 2, 3)
af3 = AFArray(a3)
a4 = randn(Float32, 2, 3)
af4 = AFArray(a4)

for op in [:.>, :.>=, :.==, :.!=, :.<, :.<=]
    @assert @eval all($op(a3, a4) .== Array(@inferred $op(af3, af4)))
    @assert @eval all($op(0, a4) .== Array(@inferred $op(0, af4)))
    @assert @eval all($op(a3, 0) .== Array(@inferred $op(af3, 0)))
end

c1 = rand()

for op in [:+, :-, :*, :.+, :.-, :.*, :./]
    @assert @eval sum($op(c1, a2)) ≈ sum(@inferred $op(c1, af2))
end

for op in [:+, :-, :*, :/, :.^, :.+, :.-, :.*, :./]
    @assert @eval sum($op(a1, c1)) ≈ sum(@inferred $op(af1, c1))
end

as = Int32[-2 -1 0 1 2]
@test all(AFArray(sign.(as)) == sign(AFArray(as)))

@test typeof(@inferred fft_r2c(AFArray(rand(Float32, 10)), 1.0, 0)) == AFArray{Complex{Float32}, 1}
