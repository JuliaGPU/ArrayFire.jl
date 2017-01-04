using ArrayFire
using Base.Test

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
@test eltype(arr6) == eltype(arr7)
@test ndims(arr6) == ndims(arr7)
@test size(arr6) == size(arr7)
@test any(isnan(arr6)) == false
@test any(isinf(arr6)) == false
@test any(any(isnan(arr6), 1)) == false
@test any(any(isinf(arr6), 1)) == false
@test all(isnan(arr6)) == false
@test all(isinf(arr6)) == false
@test all(all(isnan(arr6), 1)) == false
@test all(all(isinf(arr6), 1)) == false
@test sum(arr5) == 21
@test sum(arr6) == 3.
arr9 = AFArray{Complex{Float64},2}([1.0+3im 2. 3.; 4 5 6])
@test sum(arr9) == 21+3im
