using ArrayFire
using Base.Test

@test afinfo() == nothing
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
