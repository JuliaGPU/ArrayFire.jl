using ArrayFire
using Base.Test

@test afinfo() == nothing
@test unsafe_string(err_to_string(Cuint(0))) == "Success"
@test_throws ErrorException set_device(-5)
@test get_manual_eval_flag() == false
@test set_manual_eval_flag(true) == nothing
@test get_manual_eval_flag() == true
arr1 = @inferred AFArray{Int,1}([1, 2])
arr2 = @inferred copy(arr1)
arr3 = @inferred deepcopy(arr2)
@test typeof(arr3) == typeof(arr2)
