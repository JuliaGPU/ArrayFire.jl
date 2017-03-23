a1 = AFArray(rand(3, 2))
b1 = AFArray([1. 2.])
a1[1, :] = b1
@test all(a1[1, :] == b1)
