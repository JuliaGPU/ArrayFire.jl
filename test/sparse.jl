@testset "Sparse" begin
A = AFArray(rand(10,10));
@test !issparse(A)

# build an identity matrix
n = 10;

Ie = AFArray(eye(10))
Aid = sparse(Ie)

@test issparse(Aid)

Ie2 = full(Aid)

@test !issparse(Ie2)

@test all(Ie2 == Ie)

e = sprand(10, 20, 0.2)
c = randn(20)
a1 = AFArray(e)
c1 = AFArray(c)
@test issparse(a1)

e2 = SparseMatrixCSC(a1)
@test e2 == e

@test sum(e * c) ≈ sum(a1 * c1)
end
