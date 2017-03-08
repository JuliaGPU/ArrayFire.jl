A = AFArray(rand(10,10));
@test !issparse(A)

# build an identity matrix
n = 10;

Ie = AFArray(eye(10))
Aid = sparse(Ie)

@test issparse(Aid)

e2 = SparseMatrixCSC(Aid)

@test e2 == sparse(eye(10))

Ie2 = full(Aid)

@test !issparse(Ie2)

@test all(Ie2 == Ie)

I = AFArray(vec(collect(1:n)))
J = AFArray(vec(collect(1:n)))
V = AFArray(ones(n))
Aid2 = create_sparse_array(n, n, V, I, J, AF_STORAGE_CSR)

@test issparse(Aid2)
