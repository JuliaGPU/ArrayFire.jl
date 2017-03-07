A = AFArray(rand(10,10));
@test !issparse(A)

# build an identity matrix
n = 10;

Ie = AFArray(eye(10))
Aid = sparse(Ie)

@test issparse(Aid)

I = AFArray(vec(collect(1:n)))
J = AFArray(vec(collect(1:n)))
V = AFArray(ones(n))
Aid2 = create_sparse_array(n, n, V, I, J, AF_STORAGE_CSR)

@test issparse(Aid2)
