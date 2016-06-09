using ArrayFire
function nmf{T<:Real}(A::AbstractMatrix{T}, k::Integer, n::Int=10)
    W = rand(T, size(A,1), k)
    H = rand(T, k, size(A,2))
    for j = 1:n
        H .*= (W'*A) ./ ((W'*W)*H)
        W .*= (A*H') ./ (W*(H*H'))
    end
    return W, H
end

function nmf_gpu{T<:Real}(A::AFArray{T}, k::Integer, n::Int=10)
    W = rand(AFArray{T}, size(A,1), k)
    H = rand(AFArray{T}, k, size(A,2))
    for j = 1:n
        H = H .* (W'*A) ./ ((W'*W)*H)
        W = W .* (A*H') ./ (W*(H*H'))
    end
    sync()
    return W, H
end

info("Warmup")
a = rand(Float32,10,10)
ad = AFArray(a)
nmf(a,2)
nmf_gpu(ad,2)
for i = 1:10
    dim = i*1000
    a = rand(Float32,dim,dim)
    ad = AFArray(a)
    info("Start Timing Dim = $dim !")
    t1 = @elapsed nmf(a,2)
    t2 = @elapsed nmf_gpu(ad,2)
    println("tcpu = $t1")
    println("tgpu = $t2")
    println("Speedup = $(t1/t2)")
    gc()
end
#=a = rand(1000,1000)
ad = AFArray(a)
W = rand(AFArray{Float64}, size(a,1), 2)
H = rand(AFArray{Float64}, 2, size(a,2))
W1,H1 = nmf(a, Array(W), Array(H), 2)
W2,H2 = nmf_gpu(ad, W, H, 2)
@show sum(W1-Array(W2))
@show sum(H1-Array(H2))=#
