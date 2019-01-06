using ArrayFire, LinearAlgebra, Logging

function warmup()
    a = rand(Float32, 10,10)
    ad = AFArray(a)
    ad * ad
    a * a
    fft(a)
    fft(ad)
    cholesky(a*a')
    cholesky(ad*ad')
    rand(AFArray{Float32}, 10, 10)
    sort(vec(a))
    sort(vec(ad))
end

function matmul(a::AFArray, b::AFArray)
    for i = 1:10
        r = a * b
        sync(r)
    end
end

function fast_fourier(a::AFArray)
    for i = 1:10
        r = fft(a)
        sync(r)
    end
end

function chol(a::AFArray)
    for i = 1:10
        r,s = cholesky(a)
        sync(r)
    end
end

function random()
    for i = 1:10
        r = rand(AFArray{Float32}, 5000, 5000)
        sync(r)
    end
end

function sorting(a::AFArray)
    for i = 1:10
        r = sort(a)
        sync(r)
    end
end

function benchmark()

    warmup()
    @info("Warmup done!")
    GC.gc()

    a = rand(Float32, 2000, 2000)
    ad = AFArray(a)

    #Matmul
    @info("Matmul")
    t1 = @elapsed a * a
    t2 = @elapsed matmul(ad, ad)
    println("Time (CPU): $t1")
    println("Time (GPU): $(t2/10)")
    GC.gc()

    #FFT
    @info("FFT")
    t1 = @elapsed fft(a)
    t2 = @elapsed fast_fourier(ad)
    println("Time (CPU): $t1")
    println("Time (GPU): $(t2/10)")

    #Cholesky
    b = a * a' + 2000 * Matrix{Float32}(I, 2000, 2000)
    bd = AFArray(b)
    @info("Cholesky")
    t1 = @elapsed cholesky(b)
    t2 = @elapsed chol(bd)
    println("Time (CPU): $t1")
    println("Time (GPU): $(t2/10)")

    #Rand
    @info("Rand")
    t1 = @elapsed rand(Float32, 5000, 5000)
    t2 = @elapsed random()
    println("Time (CPU): $t1")
    println("Time (GPU): $(t2/10)")

    #Vecsort
    c = rand(Float32, 10^6)
    cd = AFArray(c)
    @info("Vec sort")
    t1 = @elapsed sort(c)
    t2 = @elapsed sorting(cd)
    println("Time (CPU): $t1")
    println("Time (GPU): $(t2/10)")

end

benchmark()
