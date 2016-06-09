using ArrayFire
function warmup()
    a = rand(10,10)
    ad = AFArray(a)
    ad * ad
    a * a
    fft(a)
    fft(ad)
    chol(a*a')
    chol(ad*ad') 
    rand(AFArray{Float64}, 10, 10)
    sort(vec(a))
    sort(vec(ad))
end

function matmul(a::AFArray, b::AFArray)
    for i = 1:10
        r = a * b
        ArrayFire.eval(r)
    end
    sync()
end

function fast_fourier(a::AFArray)
    for i = 1:10
        r = fft(a)
        ArrayFire.eval(r)
    end
    sync()
end

function cholesky(a::AFArray)
    for i = 1:10
        r = chol(a)
        ArrayFire.eval(r)
    end
    sync()
end

function random()
    for i = 1:10
        r = rand(AFArray{Float32}, 5000, 5000)
        ArrayFire.eval(r)
    end
    sync()
end

function sorting(a::AFArray)
    for i = 1:10
        r = sort(a)
        ArrayFire.eval(r)
    end
    sync()
end

function benchmark()
    
    warmup()
    info("Warmup done!")
    gc()

    a = rand(Float32, 2000, 2000) 
    ad = AFArray(a)

    #Matmul
    info("Matmul")
    t1 = @elapsed a * a
    t2 = @elapsed matmul(ad, ad)
    println("Time (CPU): $t1") 
    println("Time (GPU): $(t2/10)") 
    gc()

    #FFT
    info("FFT")
    t1 = @elapsed fft(a)
    t2 = @elapsed fast_fourier(ad)
    println("Time (CPU): $t1") 
    println("Time (GPU): $(t2/10)") 

    #Cholesky
    b = a * a' 
    bd = AFArray(b)
    info("Cholesky")
    t1 = @elapsed chol(b)
    t2 = @elapsed cholesky(bd)
    println("Time (CPU): $t1") 
    println("Time (GPU): $(t2/10)") 

    #Rand
    info("Rand")
    t1 = @elapsed rand(Float32, 5000, 5000)
    t2 = @elapsed random()
    println("Time (CPU): $t1") 
    println("Time (GPU): $(t2/10)") 
    
    #Vecsort
    c = rand(Float32, 10^6)
    cd = AFArray(c)
    info("Vec sort")
    t1 = @elapsed sort(c)
    t2 = @elapsed sorting(cd)
    println("Time (CPU): $t1") 
    println("Time (GPU): $(t2/10)") 

end

benchmark()
