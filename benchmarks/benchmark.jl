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

function benchmark()
    
    warmup()
    info("Warmup done!")

    a = rand(1000, 1000) 
    ad = AFArray(a)

    #Matmul
    info("Matmul")
    t1 = @elapsed a * a
    t2 = @elapsed ad * ad
    println("Time (CPU): $t1") 
    println("Time (GPU): $t2") 

    #FFT
    info("FFT")
    t1 = @elapsed fft(a)
    t2 = @elapsed fft(ad)
    println("Time (CPU): $t1") 
    println("Time (GPU): $t2") 

    #Cholesky
    b = a * a'
    bd = AFArray(b)
    info("Cholesky")
    t1 = @elapsed chol(b)
    t2 = @elapsed chol(bd)
    println("Time (CPU): $t1") 
    println("Time (GPU): $t2") 

    #Rand
    info("Rand")
    t1 = @elapsed rand(1000, 1000)
    t2 = @elapsed rand(AFArray{Float64}, 1000, 1000)
    println("Time (CPU): $t1") 
    println("Time (GPU): $t2") 
    
    #Vecsort
    c = rand(10^6)
    cd = AFArray(c)
    info("Vec sort")
    t1 = @elapsed sort(c)
    t2 = @elapsed sort(cd)
    println("Time (CPU): $t1") 
    println("Time (GPU): $t2") 

end

benchmark()
