#This example has been adopted from https://github.com/IntelLabs/ParallelAccelerator.jl/blob/master/examples/black-scholes/black-scholes.jl

using ArrayFire, Printf, Random

function blackscholes_serial(sptprice, strike, rate, volatility, time)
    logterm = log10( sptprice / strike)
    powterm = .5f0 * volatility * volatility
    den = volatility * sqrt(time)
    d1 = (((rate + powterm) * time) + logterm) / den
    d2 = d1 - den
    NofXd1 = cndf2(d1)
    NofXd2 = cndf2(d2)
    futureValue = strike * exp(- rate * time)
    c1 = futureValue * NofXd2
    call_ = sptprice * NofXd1 - c1
    put  = call_ - futureValue + sptprice
end

@inline function cndf2(in)
    out = 0.5f0 + 0.5f0 * erf(0.707106781f0 * in)
    return out
end

function runs(iterations)
    sptprice   = Float32[ 42.0 for i = 1:iterations ]
    initStrike = Float32[ 40.0 + (i / iterations) for i = 1:iterations ]
    rate       = Float32[ 0.5 for i = 1:iterations ]
    volatility = Float32[ 0.2 for i = 1:iterations ]
    time       = Float32[ 0.5 for i = 1:iterations ]

    sptprice_gpu = AFArray(sptprice)
    initStrike_gpu = AFArray(initStrike)
    rate_gpu = AFArray(rate)
    volatility_gpu = AFArray(volatility)
    time_gpu = AFArray(time)

    t1 = @elapsed put1 = blackscholes_serial.(sptprice, initStrike, rate, volatility, time)
    #    println("Serial checksum: ", sum(put1))
    t2 = @elapsed put2 = sync(blackscholes_serial.(sptprice_gpu, initStrike_gpu, rate_gpu, volatility_gpu, time_gpu))
    #    println("Parallel checksum: ", sum(put2))
    @assert sum(put1) ≈ sum(put2)
    return t1, t2
end

function driver(pwr)
    iterations = 10^pwr
    Random.seed!(0)
    blackscholes_serial.(Float32[], Float32[], Float32[], Float32[], Float32[])
    blackscholes_serial.(AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]),
                         AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]))
    tserial, tp1 = runs(iterations)
    tserial, tp2 = runs(iterations)
    tserial, tp3 = runs(iterations)
    tparallel = (tp1 + tp2 + tp3) / 3
    #    println("Time taken for CPU = $tserial")
    #    println("Time taken for GPU = $tparallel")
    println("10^$(pwr) options:")
    @printf("    Speedup = %.2f\n", tserial / tparallel)
    @printf("    CPU rate = 10^%.2f opts/sec\n", log10(iterations / tserial))
    @printf("    GPU rate = 10^%.2f opts/sec\n", log10(iterations / tparallel))
end

for k = 3:8
    driver(k)
end
