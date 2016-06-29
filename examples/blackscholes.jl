#This example has been adopted from https://github.com/IntelLabs/ParallelAccelerator.jl/blob/master/examples/black-scholes/black-scholes.jl

using ArrayFire

function blackscholes_serial(sptprice::AbstractArray{Float32},
                           strike::AbstractArray{Float32},
                           rate::AbstractArray{Float32},
                           volatility::AbstractArray{Float32},
                           time::AbstractArray{Float32})
    logterm = log10( sptprice ./ strike)
    powterm = .5f0 .* volatility .* volatility
    den = volatility .* sqrt(time)
    d1 = (((rate .+ powterm) .* time) .+ logterm) ./ den
    d2 = d1 .- den
    NofXd1 = cndf2(d1)
    NofXd2 = cndf2(d2)
    futureValue = strike .* exp(- rate .* time)
    c1 = futureValue .* NofXd2
    call_ = sptprice .* NofXd1 .- c1
    put  = call_ .- futureValue .+ sptprice
end

@inline function cndf2(in::AbstractArray{Float32})
    out = 0.5f0 .+ 0.5f0 .* erf(0.707106781f0 .* in)
    return out
end

function run(iterations)
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

    tic()
    put1 = blackscholes_serial(sptprice, initStrike, rate, volatility, time)
    t1 = toq()
    println("Serial checksum: ", sum(put1))
    tic()
    put2 = blackscholes_serial(sptprice_gpu, initStrike_gpu, rate_gpu, volatility_gpu, time_gpu)
    put2host = Array(put2)
    t2 = toq()
    println("Parallel checksum: ", sum(put2))
    return t1, t2
end

function driver()
    srand(0)
    tic()
    iterations = 10^7
    blackscholes_serial(Float32[], Float32[], Float32[], Float32[], Float32[])
    blackscholes_serial(AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]), 
                        AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]), AFArray(Float32[1., 2.]))
    println("SELFPRIMED ", toq())
    tserial, tparallel = run(iterations)
    println("Time taken for CPU = $tserial")
    println("Time taken for GPU = $tparallel")
    println("Speedup = $(tserial / tparallel)")
    println("CPU rate = ", iterations / tserial, " opts/sec")
    println("GPU rate = ", iterations / tparallel, " opts/sec")
end
driver()
