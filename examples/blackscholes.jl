#This example has been adopted from https://github.com/IntelLabs/ParallelAccelerator.jl/blob/master/examples/black-scholes/black-scholes.jl

function blackscholes_serial(sptprice::Array{Float64,1},
                           strike::Array{Float64,1},
                           rate::Array{Float64,1},
                           volatility::Array{Float64,1},
                           time::Array{Float64,1})
    logterm = log10(sptprice ./ strike)
    powterm = .5 .* volatility .* volatility
    den = volatility .* sqrt(time)
    d1 = (((rate .+ powterm) .* time) .+ logterm) ./ den
    d2 = d1 .- den
    NofXd1 = cndf2(d1)
    NofXd2 = cndf2(d2)
    futureValue = strike .* exp(- rate .* time)
    c1 = futureValue .* NofXd2
    call = sptprice .* NofXd1 .- c1
    put  = call .- futureValue .+ sptprice
end

function blackscholes_gpu(sptprice::AFArray{Float64,1},
                           strike::AFArray{Float64,1},
                           rate::AFArray{Float64,1},
                           volatility::AFArray{Float64,1},
                           time::AFArray{Float64,1})
    logterm = log10(sptprice ./ strike)
    powterm = .5 .* volatility .* volatility
    den = volatility .* sqrt(time)
    d1 = (((rate .+ powterm) .* time) .+ logterm) ./ den
    d2 = d1 .- den
    NofXd1 = cndf2_gpu(d1)
    NofXd2 = cndf2_gpu(d2)
    futureValue = strike .* exp(- rate .* time)
    c1 = futureValue .* NofXd2
    call = sptprice .* NofXd1 .- c1
    put  = call .- futureValue .+ sptprice
end

function blackscholes_devec(sptprice::Vector{Float64}, strike::Vector{Float64}, rate::Vector{Float64}, volatility::Vector{Float64}, time::Vector{Float64})
    sqt = sqrt(time)
    put = similar(strike)
    for i = 1:size(sptprice, 1)
        logterm = log10(sptprice[i] / strike[i])
        powterm = 0.5 * volatility[i] * volatility[i]
        den = volatility[i] * sqt[i]
        d1 = (((rate[i] + powterm) * time[i]) + logterm) / den
        d2 = d1 - den
        NofXd1 = 0.5 + 0.5 * erf(0.707106781 * d1)
        NofXd2 = 0.5 + 0.5 * erf(0.707106781 * d2)
        futureValue = strike[i] * exp(-rate[i] * time[i])
        c1 = futureValue * NofXd2
        call = sptprice[i] * NofXd1 - c1
        put[i] = call - futureValue + sptprice[i]
    end
    put
end

@inline function cndf2(in::Array{Float64,1})
    out = 0.5 .+ 0.5 .* erf(0.707106781 .* in)
    return out
end

@inline function cndf2_gpu(in::AFArray{Float64,1})
    out = 0.5 .+ 0.5 .* erf(0.707106781 .* in)
    return out
end

function run(iterations)
    sptprice   = Float64[ 42.0 for i = 1:iterations ]
    initStrike = Float64[ 40.0 + (i / iterations) for i = 1:iterations ]
    rate       = Float64[ 0.5 for i = 1:iterations ]
    volatility = Float64[ 0.2 for i = 1:iterations ]
    time       = Float64[ 0.5 for i = 1:iterations ]

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
    put2 = blackscholes_gpu(sptprice_gpu, initStrike_gpu, rate_gpu, volatility_gpu, time_gpu)
    t2 = toq()
    println("Parallel checksum: ", sum(put2))
    tic()
    put3 = blackscholes_devec(sptprice, initStrike, rate, volatility, time)
    t3 = toq()
    println("Devectorized checksum: ", sum(put2))
    return t1, t2, t3
end

function driver()
    srand(0)
    tic()
    iterations = 10^7
    blackscholes_serial(Float64[], Float64[], Float64[], Float64[], Float64[])
    blackscholes_devec(Float64[], Float64[], Float64[], Float64[], Float64[])
    blackscholes_gpu(AFArray(Float64[1.]), AFArray(Float64[1.]), AFArray(Float64[1.]), AFArray(Float64[1.]), AFArray(Float64[1.]))
    println("SELFPRIMED ", toq())
    tserial, tparallel, tdevec = run(iterations)
    println("Time taken for vectorized = $tserial")
    println("Time taken for devectorized = $tdevec")
    println("Time taken for parallel = $tparallel")
    println("Speedup over vectorized = $(tserial / tparallel)")
    println("Speedup over devectorized = $(tdevec / tparallel)")
    println("Vectorized  rate = ", iterations / tserial, " opts/sec")
    println("Devectorized  rate = ", iterations / tdevec, " opts/sec")
    println("Parallel rate = ", iterations / tparallel, " opts/sec")
end
driver()
