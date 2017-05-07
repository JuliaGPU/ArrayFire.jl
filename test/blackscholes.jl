#This example has been adopted from https://github.com/IntelLabs/ParallelAccelerator.jl/blob/master/examples/black-scholes/black-scholes.jl

function blackscholes_serial(sptprice,
                             strike,
                             rate,
                             volatility,
                             time)
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

    put1 = blackscholes_serial.(sptprice, initStrike, rate, volatility, time)
    put2 = blackscholes_serial.(sptprice_gpu, initStrike_gpu, rate_gpu, volatility_gpu, time_gpu)
    afeval(put2)
    @test sum(put1) ≈ sum(put2)
    put2
end

@inferred runs(10^3)
