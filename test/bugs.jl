a = rand(Float32, 10, 10)
ad = AFArray(a)

# Indexing - issue #96
let
    b = Float32[1.,2.,3.]
    bd = AFArray(b)
    ind = AFArray([false, true, true])
    @test Array(bd[ind]) == Float32[2., 3.]
    ind[3] = false
    @test Array(bd[ind]) == Float32[2.]
end

# Broadcast
if VERSION >= v"0.5.0"
    @test sin(ad) == sin.(ad)
end

# Sign - issue #109
if VERSION >= v"0.5.0"
    let
        a = randn(Float32, 10)
        ad = AFArray(a)
        @test Array(signbit.(ad)) == signbit.(a)
        @test Array(sign.(ad)) == sign.(a)
    end
end

# Indexing - issue #115
let
    x = rand(Float32, 3,3)
    xd = AFArray(x)
    y = x[:, [1,3]]
    yd = xd[:, [1,3]]
    @test x ≈ Array(xd)
end

# Return types of sum reduction
let
    for T in (Float32, Complex{Float32}, Int32,
              Int64, UInt32, UInt8, UInt64, Bool)
        s = AFArray(rand(T, 10))
        if T == UInt8
            @test typeof(sum(s)) == UInt32
        elseif T == Bool
            @test typeof(sum(s)) == Int64
        else
            @test typeof(sum(s)) == T
        end
    end
end

# Issue #125
let
    a = rand(Float32, 5, 5)
    ad = AFArray(a)
    if VERSION >= v"0.6.0"
        for f in (>, >=, <, <=, ==)
            for val in (0.5f0, 1)
                b = broadcast(f, val, a)
                bd = broadcast(f, val, ad)
                @test b ≈ Array(bd)
                b = broadcast(f, a, val)
                bd = broadcast(f, ad, val)
                @test b ≈ Array(bd)
            end
        end
    else
        for f in (.>, .>=, .<, .<=, .==)
            for val in (0.5f0, 1)
                b = f(val, a)
                bd = f(val, ad)
                @test b ≈ Array(bd)
                b = f(a, val)
                bd = f(ad, val)
                @test b ≈ Array(bd)
            end
        end
    end
end

# Issue #131
let
    set_device(0)
    for i = 1:10
        a = rand(AFArray{Float32}, 10)
        a + 1
    end
    sync(0)
end

# Issue #133
let
    for sz in ((10,), (10, 10), (10, 10, 10))
        a = rand(Complex{Float32}, sz...)
        ad = AFArray(a)
        fft!(a)
        fft!(ad)
        @test a ≈ Array(ad)
    end
end

# Issue #143
let
    AD = AFArray(ones(Float32,3,3)+im*ones(Float32,3,3))
    AH = Array(AD)
    @test Array(abs(AD)) ≈ abs.(AH)
end

# Issue #103
let
    a = rand(AFArray{Complex64}, 2, 2)
    a[1,1]  = Float32(1)+Float32(2)im
    @test a[1,1] == 1.0f0 + 2.0f0im
end

# Issue https://github.com/gaika/ArrayFire.jl/issues/26
@test im * ad == complex(zeros(ad), ad)
