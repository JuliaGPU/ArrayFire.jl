for T in (Float32, Float64)
    a1 = rand(Complex{T}, 10)
    b1 = fft(a1)
    af1 = AFArray(a1)
    bf1 = fft(af1, 1., 0)
    cf1 = ifft(bf1, 1 / length(bf1), 0)
    @test eltype(bf1) == eltype(cf1) == Complex{T}
    @test Array(bf1) ≈ b1
    @test Array(cf1) ≈ a1

    fft!(af1, 1.0)
    @test eltype(af1) == Complex{T}
    @test Array(af1) ≈ b1

    ifft!(af1, 1 / length(af1))
    @test eltype(af1) == Complex{T}
    @test Array(af1) ≈ a1


    a2 = rand(Complex{T}, 10, 11)
    b2 = fft(a2)
    af2 = AFArray(a2)
    bf2 = fft2(af2, 1., 0, 0)
    cf2 = ifft2(bf2, 1 / length(bf2), 0, 0)
    @test eltype(bf2) == eltype(cf2) == Complex{T}
    @test Array(bf2) ≈ b2
    @test Array(cf2) ≈ a2

    fft2!(af2, 1.0)
    @test eltype(af2) == Complex{T}
    @test Array(af2) ≈ b2

    ifft2!(af2, 1 / length(af2))
    @test eltype(af2) == Complex{T}
    @test Array(af2) ≈ a2


    a3 = rand(Complex{T}, 10, 11, 12)
    b3 = fft(a3)
    af3 = AFArray(a3)
    bf3 = fft3(af3, 1., 0, 0, 0)
    cf3 = ifft3(bf3, 1 / length(bf3), 0, 0, 0)
    @test eltype(bf3) == eltype(cf3) == Complex{T}
    @test Array(bf3) ≈ b3
    @test Array(cf3) ≈ a3

    fft3!(af3, 1.0)
    @test eltype(af3) == Complex{T}
    @test Array(af3) ≈ b3

    ifft3!(af3, 1 / length(af3))
    @test eltype(af3) == Complex{T}
    @test Array(af3) ≈ a3


    ar1 = rand(T, 10)
    br1 = rfft(ar1)
    arf1 = AFArray(ar1)
    brf1 = fft_r2c(arf1, 1., 0)
    crf1 = fft_c2r(brf1, 1 / length(arf1), isodd(size(arf1, 1)))
    @test eltype(brf1) == Complex{T}
    @test eltype(crf1) == T
    @test Array(brf1) ≈ br1
    @test Array(crf1) ≈ ar1


    ar2 = rand(T, 10, 11)
    arf2 = AFArray(ar2)
    brf2 = fft2_r2c(arf2, 1., 0, 0)
    crf2 = fft2_c2r(brf2, 1 / length(arf2), isodd(size(arf2, 1)))
    @test eltype(brf2) == Complex{T}
    @test eltype(crf2) == T
    @test Array(brf2) ≈ rfft(ar2)
    @test Array(crf2) ≈ ar2


    ar3 = rand(T, 11, 12, 10)
    arf3 = AFArray(ar3)
    brf3 = fft3_r2c(arf3, 1., 0, 0, 0)
    crf3 = fft3_c2r(brf3, 1 / length(arf3), isodd(size(arf3, 1)))
    @test eltype(brf3) == Complex{T}
    @test eltype(crf3) == T
    @test Array(brf3) ≈ rfft(ar3)
    @test Array(crf3) ≈ ar3
end
