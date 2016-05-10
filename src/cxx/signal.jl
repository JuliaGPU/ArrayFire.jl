#Signal Processing Functions

#Consts
AF_CONV_EXPAND = icxx"AF_CONV_EXPAND;"

# Fourier Transforms
import Base: fft, ifft
fft{T}(a::AFVector{T}) = AFArray{Complex{T}}(af_fft(a))
fft{T}(a::AFMatrix{T}) = AFArray{Complex{T}}(af_fft2(a))
fft{T}(a::AFAbstractArray{T,3}) = AFArray{Complex{T}}(af_fft3(a))

fft{T<:Complex}(a::AFVector{T}) = AFArray{T}(af_fft(a))
fft{T<:Complex}(a::AFMatrix{T}) = AFArray{T}(af_fft2(a))
fft{T<:Complex}(a::AFAbstractArray{T,3}) = AFArray{T}(af_fft3(a))

ifft{T}(a::AFVector{T}) = AFArray{Complex{T}}(af_ifft(a))
ifft{T}(a::AFMatrix{T}) = AFArray{Complex{T}}(af_ifft2(a))
ifft{T}(a::AFAbstractArray{T,3}) = AFArray{Complex{T}}(af_ifft3(a))

ifft{T<:Complex}(a::AFVector{T}) = AFArray{T}(af_ifft(a))
ifft{T<:Complex}(a::AFMatrix{T}) = AFArray{T}(af_ifft2(a))
ifft{T<:Complex}(a::AFAbstractArray{T,3}) = AFArray{T}(af_ifft3(a))
		
#FIR Filter
fir{T}(a::AFAbstractArray{T}, x::AFAbstractArray{T}) = AFArray{T}(af_fir(a,x))

#IIR Filter
iir{T}(ff::AFAbstractArray{T}, fb::AFAbstractArray{T}, a::AFAbstractArray{T}) = AFArray{T}(af_iir(ff, fb, a))

#Convolutions

import Base: conv, conv2

conv{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_fftConvolve(a, b, AF_CONV_EXPAND))
conv2{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_fftConvolve2(a,b, AF_CONV_EXPAND))
