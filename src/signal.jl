#Signal Processing Functions

#Consts
AF_CONV_EXPAND = icxx"AF_CONV_EXPAND;"

# Fourier Transforms
import Base: fft, ifft
function fft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(af_fft(a))
	elseif ndims(a) == 2
		AFArray{Complex{T}}(af_fft2(a))
    elseif ndims(a) == 3
        AFArray{Complex{T}}(af_fft3(a))
    end
end

function fft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(af_fft(a))
	elseif ndims(a) == 2
		AFArray{T}(af_fft2(a))
	elseif ndims(a) == 3
        AFArray{T}(af_fft3(a))
    end
end

function ifft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(af_ifft(a))
	elseif ndims(a) == 2
		AFArray{Complex{T}}(af_ifft2(a))
    elseif ndims(a) == 3
        AFArray{Complex{T}}(af_ifft3(a))
    end
end

function ifft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(af_ifft(a))
	elseif ndims(a) == 2
		AFArray{T}(af_ifft2(a))
    elseif ndims(a) == 3
        AFArray{T}(af_ifft3(a))
    end
end
		
#FIR Filter
fir{T}(a::AFAbstractArray{T}, x::AFAbstractArray{T}) = AFArray{T}(af_fir(a,x))

#IIR Filter
iir{T}(ff::AFAbstractArray{T}, fb::AFAbstractArray{T}, a::AFAbstractArray{T}) = AFArray{T}(af_iir(ff, fb, a))

#Convolutions

import Base: conv, conv2

conv{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_fftConvolve(a, b, AF_CONV_EXPAND))
conv2{T}(a::AFAbstractArray{T}, b::AFAbstractArray{T}) = AFArray{T}(af_fftConvolve2(a,b, AF_CONV_EXPAND))
