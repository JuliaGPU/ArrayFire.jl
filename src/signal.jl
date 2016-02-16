#Signal Processing Functions

# Fourier Transforms
import Base: fft, ifft
function fft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{Complex{T}}(icxx"fft2($a);")
    elseif ndims(a) == 3
        AFArray{Complex{T}}(icxx"fft3($a);")
    end
end

function fft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{T}(icxx"fft2($a);")
	elseif ndims(a) == 3
        AFArray{T}(icxx"fft3($a);")
    end
end

function ifft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(icxx"ifft($a);")
	elseif ndims(a) == 2
		AFArray{Complex{T}}(icxx"ifft2($a);")
    elseif ndims(a) == 3
        AFArray{Complex{T}}(icxx"ifft3($a);")
    end
end

function ifft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(icxx"ifft($a);")
	elseif ndims(a) == 2
		AFArray{T}(icxx"ifft2($a);")
    elseif ndims(a) == 3
        AFArray{T}(icxx"ifft3($a);")
    end
end
		
#FIR Filter
fir{T}(a::AFAbstractArray{T}, x::AFAbstractArray{T}) = AFArray{T}(icxx"af::fir($a, $x);")

#IIR Filter
iir{T}(ff::AFAbstractArray{T}, fb::AFAbstractArray{T}, a::AFAbstractArray{T}) = AFArray{T}(icxx"af::iir($b, $a, $x);")
