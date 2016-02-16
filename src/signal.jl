#Signal Processing Functions

# Fourier Transforms
# TODO: Multidimensional
import Base: fft
function fft{T}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{Complex{T}}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{Complex{T}}(icxx"fft2($a);")
	end
end
		
function fft{T<:Complex}(a::AFAbstractArray{T})
	if ndims(a) == 1
		AFArray{T}(icxx"fft($a);")
	elseif ndims(a) == 2
		AFArray{T}(icxx"fft2($a);")
	end
end
