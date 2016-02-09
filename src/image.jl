
#Image
export loadImage, saveImage, rotate, scale, transform, skew, translate

"Load image as an AFArray"
function loadImage(path::AbstractString; color = false)
	return AFArray{Float32}(icxx"loadImage($path, $color);")
end

"Save AFArray as image"
function saveImage(path::AbstractString, a::AFAbstractArray)
	icxx"saveImage($path, $a);"
end

"Rotate image by theta radians"
function rotate{T}(a::AFAbstractArray{T}, theta)
	AFArray{T}(icxx"rotate($a, $theta);")
end

"Scale image in either dimension"
function scale{T}(a::AFAbstractArray{T}, scale1::Real, scale2::Real)
	AFArray{T}(icxx"scale($a, $scale1, $scale2);")
end

"Skew an image in either dimension by angles"
function skew{T}(a::AFAbstractArray{T}, skew1::Real, skew2::Real)
	AFArray{T}(icxx"skew($a, $skew1, $skew2);")
end

"Translate a given image in either dimension"
function translate{T}(a::AFAbstractArray{T}, trans1::Real, trans2::Real)
	AFArray{T}(icxx"translate($a, $trans1, $trans2);")
end

function transform{T}(a::AFAbstractArray{T}, b::AFAbstractArray)
	AFArray{T}(icxx"transform($a, $b);")
end
