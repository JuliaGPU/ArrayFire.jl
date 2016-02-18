
#Image
export loadImage, saveImage, rotate, scale, transform, skew, translate,
        regions, rgb2gray, gray2rgb, rgb2hsv, hsv2rgb, rgb2ycbcr, ycbcr2rgb

#Image Utils
"Load image as an AFArray"
function loadImage(path::AbstractString; color = false)
	return AFArray{Float32}(icxx"loadImage($path, $color);")
end

"Save AFArray as image"
function saveImage(path::AbstractString, a::AFAbstractArray)
	icxx"saveImage($path, $a);"
end

#Image Transformation
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

#Image labelling
regions{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::regions($a.as(b8));")

#Colorspace conversions
gray2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::gray2rgb($a);")
hsv2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::hsv2rgb($a);")
rgb2gray{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::rgb2gray($a);")
rgb2hsv{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::rgb2hsv($a);")
rgb2ycbcr{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::rgb2ycbcr($a);")
ycbcr2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::ycbcr2rgb($a);")

#Filters
export SAT, bilateral, maxfilt, medfilt, minfilt, sobel, meanShift

"Summed Area Tables"
SAT{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::sat($a);")

"Apply bilateral filter on image"
function bilateral{T}(a::AFAbstractArray{T}, spatial_sigma::Real, chromatic_sigma::Real)
    AFArray{T}(icxx"af::bilateral($a, $spatial_sigma, $chromatic_sigma)")
end

"For every pixel, find max value from a window"
maxfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::maxfilt($a);")

"For every pixel, find median value from a window"
medfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::medfilt($a);")

"For every pixel, find minimum value from a window"
minfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::minfilt($a);")

"Sobel operator on image"
sobel{T}(a::AFAbstractArray{T}) = AFArray{T}(icxx"af::sobel($a);")

"Meanshift filter"
function meanShift{T}(a::AFAbstractArray{T}, spatial_sigma::Real, chromatic_sigma::Real, iter::Integer = 100)
    AFArray{T}(icxx"af::meanShift($a, $spatial_sigma, $chromatic_sigma, $iter)")
end

#Histogram
import Base: hist
export histEqual
hist{T}(a::AFAbstractArray{T}, nbins::Integer) = AFArray{T}(icxx"af::histogram($a, $nbins);")
hist{T}(a::AFAbstractArray{T}, nbins::Integer, minval::Real, maxval::Real) = AFArray{T}(icxx"af::histogram($a, $nbins, $minval, $maxval);")

"Histogram equalization of input image"
histEqual{T}(a::AFAbstractArray{T}, hist::AFAbstractArray{T}) = AFArray{T}(icxx"af::histEqual($a, $hist)")

#Morphological ops
export dilate, dilate3, erode, erode3

dilate{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::dilate($a, $mask);")
dilate3{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::dilate3($a, $mask);")
erode{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::erode($a, $mask);")
erode3{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::erode3($a, $mask);")

#Computer Vision

#Constants
const AF_SAD = icxx"AF_SAD;"
const AF_ZSAD = icxx"AF_ZSAD;"
const AF_LSAD = icxx"AF_LSAD;"
const AF_SSD = icxx"AF_SSD;"
const AF_ZSSD = icxx"AF_ZSSD;"
const AF_LSSD = icxx"AF_LSSD;"
const AF_NCC = icxx"AF_NCC;"
const AF_ZNCC = icxx"AF_ZNCC;"
const AF_SHD = icxx"AF_SHD;"

export matchTemplate, dog

function matchTemplate{T}(searchImg::AFAbstractArray{T}, template::AFAbstractArray{T}, matchType = AF_SAD)
    AFArray{T}(icxx"af::matchTemplate($searchImg, $template, $matchType);")
end

DiffOfGaussians{T}(img::AFAbstractArray{T}, radius1::Integer, radius2::Integer) = AFArray{T}(icxx"dog($img, $radius1, $radius2);")
