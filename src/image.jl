
#Image
export loadImage, saveImage, rotate, scale, transform, skew, translate,
        regions, rgb2gray, gray2rgb, rgb2hsv, hsv2rgb, rgb2ycbcr, ycbcr2rgb

#Image Utils
"Load image as an AFArray"
function loadImage(path::AbstractString; color = false)
	return AFArray{Float32}(af_loadImage(path, color))
end

"Save AFArray as image"
function saveImage(path::AbstractString, a::AFAbstractArray)
    af_saveImage(path, a)
end

#Image Transformation
"Rotate image by theta radians"
function rotate{T}(a::AFAbstractArray{T}, theta)
	AFArray{T}(af_rotate(a,theta))
end

"Scale image in either dimension"
function scale{T}(a::AFAbstractArray{T}, scale1::Real, scale2::Real)
	AFArray{T}(af_scale(a,scale1,scale2))
end

"Skew an image in either dimension by angles"
function skew{T}(a::AFAbstractArray{T}, skew1::Real, skew2::Real)
	AFArray{T}(af_skew(a,b,c))
end

"Translate a given image in either dimension"
function translate{T}(a::AFAbstractArray{T}, trans1::Real, trans2::Real)
	AFArray{T}(af_translate(a,trans1, trans2))
end

function transform{T}(a::AFAbstractArray{T}, b::AFAbstractArray)
	AFArray{T}(af_transform(a,b))
end

#Image labelling
regions{T}(a::AFAbstractArray{T}) = AFArray{T}(af_regions(a))

#Create Gaussian Kernel
gaussiankernel(rows::Integer, cols::Integer) = AFArray{Float32}(af_gaussianKernel(rows,cols))

#Colorspace conversions
gray2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(af_gray2rgb(a))
hsv2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(af_hsv2rgb(a))
rgb2gray{T}(a::AFAbstractArray{T}) = AFArray{T}(af_rgb2gray(a))
rgb2hsv{T}(a::AFAbstractArray{T}) = AFArray{T}(af_rgb2hsv(a))
rgb2ycbcr{T}(a::AFAbstractArray{T}) = AFArray{T}(rgb2ycbcr(a))
ycbcr2rgb{T}(a::AFAbstractArray{T}) = AFArray{T}(ycbcr2rgb(a))
colorspace{T}(a::AFAbstractArray{T}, from::Cxx.CppEnum{:af_cspace_t}, to::Cxx.CppEnum{:af_cspace_t}) = af_colorspace(a, to, from)

#Filters
export SAT, bilateral, maxfilt, medfilt, minfilt, sobel, meanShift

"Summed Area Tables"
SAT{T}(a::AFAbstractArray{T}) = AFArray{T}(af_SAT(a))

"Apply bilateral filter on image"
function bilateral{T}(a::AFAbstractArray{T}, spatial_sigma::Real, chromatic_sigma::Real)
    AFArray{T}(af_bilateral(a, spatial_sigma, chromatic_sigma))
end

"For every pixel, find max value from a window"
maxfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(af_maxfilt(a))

"For every pixel, find median value from a window"
medfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(af_medfilt(a))

"For every pixel, find minimum value from a window"
minfilt{T}(a::AFAbstractArray{T}) = AFArray{T}(af_minfilt(a))

"Sobel operator on image"
sobel{T}(a::AFAbstractArray{T}) = AFArray{T}(af_sobel(a))

"Meanshift filter"
function meanShift{T}(a::AFAbstractArray{T}, spatial_sigma::Real, chromatic_sigma::Real, iter::Integer = 100)
    AFArray{T}(af_meanShift(a, spatial_sigma, chromatic_sigma, iter))
end

#Histogram
import Base: hist
export histEqual
hist{T}(a::AFAbstractArray{T}, nbins::Integer) = AFArray{T}(af_histogram(a, nbins))
function hist{T}(a::AFAbstractArray{T}, nbins::Integer, minval::Real, maxval::Real)
    AFArray{T}(af_histogram(a, nbins, minval, maxval))
end

"Histogram equalization of input image"
histEqual{T}(a::AFAbstractArray{T}, hist::AFAbstractArray{T}) = AFArray{T}(af_histEqual(a, hist))

#Morphological ops
export dilate, dilate3, erode, erode3

dilate{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::dilate($a, $mask);")
dilate3{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::dilate3($a, $mask);")
erode{T}(a::AFAbstractArray{T}, mask::AFAbstractArray) = AFArray{T}(icxx"af::erode($a, $mask);")
erode3{T}(a::AFAbstractArray{T}, mask::AFAbstractArray{T}) = AFArray{T}(icxx"af::erode3($a, $mask);")

#Computer Vision

#Constants
export matchTemplate, DiffOfGaussians

function matchTemplate{T}(searchImg::AFAbstractArray{T}, template::AFAbstractArray{T}, matchType = AF_SAD)
    AFArray{T}(af_matchTemplate(searchImg, template, matchType))
end

DiffOfGaussians{T}(img::AFAbstractArray{T}, radius1::Integer, radius2::Integer) = AFArray{T}(dog(img, radius1, radius2))

import Base:show 
show(io::IO, f::AFFeatures) = show(io, f.feat)
show(io::IO, f::vcpp"af::features") = print(io, "ArrayFire Feature")

function getX(f::AFFeatures) 
    out = af_getX(f)
    AFArray{backend_eltype(out)}(out)
end
function getY(f::AFFeatures) 
    out = af_getY(f)
    AFArray{backend_eltype(out)}(out)
end
function getScore(f::AFFeatures) 
    out = af_getScore(f)
    AFArray{backend_eltype(out)}(out)
end 
function getSize(f::AFFeatures) 
    out = af_getSize(f)
    AFArray{backend_eltype(out)}(out)
end 
function getOrientation(f::AFFeatures)
    out = af_getOrientation(f)
    AFArray{backend_eltype(out)}(out)
end

function ORB(img::AFAbstractArray; fast_thr = 20.0, max_feat = 400, scl_fctr = 1.5, levels = 4, blur_img = false)
    out_features = AFFeatures()
    out_desc = AFArray()
    af_orb(out_features, out_desc, img, fast_thr, max_feat, scl_fctr, levels, blur_img)
    AFFeatures(out_features), AFArray{backend_eltype(out_desc)}(out_desc)
end 

function sift(img::AFAbstractArray; n_layers = 3, constant_thr = 0.04, 
                edge_thr = 10.0, init_sigma = 1.6, double_input = true, 
                intensity_scale = 0.00390625, feature_ratio = 0.05)
    out_features = AFFeatures()
    out_desc = AFArray()
    af_sift(out_features, out_desc, img, n_layers, constant_thr, 
            edge_thr, init_sigma, double_input, intensity_scale, feature_ratio)
    AFFeatures(out_features), AFArray{backend_eltype(out_desc)}(out_desc)
end

function fast(img::AFAbstractArray; thr = 20.0, arc_length = 9, non_max = true, feature_ratio = 0.05, edge = 3)
    AFFeatures(af_fast(img, thr, arc_length, non_max, feature_ratio, edge))
end

function harris(img::AFAbstractArray; max_corners = 500, min_response = 1e5, sigma = 1.0, block_size = 0, k_thr = 0.04)
    AFFeatures(af_harris(img, max_corners, min_response, sigma, block_size, k_thr))
end

function susan(img::AFAbstractArray; radius = 3, diff_thr = 32.0, geom_thr = 10.0, feature_ratio = 0.05, edge = 3)
    AFFeatures(af_susan(img, radius, diff_thr, geom_thr, feature_ratio, edge))
end

function hammingMatcher(query::AFAbstractArray, train::AFAbstractArray; dim = 1, n_dist = 1)
    out_idx = AFArray()
    out_dist = AFArray()
    af_hammingMatcher(out_idx, out_dist, query, train, dim-1, dist)
    AFArray{backend_eltype(out_idx)}(out_idx), AFArray{backend_eltype(out_dist)}(out_dist)
end

function nearnestNeighbour(query::AFAbstractArray, train::AFAbstractArray; dim = 1, n_dist = 1, dist_type = AF_SSD)
    out_idx = AFArray()
    out_dist = AFArray()
    af_nearestNeighbour(out_idx, out_dist, query, train, dim, n_dist, dist_type)
    AFArray{backend_eltype(out_idx)}(out_idx), AFArray{backend_eltype(out_dist)}(out_dist)
end
