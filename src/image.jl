### Image Processing

import Base: hist, scale

# Export methods

export  loadImage, 
        saveImage, 
        isImageIOAvailable, 
        colorspace, 
        gray2rgb, 
        rgb2gray, 
        rgb2hsv, 
        rgb2ycbcr, 
        ycbcr2rgb, 
        hsv2rgb,
        regions,
        SAT,
        bilateral,
        maxfilt,
        meanshift,
        medfilt,
        minfilt,
        sobel,
        histequal,
        resize,
        rotate,
        skew,
        transform, 
        transformCoordinates,
        translate,
        dilate,
        erode,
        dilate3d,
        erode3d,
        gaussiankernel

# Export constants

export  AF_GRAY,
        AF_RGB,
        AF_HSV,
        AF_YCbCr,
        AF_YCC_601,
        AF_YCC_709,
        AF_YCC_2020,
        AF_CONNECTIVITY_4,
        AF_CONNECTIVITY_8,
        AF_PAD_ZERO,
        AF_PAD_SYM,
        AF_INTERP_NEAREST,
        AF_INTERP_LINEAR,
        AF_INTERP_BILINEAR,
        AF_INTERP_CUBIC,
        AF_INTERP_LOWER

# Colorspaces and Constants

AF_GRAY = 0 
AF_RGB = 1  
AF_HSV = 2      
AF_YCbCr = 3     

AF_YCC_601 = 601  
AF_YCC_709 = 709  
AF_YCC_2020 = 2020  

AF_CONNECTIVITY_4 = 4
AF_CONNECTIVITY_8 = 8

AF_PAD_ZERO = 0
AF_PAD_SYM = 1

AF_INTERP_NEAREST = 0 
AF_INTERP_LINEAR = 1
AF_INTERP_BILINEAR = 2
AF_INTERP_CUBIC = 3
AF_INTERP_LOWER = 4

# Read and write images

function loadImage(path::AbstractString; color = false)
    out = new_ptr()
    af_load_image(out, path, color)
    AFArray{backend_eltype(out[])}(out[])
end

function saveImage(path::AbstractString, a::AFArray)
    af_save_image(path, a)
end    

function isImageIOAvailable()
    out = Base.Ref{Bool}(0)
    af_is_image_io_available(out)
    out[]
end

# Colorspace conversions

function colorspace(img::AFArray, to::Int, from::Int)
    out = new_ptr()
    af_colorspace(out, img, to, from)
    AFArray{backend_eltype(out[])}(out[])
end

function gray2rgb(img::AFArray; rFactor::Cdouble = 1., gFactor::Cdouble = 1., bFactor::Cdouble = 1.)
    out = new_ptr()
    af_gray2rgb(out, img, rFactor, gFactor, bFactor)
    AFArray{backend_eltype(out[])}(out[])
end 

function hsv2rgb(img::AFArray)
    out = new_ptr()
    af_hsv2rgb(out, img)
    AFArray{backend_eltype(out[])}(out[])
end

function rgb2gray(img::AFArray; rPercent::Cfloat = 1., gPercent::Cfloat = 1., bFactor::Cfloat = 1.)
    out = new_ptr()
    af_rgb2gray(out, img)
    AFArray{backend_eltype(out[])}(out[])
end

function rgb2hsv(img::AFArray)
    out = new_ptr()
    af_rgb2hsv(out, img)
    AFArray{backend_eltype(out[])}(out[])
end

function rgb2ycbcr(img::AFArray; std::Int = AF_YCC_601)
    out = new_ptr()
    af_rgb2ycbcr(out, img, std)
    AFArray{backend_eltype(out[])}(out[])
end

function ycbcr2rgb(img::AFArray; std::Int = AF_YCC_601)
    out = new_ptr()
    af_ycbcr2rgb(out, img, std)
    AFArray{backend_eltype(out[])}(out[])
end

# Regions

function regions(a::AFArray; conn = AF_CONNECTIVITY_4, typ = Float32)
    out = new_ptr()
    af_regions(out, a, conn, typ)
    AFArray{typ}(out[])
end 

# Filters

function SAT(a::AFArray)
    out = new_ptr()
    af_sat(out, a)
    AFArray{backend_eltype(out[])}(out[])
end

function bilateral(a::AFArray, spatial_sigma::Real, chromatic_sigma::Real; color = false)
    out = new_ptr()
    af_bilateral(out, a, spatial_sigma, chromatic_sigma, color)
    AFArray{backend_eltype(out[])}(out[])
end

function maxfilt(a::AFArray; wind_length = 3, wind_width = 3, edge_pad = AF_PAD_ZERO)
    out = new_ptr()
    af_maxfilt(out, a, wind_length, wind_width, edge_pad)
    AFArray{backend_eltype(out[])}(out[])
end

function meanshift(a::AFArray, spatial_sigma::Real, chromatic_sigma::Real, iter::Int; color = false)
    out = new_ptr()
    af_meanshift(out, a, spatial_sigma, chromatic_sigma, Cuint(iter), color)
    AFArray{backend_eltype(out[])}(out[])
end

function medfilt(a::AFArray; wind_length = 3, wind_width = 3, edge_pad = AF_PAD_ZERO)
    out = new_ptr()
    af_medfilt(out, a, wind_length, wind_width, edge_pad)
    AFArray{backend_eltype(out[])}(out[])
end

function minfilt(a::AFArray; wind_length = 3, wind_width = 3, edge_pad = AF_PAD_ZERO)
    out = new_ptr()
    af_minfilt(out, a, wind_length, wind_width, edge_pad)
    AFArray{backend_eltype(out[])}(out[])
end

function sobel(a::AFArray; ker_size = 3)
    dx, dy = new_ptr(), new_ptr()
    af_sobel_operator(dx, dy, a, Cuint(ker_size))
    AFArray{backend_eltype(dx[])}(dx[]), AFArray{backend_eltype(dy[])}(dy[])
end

# Histograms

function histequal(a::AFArray, hist::AFArray)
    out = new_ptr()
    af_hist_equal(out, a, hist)
    AFArray{backend_eltype(out[])}(out[])
end

function hist(a::AFArray, nbins::Int, minval::Cdouble, maxval::Cdouble)
    out = new_ptr()
    af_histogram(out, a, Cuint(nbins), minval, maxval)
    AFArray{backend_eltype(out[])}(out[])
end

# Image Transformations

function resize(a::AFArray, dim1::Int, dim2::Int; 
                method = AF_INTERP_NEAREST)
    out = new_ptr()
    af_resize(out, a, dim1, dim2, method)
    AFArray{backend_eltype(out[])}(out[])
end

function rotate(a::AFArray, theta::Real; crop::Bool = true, 
                method = AF_INTERP_NEAREST)
    out = new_ptr()
    af_rotate(out, a, theta, crop, method)
    AFArray{backend_eltype(out[])}(out[])
end 

function scale(a::AFArray, scale0::Real, scale1::Real; 
                dim1 = 0, dim2 = 0, method = AF_INTERP_NEAREST)
    out = new_ptr()
    af_scale(out, a, scale0, scale1, dim1, dim2, method)
    AFArray{backend_eltype(out[])}(out[])
end 

function skew(a::AFArray, skew1::Real, skew2::Real; 
                dim1 = 0, dim2 = 0, method = AF_INTERP_NEAREST, 
                inverse = true)
    out = new_ptr()
    af_skew(out, a, skew1, skew2, dim1, dim2, method, inverse)
    AFArray{backend_eltype(out[])}(out[])
end 

function transform(a::AFArray, transform::AFArray; 
                    dim1 = 0, dim2 = 0, method = AF_INTERP_NEAREST, 
                    inverse = true) 
    out = new_ptr()
    af_transform(out, a, transform, dim1, dim2, method, inverse)
    AFArray{backend_eltype(out[])}(out[])
end 

function transformCoordinates(a::AFArray, d1::Real, d2::Real)
    out = new_ptr()
    af_transform_coordinates(out, a, d1, d2)
    AFArray{backend_eltype(out[])}(out[])
end 

function translate(a::AFArray, trans1::Real, trans2::Real; 
                    dim1 = 0, dim2 =0, method = AF_INTERP_NEAREST)
    out = new_ptr()
    af_translate(out, a, trans1, trans2, dim1, dim2, method)
    AFArray{backend_eltype(out[])}(out[])
end 

# Morphological Operations

for (op, fn) in ((:dilate, :af_dilate), (:dilate3d, :af_dilate3d),
                (:erode, :af_erode), (:erode3d, :af_erode3d))

    @eval function ($op)(a::AFArray, mask::AFArray)
        out = new_ptr()
        eval($fn)(out, a, mask)
        AFArray{backend_eltype(out[])}(out[])
    end

end

# Gaussian Kernel

function gaussiankernel(rows::Int, cols::Int; sigma_r = 0., sigma_c = 0.)
    out = new_ptr()
    af_gaussian_kernel(out, rows, cols, sigma_r, sigma_c)
    AFArray{backend_eltype(out[])}(out[])
end
