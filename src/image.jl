### Image Processing

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
        sobel

# Export constants

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
    out = new_ptr()
    af_sobel_operator(dx, dy, img, Cuint(ker_size))
    AFArray{backend_eltype(out[])}(out[])
end
