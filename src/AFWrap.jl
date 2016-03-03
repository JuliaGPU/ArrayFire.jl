#This wraps all the basic ArrayFire routines as is.

#Math functions
import Base: +, -

# Resolve conflicts
+(x::AFAbstractArray{Bool},y::Bool) = AFArray{Bool}(@cxx +(x.array,y))
+(y::Bool,x::AFAbstractArray{Bool}) = AFArray{Bool}(@cxx +(y,x.array))
-(x::AFAbstractArray{Bool},y::Bool) = AFArray{Bool}(@cxx -(x.array,y))
-(y::Bool,x::AFAbstractArray{Bool}) = AFArray{Bool}(@cxx -(y,x.array))


for (op,cppop) in ((:+,:+),(:(.+),:+),(:-,:-),(:(.-),:-),(:.*,:*),(:./,:/),(:.>>,:>>),(:.<<,:<<))
    @eval function Base.($(quot(op))){T,S}(x::AFAbstractArray{T}, y::AFAbstractArray{S})
        a1 = x.array
        a2 = y.array
        AFArray{af_promote(T,S)}(@cxx ($(cppop))(a1,a2))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(x::AFAbstractArray{T}, y::S)
        a = x.array
        # This is special behavior hardcoded in arrayfire for real floats
        ST = S <: AbstractFloat ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(a, y))
    end
    @eval function Base.($(quot(op))){T,S<:Number}(y::S, x::AFAbstractArray{T})
        a = x.array
        # This is special behavior hardcoded in arrayfire for real floats
        ST = S <: AbstractFloat ? T : S
        AFArray{af_promote(T,ST)}(@cxx ($(cppop))(y, a))
    end
end

#Functions to generate arrays
af_constant(a::Real,b::vcpp"af::dim4", c::Cxx.CppEnum{:af_dtype}) = icxx"af::constant($a,$b,$c);"
af_randu(a::vcpp"af::dim4", b::Cxx.CppEnum{:af_dtype}) = icxx"af::randu($a,$b);"
af_randn(a::vcpp"af::dim4", b::Cxx.CppEnum{:af_dtype}) = icxx"af::randn($a,$b);"
af_identity(a::vcpp"af::dim4", b::Cxx.CppEnum{:af_dtype}) = icxx"af::identity($a,$b);"
af_diag(a::AFAbstractArray, b::Integer) = icxx"af::diag($a,$b);"
af_getSeed() = icxx"af::getSeed();"
af_range(a::vcpp"af::dim4", b::Integer, c::Cxx.CppEnum{:af_dtype}) = icxx"af::range($a,$b,$c);"
af_setSeed(a::Integer) = icxx"af::setSeed($a);"
af_iota(a::vcpp"af::dim4", b::vcpp"af::dim4", c::Cxx.CppEnum{:af_dtype}) = icxx"af::iota($a,$b,$c);"

#Function to modify or reorg arrays
af_tile(a,b) = icxx"af::tile($a,$b);"
af_join(a,b,c) = icxx"af::join($a,$b,$c);"
af_join(a,b,c,d) = icxx"af::join($a,$b,$c,$d);"
af_join(a,b,c,d,el) = icxx"af::join($a,$b,$c,$d,$el);"
af_moddims(a,b) = icxx"af::moddims($a,$b);"
af_reorder(a,b,c,d,el) = icxx"af::reorder($a,$b,$c,$d,$el);"
af_replace(a,b,c) = icxx"af::replace($a,$b,$c);"
af_select(a,b) = icxx"af::select($a,$b,$c);"
af_shift(a,b,c,d,el) = icxx"af::shift($a,$b,$c,$d,$el);"

#Some more array methods
af_complex_conj(a::AFAbstractArray) = icxx"$a.H();"
af_as(a::AFAbstractArray, to::Cxx.CppEnum{:af_dtype}) = icxx"(af::array)($a.as($to));"
af_copy(a::AFAbstractArray) = icxx"$a.copy();"
af_bytes(a::AFAbstractArray) = icxx"$a.bytes();"
af_isempty(a::AFAbstractArray) = icxx"$a.isempty();"
af_isscalar(a::AFAbstractArray) = icxx"$a.isscalar();"
af_isvector(a::AFAbstractArray) = icxx"$a.isvector();"
af_isrow(a::AFAbstractArray) = icxx"$a.isrow();"
af_iscolumn(a::AFAbstractArray) = icxx"$a.iscolumn();"
af_iscomplex(a::AFAbstractArray) = icxx"$a.iscomplex();"
af_isreal(a::AFAbstractArray) = icxx"$a.isreal();"
af_isdouble(a::AFAbstractArray) = icxx"$a.isdouble();"
af_issingle(a::AFAbstractArray) = icxx"$a.issingle();"
af_isrealfloating(a::AFAbstractArray) = icxx"$a.isrealfloating();"
af_isfloating(a::AFAbstractArray) = icxx"$a.isfloating();"
af_isinteger(a::AFAbstractArray) = icxx"$a.isinteger();"
af_isbool(a::AFAbstractArray) = icxx"$a.isbool();"

#Helper functions
af_isNaN(a::AFAbstractArray) = icxx"af::isNaN($a);"
af_isinf(a::AFAbstractArray) = icxx"af::isInf($a);"
af_iszero(a::AFAbstractArray) = icxx"af::iszero($a);"

#ArrayFire info functions
af_info() = icxx"af::info();"
af_isDoubleAvailable() = icxx"af::isDoubleAvailable();"
af_sync() = icxx"af::sync();"
af_sync(a) = icxx"af::sync($a);"
af_getDevice() = icxx"af::getDevice();"
af_setDevice(a) = icxx"af::setDevice($a);"

#Loading and Saving Arrays
af_saveArray(a::AFAbstractArray, b::AbstractString, c::AbstractString) = icxx"af::saveArray($a,$b,$c);"
af_readArray(a::AbstractString, b::AbstractString) = icxx"af::readArray($a,$b);"

#Numeric functions
af_abs(a::AFAbstractArray) = icxx"af::abs($a);"
af_max(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"af::max($a, $b);"
af_min(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"af::min($a, $b);"
af_arg(a::AFAbstractArray) = icxx"af::arg($a);"
af_ceil(a::AFAbstractArray) = icxx"af::cell($a);"
af_floor(a::AFAbstractArray) = icxx"af::floor($a);"
af_hypot(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::hypot($a, $b);"
af_mod(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::mod($a, $b);"
of_rem(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::rem($a, $b);"
af_round(a::AFAbstractArray) = icxx"af::round($a);"
af_sign(a::AFAbstractArray) = icxx"af::sign($a);"
af_trunc(a::AFAbstractArray) = icxx"af_trunc($a);"

#Logical ops
af_equals(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a == $b;"
af_gt(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a > $b;"
af_ge(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a >= $b;"
af_bitand(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a & $b;"
af_and(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a && $b;"
af_bitor(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a | $b;"
af_bitxor(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a ^ $b;"
af_le(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a <= $b;"
af_lt(a::Union{Real,AFAbstractArray}, b::Union{Real,AFAbstractArray}) = icxx"$a < $b;"
af_neg(a::AFAbstractArray) = icxx"-$a;"
af_neq(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a != $b;"
af_not(a::AFAbstractArray) = icxx"!$a;"
af_or(a::AFAbstractArray, b::AFAbstractArray) = icxx"$a || $b;"

#Complex functions
af_complex(a, b) = icxx"af::complex($a,$b);"
af_complex(a) = icxx"af::complex($a);"
af_conjg(a) = icxx"af::conjg($a);"
af_imag(a) = icxx"af::imag($a);"
af_real(a) = icxx"af::real($a);"

#Expontential functions
af_pow(a,b) = icxx"af::pow($a, $b);"
af_root(a,b) = icxx"af::root($a, $b);"

#Linear algebra
const AF_MAT_NONE = icxx"AF_MAT_NONE;"
const AF_MAT_CTRANS = icxx"AF_MAT_CTRANS;"

af_dot(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::dot($a, $b);"

af_matmul(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::matmul($a, $b);"
af_matmul3(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray) = icxx"af::matmul($a, $b, $c);"
af_matmul4(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::AFAbstractArray) = icxx"af::matmul($a, $b, $c, $d);"
af_matmul_flags(a::AFAbstractArray, b::AFAbstractArray, flag1::Cxx.CppEnum{:af_mat_prop}, flag2::Cxx.CppEnum{:af_mat_prop}) = icxx"af::matmul($a,$b,$flag1,$flag2);"
af_matmulNT(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::matmulNT($a,$b);"
af_matmulTN(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::matmulTN($a,$b);"
af_matmulTT(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::matmulTT($a,$b);"

af_transpose(a::AFAbstractArray) = icxx"af::transpose($a);"
af_ctranspose(a::AFAbstractArray) = icxx"af::transpose($a, true);"
af_transposeInPlace(a::AFAbstractArray) = icxx"af::transposeInPlace($a);"
af_ctransposeInPlace(a::AFAbstractArray) = icxx"af:transposeInPlace($a, true);"

af_solve(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::solve($a,$b);"
af_solveLU(a::AFAbstractArray, b::AFAbstractArray) = icxx"af_solveLU($a,$b);"

af_cholesky(out::vcpp"af::array", a::AFMatrix, flag::Bool) = icxx"af::cholesky($out,$a,$flag);"
af_choleskyInPlace(a::AFMatrix, flag::Bool) = icxx"af::choleskyInPlace($a,$flag);"
af_lu(l::vcpp"af::array", u::vcpp"af::array", p::vcpp"af::array", a::AFMatrix) = icxx"af::lu($l,$u,$p,$a);"
af_qr(q::vcpp"af::array", r::vcpp"af::array", tau::vcpp"af::array", a::AFMatrix) = icxx"af::qr($q,$r,$tau,$a);"
af_svd(u::vcpp"af::array", s::vcpp"af::array", vt::vcpp"af::array", a::AFMatrix) = icxx"af::svd($u,$s, $vt,$a);"
af_upper(a::Union{AFMatrix, vcpp"af::array"}) = icxx"af::upper($a);"
af_lower(a::Union{AFMatrix, vcpp"af::array"}) = icxx"af::lower($a);"

af_det(a::AFAbstractArray) = icxx"af::det<float>($a);"
af_inverse(a::AFAbstractArray) = icxx"af::inverse($a);"
af_norm(a::AFAbstractArray) = icxx"af::norm($a);"
af_rank(a::AFAbstractArray) = icxx"af::rank($a);"

#Signal 
af_fft(a::AFVector) = icxx"af::fft($a);"
af_fft2(a::AFMatrix) = icxx"af::fft2($a);"
af_fft3(a::AFAbstractArray) = icxx"af::fft3($a);"
af_ifft(a::AFVector) = icxx"af::ifft($a);"
af_ifft2(a::AFMatrix) = icxx"af::ifft2($a);"
af_ifft3(a::AFAbstractArray) = icxx"af::ifft3($a);"

af_fir(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::fir($a, $b);"
af_iir(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray) = icxx"af::iir($a, $b, $c);"

const AF_CONV_EXPAND = icxx"AF_CONV_EXPAND;"

af_convolve(a::AFAbstractArray, b::AFAbstractArray) = icxx"af_convolve($a,$b);"
af_convolve1(a::AFAbstractArray, b::AFAbstractArray) = icxx"af_convolve1($a,$b);"
af_convolve2(a::AFAbstractArray, b::AFAbstractArray) = icxx"af_convolve2($a,$b);"
af_convolve3(a::AFAbstractArray, b::AFAbstractArray) = icxx"af_convolve3($a,$b);"
af_fftConvolve(a::AFAbstractArray, b::AFAbstractArray, flag::Cxx.CppEnum{:af_conv_mode}) = icxx"af::fftConvolve($a,$b,$flag);"
af_fftConvolve2(a::AFAbstractArray, b::AFAbstractArray, flag::Cxx.CppEnum{:af_conv_mode}) = icxx"af::fftConvolve2($a,$b,$flag);"
af_fftConvolve3(a::AFAbstractArray, b::AFAbstractArray, flag::Cxx.CppEnum{:af_conv_mode}) = icxx"af::fftConvolve3($a,$b,$flag);"

af_approx1(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::approx1($a,$b);"
af_approx2(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray) = icxx"af::approx2($a,$b,$c);"

#Statistics
af_mean(a::AFAbstractArray) = icxx"af::mean<float>($a);"
af_mean(a::AFAbstractArray, b::Integer) = icxx"af::mean($a, $b);"
af_mean(a::AFAbstractArray, b::AFAbstractArray, c::Integer) = icxx"af::mean($a, $b, $c);"
af_median(a::AFAbstractArray) = icxx"af::median<float>($a);"
af_median(a::AFAbstractArray, b::Integer) = icxx"af::median($a,$b);"
af_stdev(a::AFAbstractArray) = icxx"af::stdev<float>($a);"
af_stdev(a::AFAbstractArray, b::Integer) = icxx"af::stdev($a,$b);"
af_var(a::AFAbstractArray) = icxx"af::var<float>($a);"
af_var(a::AFAbstractArray, b::Integer) = icxx"af::var($a,$b);"
af_var(a::AFAbstractArray, b::AFAbstractArray, c::Integer) = icxx"af::var($a,$b,$c);"
af_cov(a::AFAbstractArray, b::Integer) = icxx"af::cov($a,$b);"
af_corrcoef(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::corrcoef($a, $b);"

#Vector algorithms
af_sum(a::AFAbstractArray) = icxx"af::sum<float>($a);"
af_sum(a::AFAbstractArray, b::Integer) = icxx"af::sum($a,$b);"
af_product(a::AFAbstractArray) = icxx"af::product<float>($a);"
af_product(a::AFAbstractArray, b::Integer) = icxx"af::product($a,$b);"
af_max(a::AFAbstractArray) = icxx"af::max<float>($a);"
af_max(a::AFAbstractArray, b::Real) = icxx"af::max($a,$b);"
af_max_dim(a::AFAbstractArray, b::Integer) = icxx"""
    int dim = $b;
    af::max($a,dim);
    """
af_min(a::AFAbstractArray) = icxx"af::min<float>($a);"
af_min(a::AFAbstractArray, b::Real) = icxx"af::min($a,$b);"
af_min_dim(a::AFAbstractArray, b::Integer) = icxx"""
    int dim = $b;
    af::min($a,dim);
"""
af_anyTrue(a::AFAbstractArray) = icxx"af::anyTrue<bool>($a);"
af_anyTrue(a::AFAbstractArray, b::Integer) = icxx"af::anyTrue($a,$b);"
af_allTrue(a::AFAbstractArray) = icxx"af::allTrue<bool>($a);"
af_allTrue(a::AFAbstractArray, b::Integer) = icxx"af::allTrue($a,$b);"
af_count(a::AFAbstractArray) = icxx"af::count<int>($a);"
af_count(a::AFAbstractArray, b::Integer) = icxx"af::count($a,$b);"
af_accum(a::AFAbstractArray, b::Integer) = icxx"af::accum($a,$b);"
af_where(a::AFAbstractArray) = icxx"af::where($a);"
af_diff(a::AFAbstractArray, b::Integer) = icxx"af::diff($a,$b);"
af_diff2(a::AFAbstractArray, b::Integer) = icxx"af::diff2($a,$b);"
af_grad(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray) = icxx"af::grad($a,$b,$c);"
af_sort(a::AFAbstractArray, b::Integer, flag::Bool) = icxx"af::sort($a,$b,$flag);"
af_sort(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::Integer,flag::Bool) = icxx"af::sort($a,$b,$c,$d,$flag);"
af_sort(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::AFAbstractArray, el::Integer,flag::Bool) = icxx"af::sort($a,$b,$c,$d,$flag);"
af_setIntersect(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::setIntersect($a,$b);"
af_setUnion(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::setUnion($a,$b);"
af_setUnique(a::AFAbstractArray) = icxx"af::setUnique($a);"
af_flip(a::AFAbstractArray, b::Integer) = icxx"af::flip($a,$b);"
af_flat(a::AFAbstractArray) = icxx"af::flat($a);"

#Image
af_loadImage(a::AbstractString, b::Bool) = icxx"af::loadImage($a,$b);"
af_saveImage(a::AbstractString, b::AFAbstractArray) = icxx"af::saveImage($a,$b);"
af_rotate(a::AFAbstractArray, b::Real) = icxx"af::rotate($a,$b);"
af_scale(a::AFAbstractArray, b::Real, c::Real) = icxx"af::scale($a,$b,$c);"
af_skew(a::AFAbstractArray, b::Real, c::Real) = icxx"af::skew($a,$b,$c);"
af_translate(a::AFAbstractArray, b::Real, c::Real) = icxx"af::translate($a,$b,$c);"
af_transform(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::transform($a,$b);"

af_regions(a::AFAbstractArray) = icxx"$a.as(b8);"

const AF_RGB = icxx"AF_RGB;"
const AF_GRAY= icxx"AF_GRAY;"
const AF_HSV = icxx"AF_HSV;"

af_gray2rgb(a::AFAbstractArray) = icxx"af::gray2rgb($a);"
af_hsv2rgb(a::AFAbstractArray) = icxx"af::hsv2rgb($a);"
af_rgb2gray(a::AFAbstractArray) = icxx"af::rgb2gray($a);"
af_rgb2hsv(a::AFAbstractArray) = icxx"af::rgb2hsv($a);"
af_rgb2ycbcr(a::AFAbstractArray) = icxx"af::rgb2ycbcr($a);"
af_ycbcr2rgb(a::AFAbstractArray) = icxx"af::ycbcr2rgb($a);"
af_colorspace(a,b::Cxx.CppEnum{:af_cspace_t},c::Cxx.CppEnum{:af_cspace_t}) = icxx"af::colorspace($a,$b,$c);"

af_SAT(a::AFAbstractArray) = icxx"af::sat($a);"
af_bilateral(a::AFAbstractArray, b::Real, c::Real) = icxx"af::bilateral($a,$b,$c);"
af_maxfilt(a::AFAbstractArray) = icxx"af::maxfilt($a);"
af_medfilt(a::AFAbstractArray) = icxx"af::medfilt($a);"
af_minfilt(a::AFAbstractArray) = icxx"af::minfilt($a);"
af_sobel(a::AFAbstractArray) = icxx"af::sobel($a);"
af_meanShift(a::AFAbstractArray, b::Real, c::Real, d::Integer) = icxx"af::meanShift($a,$b,$c,$d);"

af_histogram(a::AFAbstractArray, b::Integer) = icxx"af::histogram($a,$b);"
af_histogram(a::AFAbstractArray, b::Integer, c::Real, d::Real) = icxx"af::histogram($a,$b,$c,$d);"
af_histEqual(a::AFAbstractArray, b::AFAbstractArray) = icxx"af::histEqual($a, $b);"

af_dilate(a::AFAbstractArray,b::AFAbstractArray) = icxx"af::dilate($a,$b);"
af_dilate3(a::AFAbstractArray,b::AFAbstractArray) = icxx"af::dilate3($a,$b);"
af_erode(a::AFAbstractArray,b::AFAbstractArray) = icxx"af::erode($a,$b);"
af_erode3(a::AFAbstractArray,b::AFAbstractArray) = icxx"af::erode3($a,$b);"


#Computer Vision
const AF_SAD = icxx"AF_SAD;"
const AF_ZSAD = icxx"AF_ZSAD;"
const AF_LSAD = icxx"AF_LSAD;"
const AF_SSD = icxx"AF_SSD;"
const AF_ZSSD = icxx"AF_ZSSD;"
const AF_LSSD = icxx"AF_LSSD;"
const AF_NCC = icxx"AF_NCC;"
const AF_ZNCC = icxx"AF_ZNCC;"
const AF_SHD = icxx"AF_SHD;"

af_getX(a::AFFeatures) = icxx"$a.getX();"
af_getY(a::AFFeatures) = icxx"$a.getY();"
af_getSize(a::AFFeatures) = icxx"$a.getSize();"
af_getScore(a::AFFeatures) = icxx"$a.getScore();"
af_getOrientation(a::AFFeatures) = icxx"$a.getOrientation();"

af_gaussiankernel(a::Integer, b::Integer) = icxx"af::gaussianKernel($a,$b);"
af_dog(a::AFAbstractArray, b::Integer, c::Integer) = icxx"af::dog($a,$b,$c);"
af_matchTemplate(a::AFAbstractArray, b::AFAbstractArray, c::Cxx.CppEnum{:af_match_type}) = icxx"af::matchTemplate($a,$b,$c);"
af_orb(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::Real, 
        el::Integer, f::Real, g::Integer, h::Bool) = icxx"af::orb($a,$b,$c,$d,$el,$f,$g,$h);"
af_sift(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, d::Integer, 
        el::Real, f::Real, g::Real, h::Bool, i::Real,j::Real) = icxx"af::sift($a,$b,$c,$d,$el,$f,$g,$h,$i,$j);"
af_fast(a::AFAbstractArray, b::Real, c::Integer, d::Bool, el::Real, f::Integer) = icxx"af::fast($a,$b,$c,$d,$el,$f);"
af_harris(a::AFAbstractArray, b::Integer, c::Real, d::Real, el::Integer, f::Real) = icxx"af::harris($a,$b,$c,$d,$el,$f);"
af_susan(a::AFAbstractArray, b::Integer, c::Real, d::Real, el::Real, f::Integer) = icxx"af::susan($a,$b,$c,$d,$el,$f);"
af_hammingMatcher(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, 
                    d::AFAbstractArray, el::Integer, f::Integer) = icxx"af::hammingMatcher($a,$b,$c,$d,$el,$f);"
af_nearestNeighbour(a::AFAbstractArray, b::AFAbstractArray, c::AFAbstractArray, 
                        d::AFAbstractArray, el::Integer, f::Integer, 
                        g::Cxx.CppEnum{:af_match_type}) = icxx"af::nearestNeighbour($a,$b,$c,$d,$el,$f,$g);"

#Graphics
af_window() = icxx"af::Window();"
af_window(a) = icxx"af::Window($a);"
af_window(a,b,c) = icxx"af::Window($a,$b,$c);"
af_setTitle(a,b) = icxx"$a.setTitle($b);"
af_setImage(a,b,c) = icxx"$a.image($b,$c);"
af_setPlot(a,b,c,d) = icxx"$a.plot($b,$c,$d);"
af_show(a) = icxx"$a.show();"
af_setSurface(a,b) = icxx"$a.surface($b);"
