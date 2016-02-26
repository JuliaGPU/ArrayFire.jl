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
af_complex_conj(a) = icxx"$a.H();"
af_as(a,to) = icxx"$a.as($to);"
af_copy(a) = icxx"$a.copy();"
af_bytes(a) = icxx"$a.bytes();"
af_isempty(a) = icxx"$a.isempty();"
af_isscalar(a) = icxx"$a.isscalar();"
af_isvector(a) = icxx"$a.isvector();"
af_isrow(a) = icxx"$a.isrow();"
af_iscolumn(a) = icxx"$a.iscolumn();"
af_iscomplex(a) = icxx"$a.iscomplex();"
af_isreal(a) = icxx"$a.isreal();"
af_isdouble(a) = icxx"$a.isdouble();"
af_issingle(a) = icxx"$a.issingle();"
af_isrealfloating(a) = icxx"$a.isrealfloating();"
af_isfloating(a) = icxx"$a.isfloating();"
af_isinteger(a) = icxx"$a.isinteger();"
af_isbool(a) = icxx"$a.isbool();"

#Helper functions
af_isNaN(a) = icxx"af::isNaN($a);"
af_isinf(a) = icxx"af::isInf($a);"
af_iszero(a) = icxx"af::iszero($a);"

#ArrayFire info functions
af_info() = icxx"af::info();"
af_isDoubleAvailable() = icxx"af::isDoubleAvailable();"
af_sync() = icxx"af::sync();"
af_sync(a) = icxx"af::sync($a);"
af_getDevice() = icxx"af::getDevice();"
af_setDevice(a) = icxx"af::setDevice($a);"

#Loading and Saving Arrays
af_saveArray(a,b,c) = icxx"af::saveArray($a,$b,$c);"
af_readArray(a,b) = icxx"af::readArray($a,$b);"

#Numeric functions
af_abs(a) = icxx"af::abs($a);"
af_max(a, b) = icxx"af::max($a, $b);"
af_min(a, b) = icxx"af::min($a, $b);"
af_arg(a) = icxx"af::arg($a);"
af_ceil(a) = icxx"af::cell($a);"
af_floor(a) = icxx"af::floor($a);"
af_hypot(a,b) = icxx"af::hypot($a, $b);"
af_mod(a,b) = icxx"af::mod($a, $b);"
of_rem(a,b) = icxx"af::rem($a, $b);"
af_round(a) = icxx"af::round($a);"
af_sign(a) = icxx"af::sign($a);"
af_trunc(a) = icxx"af_trunc($a);"

#Logical ops
af_equals(a, b) = AFArray{Bool}(icxx"$a == $b;")
af_gt(a,b) = AFArray{Bool}(icxx"$a > $b;")
af_ge(a,b) = AFArray{Bool}(icxx"$a >= $b;")
af_bitand(a,b) = AFArray{Bool}(icxx"$a & $b;")
af_and(a,b) = AFArray{Bool}(icxx"$a && $b;")
af_bitor(a,b) = AFArray{Bool}(icxx"$a | $b;")
af_bitxor(a,b) = AFArray{Bool}(icxx"$a ^ $b;")
af_le(a,b) = AFArray{Bool}(icxx"$a <= $b;")
af_lt(a,b) = AFArray{Bool}(icxx"$a < $b;")
af_neg(a) = icxx"-$a;"
af_neq(a,b) = AFArray{Bool}(icxx"$a != $b;")
af_not(a) = AFArray{Bool}(icxx"!$a;")
af_or(a,b) = AFArray{Bool}(icxx"$a || $b;")

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

af_cholesky(out::AFAbstractArray, a::AFAbstractArray, flag::AFAbstractArray) = icxx"af::cholesky($out,$a,$flag);"
af_choleskyInPlace(a::AFAbstractArray, flag) = icxx"af::choleskyInPlace($a,$flag);"
af_lu(l::AFAbstractArray, u::AFAbstractArray, p::AFAbstractArray, a::AFAbstractArray) = icxx"af::lu($l,$u,$p,$a);"
af_qr(q::AFAbstractArray, r::AFAbstractArray, tau::AFAbstractArray, a::AFAbstractArray) = icxx"af::qr($q,$r,$tau,$a);"
af_svd(u::AFAbstractArray, s::AFAbstractArray, vt::AFAbstractArray, a::AFAbstractArray) = icxx"af::svd($u,$s, $vt,$a);"
af_upper(a::AFAbstractArray) = icxx"af::upper($a);"
af_lower(a::AFAbstractArray) = icxx"af::lower($a);"

af_det(a::AFAbstractArray) = icxx"af::det<float>($a);"
af_inverse(a::AFAbstractArray) = icxx"af::inverse($a);"
af_norm(a::AFAbstractArray) = icxx"af::norm($a);"
af_rank(a::AFAbstractArray) = icxx"af::rank($a);"

#Signal 
af_fft(a) = icxx"af::fft($a);"
af_fft2(a) = icxx"af::fft2($a);"
af_fft3(a) = icxx"af::fft3($a);"
af_ifft(a) = icxx"af::ifft($a);"
af_ifft2(a) = icxx"af::ifft2($a);"
af_ifft3(a) = icxx"af::ifft3($a);"

af_fir(a,b) = icxx"af::fir($a, $b);"
af_iir(a,b,c) = icxx"af::iir($a, $b, $c);"

const AF_CONV_EXPAND = icxx"AF_CONV_EXPAND;"

af_convolve(a,b) = icxx"af_convolve($a,$b);"
af_convolve1(a,b) = icxx"af_convolve1($a,$b);"
af_convolve2(a,b) = icxx"af_convolve2($a,$b);"
af_convolve3(a,b) = icxx"af_convolve3($a,$b);"
af_fftConvolve(a,b,flag) = icxx"af::fftConvolve($a,$b,$flag);"
af_fftConvolve2(a,b,flag) = icxx"af::fftConvolve2($a,$b,$flag);"
af_fftConvolve3(a,b,flag) = icxx"af::fftConvolve3($a,$b,$flag);"

af_approx1(a,b) = icxx"af::approx1($a,$b);"
af_approx2(a,b,c) = icxx"af::approx2($a,$b,$c);"

#Statistics
af_mean(a) = icxx"af::mean<float>($a);"
af_mean(a,b) = icxx"af::mean($a, $b);"
af_mean(a,b,c) = icxx"af::mean($a, $b, $c);"
af_median(a) = icxx"af::median<float>($a);"
af_median(a,b) = icxx"af::median($a,$b);"
af_stdev(a) = icxx"af::stdev($a);"
af_stdev(a,b) = icxx"af::stdev($a,$b);"
af_var(a) = icxx"af::var($a);"
af_var(a,b) = icxx"af::var($a,$b);"
af_var(a,b,c) = icxx"af::var($a,$b,$c);"
af_cov(a,b) = icxx"af::cov($a,$b);"
af_corrcoef(a,b) = icxx"af::corrcoef($a, $b);"

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
af_anyTrue(a::AFAbtractArray) = icxx"af::anyTrue<bool>($a);"
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

af_gaussiankernel(a,b) = icxx"af::gaussianKernel($a,$b);"

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

af_getX(a) = icxx"$a.getX();"
af_getY(a) = icxx"$a.getY();"
af_getSize(a) = icxx"$a.getSize();"
af_getScore(a) = icxx"$a.getScore();"
af_getOrientation(a) = icxx"$a.getOrientation();"

af_dog(a,b,c) = icxx"af::dog($a,$b,$c);"
af_matchTemplate(a,b,c) = icxx"af::matchTemplate($a,$b,$c);"
af_orb(a,b,c,d,el,f) = icxx"af::orb($a,$b,$c,$d,$el,$f);"
af_sift(a,b,c,d,el,f,g,h,i,j) = icxx"af::sift($a,$b,$c,$d,$el,$f,$g,$h,$i,$j);"
af_fast(a,b,c,d,el,f) = icxx"af::fast($a,$b,$c,$d,$el,$f);"
af_harris(a,b,c,d,el,f) = icxx"af::harris($a,$b,$c,$d,$el,$f);"
af_susan(a,b,c,d,el,f) = icxx"af::susan($a,$b,$c,$d,$el,$f);"
af_hammingMatcher(a,b,c,d,el,f) = icxx"af::hammingMatcher($a,$b,$c,$d,$el,$f);"
af_nearestNeighbour(a,b,c,d,el,f,g) = icxx"af::nearestNeighbour($a,$b,$c,$d,$el,$f,$g);"

#Graphics
af_window() = icxx"af::Window();"
af_window(a) = icxx"af::Window($a);"
af_window(a,b,c) = icxx"af::Window($a,$b,$c);"
af_setTitle(a,b) = icxx"$a.setTitle($b);"
af_setImage(a,b,c) = icxx"$a.image($b,$c);"
af_setPlot(a,b,c,d) = icxx"$a.plot($b,$c,$d);"
af_show(a) = icxx"$a.show();"
af_setSurface(a,b) = icxx"$a.surface($b);"
