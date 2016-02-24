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

af_dot(a,b) = icxx"af::dot($a, $b);"

af_matmul(a,b) = icxx"af::matmul($a, $b);"
af_matmul3(a,b,c) = icxx"af::matmul($a, $b, $c);"
af_matmul4(a,b,c,d) = icxx"af::matmul($a, $b, $c, $d);"
af_matmul_flags(a,b,flag1,flag2) = icxx"af::matmul($a,$b,$flag1,$flag2);"
af_matmulNT(a,b) = icxx"af::matmulNT($a,$b);"
af_matmulTN(a,b) = icxx"af::matmulTN($a,$b);"
af_matmulTT(a,b) = icxx"af::matmulTT($a,$b);"

af_transpose(a) = icxx"af::transpose($a);"
af_ctranspose(a) = icxx"af::transpose($a, true);"
af_transposeInPlace(a) = icxx"af::transposeInPlace($a);"
af_ctransposeInPlace(a) = icxx"af:transposeInPlace($a, true);"

af_solve(a,b) = icxx"af::solve($a,$b);"
af_solveLU(a,b) = icxx"af_solveLU($a,$b);"

af_cholesky(out,a,flag) = icxx"af::cholesky($out,$a,$flag);"
af_choleskyInPlace(a, flag) = icxx"af::choleskyInPlace($a,$flag);"
af_lu(l,u,p,a) = icxx"af::lu($l,$u,$p,$a);"
af_qr(q,r,tau,a) = icxx"af::qr($q,$r,$tau,$a);"
af_svd(u,s,vt,a) = icxx"af::svd($u,$s, $vt,$a);"

af_det(a) = icxx"af::det<float>($a);"
af_inverse(a) = icxx"af::inverse($a);"
af_norm(a) = icxx"af::norm($a);"
af_rank(a) = icxx"af::rank($a);"

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
af_sum(a) = icxx"af::sum<float>($a);"
af_sum(a,b) = icxx"af::sum($a,$b);"
af_product(a) = icxx"af::product<float>($a);"
af_product(a,b) = icxx"af::product($a,$b);"
af_max(a) = icxx"af::max<float>($a);"
af_max(a,b) = icxx"af::max($a,$b);"
af_max_dim(a,b) = icxx"""
    int dim = $b;
    af::max($a,dim);
    """
af_min(a) = icxx"af::min<float>($a);"
af_min(a,b) = icxx"af::min($a,$b);"
af_min_dim(a,b) = icxx"""
    int dim = $b;
    af::min($a,dim);
"""
af_anyTrue(a) = icxx"af::anyTrue<bool>($a);"
af_anyTrue(a,b) = icxx"af::anyTrue($a,$b);"
af_allTrue(a) = icxx"af::allTrue<bool>($a);"
af_allTrue(a,b) = icxx"af::allTrue($a,$b);"
af_count(a) = icxx"af::count<int>($a);"
af_count(a,b) = icxx"af::count($a,$b);"
af_accum(a,b) = icxx"af::accum($a,$b);"
af_where(a) = icxx"af::where($a);"
af_diff(a,b) = icxx"af::diff($a,$b);"
af_diff2(a,b) = icxx"af::diff2($a,$b);"
af_grad(a,b,c) = icxx"af::grad($a,$b,$c);"
af_sort(a,b,flag) = icxx"af::sort($a,$b,$flag);"
af_sort(a,b,c,d,flag) = icxx"af::sort($a,$b,$c,$d,$flag);"
af_sort(a,b,c,d,el,flag) = icxx"af::sort($a,$b,$c,$d,$flag);"
af_setIntersect(a,b) = icxx"af::setIntersect($a,$b);"
af_setUnion(a,b) = icxx"af::setUnion($a,$b);"
af_setUnique(a) = icxx"af::setUnique($a);"
af_flip(a,b) = icxx"af::flip($a,$b);"
af_flat(a) = icxx"af::flat($a);"

#Image
af_loadImage(a,b) = icxx"af::loadImage($a,$b);"
af_saveImage(a,b) = icxx"af::saveImage($a,$b);"
af_rotate(a,b) = icxx"af::rotate($a,$b);"
af_scale(a,b,c) = icxx"af::scale($a,$b,$c);"
af_skew(a,b,c) = icxx"af::skew($a,$b,$c);"
af_translate(a,b,c) = icxx"af::translate($a,$b,$c);"
af_transform(a,b) = icxx"af::transform($a,$b);"

af_regions(a) = icxx"$a.as(b8);"

const AF_RGB = icxx"AF_RGB;"
const AF_GRAY= icxx"AF_GRAY;"
const AF_HSV = icxx"AF_HSV;"

af_gray2rgb(a) = icxx"af::gray2rgb($a);"
af_hsv2rgb(a) = icxx"af::hsv2rgb($a);"
af_rgb2gray(a) = icxx"af::rgb2gray($a);"
af_rgb2hsv(a) = icxx"af::rgb2hsv($a);"
af_rgb2ycbcr(a) = icxx"af::rgb2ycbcr($a);"
af_ycbcr2rgb(a) = icxx"af::ycbcr2rgb($a);"
af_colorspace(a,b,c) = icxx"af::colorspace($a,$b,$c);"

af_SAT(a) = icxx"af::sat($a);"
af_bilateral(a,b,c) = icxx"af::bilateral($a,$b,$c);"
af_maxfilt(a) = icxx"af::maxfilt($a);"
af_medfilt(a) = icxx"af::medfilt($a);"
af_minfilt(a) = icxx"af::minfilt($a);"
af_sobel(a) = icxx"af::sobel($a);"
af_meanShift(a,b,c,d) = icxx"af::meanShift($a,$b,$c,$d);"

af_histogram(a,b) = icxx"af::histogram($a,$b);"
af_histogram(a,b,c,d) = icxx"af::histogram($a,$b,$c,$d);"
af_histEqual(a,b) = icxx"af::histEqual($a, $b);"

af_dilate(a,b) = icxx"af::dilate($a,$b);"
af_dilate3(a,b) = icxx"af::dilate3($a,$b);"
af_erode(a,b) = icxx"af::erode($a,$b);"
af_erode3(a,b) = icxx"af::erode3($a,$b);"

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

af_dog(a,b,c) = icxx"af::dog($a,$b,$c);"
af_matchTemplate(a,b,c) = icxx"af::matchTemplate($a,$b,$c);"
