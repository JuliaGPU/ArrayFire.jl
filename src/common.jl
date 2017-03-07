# Automatically generated using Clang.jl wrap_c, version 0.0.0

using Compat

export AF_BACKEND_CPU, AF_BACKEND_CUDA, AF_BACKEND_DEFAULT, AF_BACKEND_OPENCL, AF_BINARY_ADD, AF_BINARY_MAX
export AF_BINARY_MIN, AF_BINARY_MUL, AF_COLORMAP_BLUE, AF_COLORMAP_COLORS, AF_COLORMAP_DEFAULT, AF_COLORMAP_HEAT
export AF_COLORMAP_MOOD, AF_COLORMAP_RED, AF_COLORMAP_SPECTRUM, AF_CONNECTIVITY_4, AF_CONNECTIVITY_8, AF_CONV_AUTO
export AF_CONV_DEFAULT, AF_CONV_EXPAND, AF_CONV_FREQ, AF_CONV_SPATIAL, AF_ERR_ARG, AF_ERR_ARR_BKND_MISMATCH
export AF_ERR_BATCH, AF_ERR_DEVICE, AF_ERR_DIFF_TYPE, AF_ERR_DRIVER, AF_ERR_INTERNAL, AF_ERR_INVALID_ARRAY
export AF_ERR_LOAD_LIB, AF_ERR_LOAD_SYM, AF_ERR_NONFREE, AF_ERR_NOT_CONFIGURED, AF_ERR_NOT_SUPPORTED, AF_ERR_NO_DBL
export AF_ERR_NO_GFX, AF_ERR_NO_MEM, AF_ERR_RUNTIME, AF_ERR_SIZE, AF_ERR_TYPE, AF_ERR_UNKNOWN, AF_FIF_BMP
export AF_FIF_EXR, AF_FIF_HDR, AF_FIF_ICO, AF_FIF_JNG, AF_FIF_JP2, AF_FIF_JPEG, AF_FIF_PNG, AF_FIF_PPM
export AF_FIF_PPMRAW, AF_FIF_PSD, AF_FIF_RAW, AF_FIF_TIFF, AF_GRAY, AF_HOMOGRAPHY_LMEDS, AF_HOMOGRAPHY_RANSAC
export AF_HSV, AF_ID, AF_INTERP_BICUBIC, AF_INTERP_BICUBIC_SPLINE, AF_INTERP_BILINEAR, AF_INTERP_BILINEAR_COSINE
export AF_INTERP_CUBIC, AF_INTERP_CUBIC_SPLINE, AF_INTERP_LINEAR, AF_INTERP_LINEAR_COSINE, AF_INTERP_LOWER
export AF_INTERP_NEAREST, AF_LSAD, AF_LSSD, AF_MARKER_CIRCLE, AF_MARKER_CROSS, AF_MARKER_NONE, AF_MARKER_PLUS
export AF_MARKER_POINT, AF_MARKER_SQUARE, AF_MARKER_STAR, AF_MARKER_TRIANGLE, AF_MAT_BLOCK_DIAG, AF_MAT_CONJ
export AF_MAT_CTRANS, AF_MAT_DIAG_UNIT, AF_MAT_LOWER, AF_MAT_NONE, AF_MAT_ORTHOG, AF_MAT_POSDEF, AF_MAT_SYM
export AF_MAT_TRANS, AF_MAT_TRI_DIAG, AF_MAT_UPPER, AF_MOMENT_FIRST_ORDER, AF_MOMENT_M00, AF_MOMENT_M01
export AF_MOMENT_M10, AF_MOMENT_M11, AF_NCC, AF_NORM_EUCLID, AF_NORM_MATRIX_1, AF_NORM_MATRIX_2, AF_NORM_MATRIX_INF
export AF_NORM_MATRIX_L_PQ, AF_NORM_VECTOR_1, AF_NORM_VECTOR_2, AF_NORM_VECTOR_INF, AF_NORM_VECTOR_P, AF_PAD_SYM
export AF_PAD_ZERO, AF_RANDOM_ENGINE_DEFAULT, AF_RANDOM_ENGINE_MERSENNE, AF_RANDOM_ENGINE_MERSENNE_GP11213
export AF_RANDOM_ENGINE_PHILOX, AF_RANDOM_ENGINE_PHILOX_4X32_10, AF_RANDOM_ENGINE_THREEFRY, AF_RANDOM_ENGINE_THREEFRY_2X32_16
export AF_RGB, AF_SAD, AF_SHD, AF_SSD, AF_STORAGE_COO, AF_STORAGE_CSC, AF_STORAGE_CSR, AF_STORAGE_DENSE
export AF_SUCCESS, AF_YCC_2020, AF_YCC_601, AF_YCC_709, AF_YCbCr, AF_ZNCC, AF_ZSAD, AF_ZSSD, afDevice
export afHost, af_array, af_backend, af_binary_op, af_border_type, af_colormap, af_connectivity, af_conv_domain
export af_conv_mode, af_cspace_t, af_dtype, af_err, af_features, af_homography_type, af_image_format, af_interp_type
export af_marker_type, af_mat_prop, af_match_type, af_moment_type, af_norm_type, af_random_engine, af_random_engine_type
export af_someenum_t, af_source, af_storage, af_window, af_ycc_std, b8, c32, c64, dim_t, f32, f64, intl
export s16, s32, s64, u16, u32, u64, u8, uintl

typealias dim_t Clonglong
typealias intl Clonglong
typealias uintl Culonglong

# begin enum af_err
typealias af_err UInt32
const AF_SUCCESS = (UInt32)(0)
const AF_ERR_NO_MEM = (UInt32)(101)
const AF_ERR_DRIVER = (UInt32)(102)
const AF_ERR_RUNTIME = (UInt32)(103)
const AF_ERR_INVALID_ARRAY = (UInt32)(201)
const AF_ERR_ARG = (UInt32)(202)
const AF_ERR_SIZE = (UInt32)(203)
const AF_ERR_TYPE = (UInt32)(204)
const AF_ERR_DIFF_TYPE = (UInt32)(205)
const AF_ERR_BATCH = (UInt32)(207)
const AF_ERR_DEVICE = (UInt32)(208)
const AF_ERR_NOT_SUPPORTED = (UInt32)(301)
const AF_ERR_NOT_CONFIGURED = (UInt32)(302)
const AF_ERR_NONFREE = (UInt32)(303)
const AF_ERR_NO_DBL = (UInt32)(401)
const AF_ERR_NO_GFX = (UInt32)(402)
const AF_ERR_LOAD_LIB = (UInt32)(501)
const AF_ERR_LOAD_SYM = (UInt32)(502)
const AF_ERR_ARR_BKND_MISMATCH = (UInt32)(503)
const AF_ERR_INTERNAL = (UInt32)(998)
const AF_ERR_UNKNOWN = (UInt32)(999)
# end enum af_err

# begin enum af_dtype
typealias af_dtype UInt32
const f32 = (UInt32)(0)
const c32 = (UInt32)(1)
const f64 = (UInt32)(2)
const c64 = (UInt32)(3)
const b8 = (UInt32)(4)
const s32 = (UInt32)(5)
const u32 = (UInt32)(6)
const u8 = (UInt32)(7)
const s64 = (UInt32)(8)
const u64 = (UInt32)(9)
const s16 = (UInt32)(10)
const u16 = (UInt32)(11)
# end enum af_dtype

# begin enum af_source
typealias af_source UInt32
const afDevice = (UInt32)(0)
const afHost = (UInt32)(1)
# end enum af_source

typealias af_array Ptr{Void}

# begin enum af_interp_type
typealias af_interp_type UInt32
const AF_INTERP_NEAREST = (UInt32)(0)
const AF_INTERP_LINEAR = (UInt32)(1)
const AF_INTERP_BILINEAR = (UInt32)(2)
const AF_INTERP_CUBIC = (UInt32)(3)
const AF_INTERP_LOWER = (UInt32)(4)
const AF_INTERP_LINEAR_COSINE = (UInt32)(5)
const AF_INTERP_BILINEAR_COSINE = (UInt32)(6)
const AF_INTERP_BICUBIC = (UInt32)(7)
const AF_INTERP_CUBIC_SPLINE = (UInt32)(8)
const AF_INTERP_BICUBIC_SPLINE = (UInt32)(9)
# end enum af_interp_type

# begin enum af_border_type
typealias af_border_type UInt32
const AF_PAD_ZERO = (UInt32)(0)
const AF_PAD_SYM = (UInt32)(1)
# end enum af_border_type

# begin enum af_connectivity
typealias af_connectivity UInt32
const AF_CONNECTIVITY_4 = (UInt32)(4)
const AF_CONNECTIVITY_8 = (UInt32)(8)
# end enum af_connectivity

# begin enum af_conv_mode
typealias af_conv_mode UInt32
const AF_CONV_DEFAULT = (UInt32)(0)
const AF_CONV_EXPAND = (UInt32)(1)
# end enum af_conv_mode

# begin enum af_conv_domain
typealias af_conv_domain UInt32
const AF_CONV_AUTO = (UInt32)(0)
const AF_CONV_SPATIAL = (UInt32)(1)
const AF_CONV_FREQ = (UInt32)(2)
# end enum af_conv_domain

# begin enum af_match_type
typealias af_match_type UInt32
const AF_SAD = (UInt32)(0)
const AF_ZSAD = (UInt32)(1)
const AF_LSAD = (UInt32)(2)
const AF_SSD = (UInt32)(3)
const AF_ZSSD = (UInt32)(4)
const AF_LSSD = (UInt32)(5)
const AF_NCC = (UInt32)(6)
const AF_ZNCC = (UInt32)(7)
const AF_SHD = (UInt32)(8)
# end enum af_match_type

# begin enum af_ycc_std
typealias af_ycc_std UInt32
const AF_YCC_601 = (UInt32)(601)
const AF_YCC_709 = (UInt32)(709)
const AF_YCC_2020 = (UInt32)(2020)
# end enum af_ycc_std

# begin enum af_cspace_t
typealias af_cspace_t UInt32
const AF_GRAY = (UInt32)(0)
const AF_RGB = (UInt32)(1)
const AF_HSV = (UInt32)(2)
const AF_YCbCr = (UInt32)(3)
# end enum af_cspace_t

# begin enum af_mat_prop
typealias af_mat_prop UInt32
const AF_MAT_NONE = (UInt32)(0)
const AF_MAT_TRANS = (UInt32)(1)
const AF_MAT_CTRANS = (UInt32)(2)
const AF_MAT_CONJ = (UInt32)(4)
const AF_MAT_UPPER = (UInt32)(32)
const AF_MAT_LOWER = (UInt32)(64)
const AF_MAT_DIAG_UNIT = (UInt32)(128)
const AF_MAT_SYM = (UInt32)(512)
const AF_MAT_POSDEF = (UInt32)(1024)
const AF_MAT_ORTHOG = (UInt32)(2048)
const AF_MAT_TRI_DIAG = (UInt32)(4096)
const AF_MAT_BLOCK_DIAG = (UInt32)(8192)
# end enum af_mat_prop

# begin enum af_norm_type
typealias af_norm_type UInt32
const AF_NORM_VECTOR_1 = (UInt32)(0)
const AF_NORM_VECTOR_INF = (UInt32)(1)
const AF_NORM_VECTOR_2 = (UInt32)(2)
const AF_NORM_VECTOR_P = (UInt32)(3)
const AF_NORM_MATRIX_1 = (UInt32)(4)
const AF_NORM_MATRIX_INF = (UInt32)(5)
const AF_NORM_MATRIX_2 = (UInt32)(6)
const AF_NORM_MATRIX_L_PQ = (UInt32)(7)
const AF_NORM_EUCLID = (UInt32)(2)
# end enum af_norm_type

# begin enum af_image_format
typealias af_image_format UInt32
const AF_FIF_BMP = (UInt32)(0)
const AF_FIF_ICO = (UInt32)(1)
const AF_FIF_JPEG = (UInt32)(2)
const AF_FIF_JNG = (UInt32)(3)
const AF_FIF_PNG = (UInt32)(13)
const AF_FIF_PPM = (UInt32)(14)
const AF_FIF_PPMRAW = (UInt32)(15)
const AF_FIF_TIFF = (UInt32)(18)
const AF_FIF_PSD = (UInt32)(20)
const AF_FIF_HDR = (UInt32)(26)
const AF_FIF_EXR = (UInt32)(29)
const AF_FIF_JP2 = (UInt32)(31)
const AF_FIF_RAW = (UInt32)(34)
# end enum af_image_format

# begin enum af_moment_type
typealias af_moment_type UInt32
const AF_MOMENT_M00 = (UInt32)(1)
const AF_MOMENT_M01 = (UInt32)(2)
const AF_MOMENT_M10 = (UInt32)(4)
const AF_MOMENT_M11 = (UInt32)(8)
const AF_MOMENT_FIRST_ORDER = (UInt32)(15)
# end enum af_moment_type

# begin enum af_homography_type
typealias af_homography_type UInt32
const AF_HOMOGRAPHY_RANSAC = (UInt32)(0)
const AF_HOMOGRAPHY_LMEDS = (UInt32)(1)
# end enum af_homography_type

# begin enum af_backend
typealias af_backend UInt32
const AF_BACKEND_DEFAULT = (UInt32)(0)
const AF_BACKEND_CPU = (UInt32)(1)
const AF_BACKEND_CUDA = (UInt32)(2)
const AF_BACKEND_OPENCL = (UInt32)(4)
# end enum af_backend

# begin enum af_someenum_t
typealias af_someenum_t UInt32
const AF_ID = (UInt32)(0)
# end enum af_someenum_t

# begin enum af_binary_op
typealias af_binary_op UInt32
const AF_BINARY_ADD = (UInt32)(0)
const AF_BINARY_MUL = (UInt32)(1)
const AF_BINARY_MIN = (UInt32)(2)
const AF_BINARY_MAX = (UInt32)(3)
# end enum af_binary_op

# begin enum af_random_engine_type
typealias af_random_engine_type UInt32
const AF_RANDOM_ENGINE_PHILOX_4X32_10 = (UInt32)(100)
const AF_RANDOM_ENGINE_THREEFRY_2X32_16 = (UInt32)(200)
const AF_RANDOM_ENGINE_MERSENNE_GP11213 = (UInt32)(300)
const AF_RANDOM_ENGINE_PHILOX = (UInt32)(100)
const AF_RANDOM_ENGINE_THREEFRY = (UInt32)(200)
const AF_RANDOM_ENGINE_MERSENNE = (UInt32)(300)
const AF_RANDOM_ENGINE_DEFAULT = (UInt32)(100)
# end enum af_random_engine_type

# begin enum af_colormap
typealias af_colormap UInt32
const AF_COLORMAP_DEFAULT = (UInt32)(0)
const AF_COLORMAP_SPECTRUM = (UInt32)(1)
const AF_COLORMAP_COLORS = (UInt32)(2)
const AF_COLORMAP_RED = (UInt32)(3)
const AF_COLORMAP_MOOD = (UInt32)(4)
const AF_COLORMAP_HEAT = (UInt32)(5)
const AF_COLORMAP_BLUE = (UInt32)(6)
# end enum af_colormap

# begin enum af_marker_type
typealias af_marker_type UInt32
const AF_MARKER_NONE = (UInt32)(0)
const AF_MARKER_POINT = (UInt32)(1)
const AF_MARKER_CIRCLE = (UInt32)(2)
const AF_MARKER_SQUARE = (UInt32)(3)
const AF_MARKER_TRIANGLE = (UInt32)(4)
const AF_MARKER_CROSS = (UInt32)(5)
const AF_MARKER_PLUS = (UInt32)(6)
const AF_MARKER_STAR = (UInt32)(7)
# end enum af_marker_type

# begin enum af_storage
typealias af_storage UInt32
const AF_STORAGE_DENSE = (UInt32)(0)
const AF_STORAGE_CSR = (UInt32)(1)
const AF_STORAGE_CSC = (UInt32)(2)
const AF_STORAGE_COO = (UInt32)(3)
# end enum af_storage

type af_seq
    _begin::Cdouble
    _end::Cdouble
    step::Cdouble
end

type af_index_t
    idx::Void
    isSeq::Bool
    isBatch::Bool
end

type af_cfloat
    real::Cfloat
    imag::Cfloat
end

type af_cdouble
    real::Cdouble
    imag::Cdouble
end

typealias af_features Ptr{Void}
typealias af_window Culonglong

type af_cell
    row::Cint
    col::Cint
    title::Cstring
    cmap::af_colormap
end

typealias af_random_engine Ptr{Void}
