# Julia wrapper for header: /usr/include/arrayfire.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


export abs, accum, acos, acosh, add, afinfo, afinit, afversion, all, all_true_all, alloc_device, and, anisotropic_diffusion
export any, any_true_all, approx1, approx2, arg, array_to_string, asin, asinh, assign_seq, atan, atan2
export atanh, bilateral, bitand, bitor, bitshiftl, bitshiftr, bitxor, canny, cbrt, ceil, cholesky_inplace
export color_space, complex, conj, convolve1, convolve2, convolve2_sep, convolve3, copy, corrcoef, cos
export cosh, count, count_all, cov, create_features, create_handle, create_indexers, create_random_engine
export create_sparse_array, create_sparse_array_from_dense, create_sparse_array_from_ptr, create_window
export delete_image_memory, destroy_window, det, device_array, device_gc, diag, diagm
export diff1, diff2, dilate, dilate3, div, dog, dot, dot_all, draw_hist, draw_image, draw_plot, draw_plot3
export draw_plot_2d, draw_plot_3d, draw_plot_nd, draw_scatter, draw_scatter3, draw_scatter_2d, draw_scatter_3d
export draw_scatter_nd, draw_surface, draw_vector_field_2d, draw_vector_field_3d, draw_vector_field_nd
export eq, erf, erfc, erode, erode3, exp, expm1, factorial, fast, flip, floor, free_device, full, gaussian_kernel
export ge, get_active_backend, get_available_backends, get_backend_count, get_backend_id, get_data_ptr
export get_dbl_support, get_default_random_engine, get_device, get_device_count, get_device_id, get_device_ptr
export get_dims, get_elements, get_features_num, get_features_orientation, get_features_score, get_features_size
export get_features_xpos, get_features_ypos, get_manual_eval_flag, get_mem_step_size, get_revision, get_seed
export gloh, gradient, gray2rgb, gt, hamming_matcher, harris, hist_equal, histogram, homography, hsv2rgb
export hypot, identity, imag, imax, imax_all, imin, imin_all, index, index_gen, inverse, iota, is_bool
export is_column, is_complex, is_double, is_empty, is_floating, is_image_io_available, is_integer, is_lapack_available
export is_locked_array, is_real, is_realfloating, is_row, is_scalar, is_single, is_vector, is_window_closed
export isinf, isnan, issparse, iszero, le, lgamma, load_image, load_image_memory, load_image_native, lock_array
export lock_device_ptr, log, log10, log1p, log2, lookup, lower, lt, lu, lu_inplace, make_seq, match_template
export matmul, max_all, maxfilt, maximum, maxof, mean_all, mean_all_weighted, mean_shift, medfilt, medfilt1
export medfilt2, median_all, min_all, minfilt, minimum, minof, mod, moments, moments_all, mul, nearest_neighbour
export neq, norm, not, or, orb, pow, pow2, print_array, print_array_gen, print_mem_info, prod, product_all
export product_nan, product_nan_all, qr, qr_inplace, random_engine_get_seed, random_engine_get_type, random_engine_set_seed
export random_engine_set_type, random_normal, random_uniform, range, rank, read_array_index, read_array_key
export read_array_key_check, real, regions, release_features, release_random_engine, rem, reorder, replace
export replace, resize, retain_features, retain_random_engine, rgb2gray, rgb2hsv, rgb2ycbcr, root, rotate
export round, sat, save_array, save_image, save_image_memory, save_image_native, scale, scan, scan_by_key
export set_axes_limits_2d, set_axes_limits_3d, set_axes_limits_compute, set_axes_titles, set_backend, set_default_random_engine_type
export set_device, set_intersect, set_manual_eval_flag, set_mem_step_size, set_position, set_seed, set_size
export set_title, set_union, set_unique, set_visibility, shift, show_window, sift, sigmoid, signbit, sin, sinh
export skew, sobel_operator, solve, solve_lu, sort_by_key, sparse_convert_to, sparse_get_col_idx, sparse_get_info
export sparse_get_nnz, sparse_get_row_idx, sparse_get_storage, sparse_get_values, sqrt, stdev_all, sub
export sum, sum_all, sum_nan, sum_nan_all, susan, svd_inplace, sync, tan, tanh, tgamma, tile, topk, transform
export transform_coordinates, translate, transpose_inplace, trunc, unlock_array, unlock_device_ptr, unwrap
export upper, var_all, var_all_weighted, wrap, write_array, ycbcr2rgb

function sum(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray!(out[])
end

function sum_nan(_in::AFArray{T,N},dim::Integer,nanval::Real) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sum_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in.arr,Cint(dim - 1),Cdouble(nanval)))
    AFArray{T,N}(out[])
end

function prod(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_product,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function product_nan(_in::AFArray{T,N},dim::Integer,nanval::Real) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_product_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in.arr,Cint(dim - 1),Cdouble(nanval)))
    AFArray{T,N}(out[])
end

function minimum(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_min,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function maximum(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_max,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function all(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_all_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function any(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_any_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function count(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_count,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function sum_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_sum_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function sum_nan_all(_in::AFArray,nanval::Real)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_sum_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in.arr,Cdouble(nanval)))
    (real[],imag[])
end

function product_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_product_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function product_nan_all(_in::AFArray,nanval::Real)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_product_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in.arr,Cdouble(nanval)))
    (real[],imag[])
end

function min_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_min_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function max_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_max_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function all_true_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_all_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function any_true_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_any_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function count_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_count_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function imin(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    _error(ccall((:af_imin,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in.arr,Cint(dim - 1)))
    (AFArray!(out[]),AFArray!(idx[]))
end

function imax(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    _error(ccall((:af_imax,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in.arr,Cint(dim - 1)))
    (AFArray!(out[]),AFArray!(idx[]))
end

function imin_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    _error(ccall((:af_imin_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in.arr))
    (real[],imag[],idx[])
end

function imax_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    _error(ccall((:af_imax_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in.arr))
    (real[],imag[],idx[])
end

function accum(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_accum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function scan(_in::AFArray{T,N},dim::Integer,op::af_binary_op,inclusive_scan::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_scan,af_lib),af_err,(Ptr{af_array},af_array,Cint,af_binary_op,Bool),out,_in.arr,Cint(dim - 1),op,inclusive_scan))
    AFArray{T,N}(out[])
end

function scan_by_key(key::AFArray,_in::AFArray,dim::Integer,op::af_binary_op,inclusive_scan::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_scan_by_key,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cint,af_binary_op,Bool),out,key.arr,_in.arr,Cint(dim - 1),op,inclusive_scan))
    AFArray!(out[])
end

function diff1(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_diff1,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function diff2(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_diff2,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim - 1)))
    AFArray{T,N}(out[])
end

function sort_by_key(keys::AFArray,values::AFArray,dim::Integer,isAscending::Bool)
    out_keys = RefValue{af_array}(0)
    out_values = RefValue{af_array}(0)
    _error(ccall((:af_sort_by_key,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,UInt32,Bool),out_keys,out_values,keys.arr,values.arr,UInt32(dim - 1),isAscending))
    (AFArray!(out_keys[]),AFArray!(out_values[]))
end

function set_unique(_in::AFArray{T,N},is_sorted::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_unique,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_sorted))
    AFArray{T,N}(out[])
end

function set_union(first::AFArray,second::AFArray,is_unique::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_union,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first.arr,second.arr,is_unique))
    AFArray!(out[])
end

function set_intersect(first::AFArray,second::AFArray,is_unique::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_intersect,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first.arr,second.arr,is_unique))
    AFArray!(out[])
end

function add(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_add,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function sub(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sub,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function mul(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_mul,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function div(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_div,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function lt(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_lt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function gt(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_gt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function le(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_le,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function ge(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_ge,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function eq(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_eq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function neq(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_neq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function and(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_and,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function or(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_or,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{Bool,batched(N1,N2)}(out[])
end

function not(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_not,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Bool,N}(out[])
end

function bitand(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitand,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function bitor(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function bitxor(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitxor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function bitshiftl(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitshiftl,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function bitshiftr(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitshiftr,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function minof(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_minof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function maxof(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_maxof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function rem(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_rem,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function mod(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_mod,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function abs(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_abs,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function arg(_in::AFArray{Complex{T},N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_arg,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function signbit(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sign,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Float32,N}(out[])
end

function round(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_round,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function trunc(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_trunc,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function floor(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_floor,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function ceil(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_ceil,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function hypot(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_hypot,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function sin(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sin,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function cos(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_cos,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function tan(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_tan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function asin(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_asin,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function acos(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_acos,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function atan(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_atan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function atan2(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_atan2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function complex(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_cplx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Complex{T},N}(out[])
end

function real(_in::AFArray{Complex{T},N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_real,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function imag(_in::AFArray{Complex{T},N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_imag,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function conj(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_conjg,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function sinh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function cosh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_cosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function tanh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_tanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function asinh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_asinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function acosh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_acosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function atanh(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_atanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function root(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_root,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function pow(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},batch::Bool) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_pow,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function pow2(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_pow2,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function exp(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_exp,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function sigmoid(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sigmoid,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function expm1(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_expm1,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function erf(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_erf,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function erfc(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_erfc,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function log(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_log,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function log1p(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_log1p,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function log10(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_log10,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function log2(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_log2,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function sqrt(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sqrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function cbrt(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_cbrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function factorial(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_factorial,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function tgamma(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_tgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function lgamma(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_lgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function iszero(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_iszero,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Bool,N}(out[])
end

function isinf(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_isinf,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Bool,N}(out[])
end

function isnan(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_isnan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Bool,N}(out[])
end

function make_seq(_begin::Real,_end::Real,step::Real)
    ccall((:af_make_seq,af_lib),af_seq,(Cdouble,Cdouble,Cdouble),Cdouble(_begin),Cdouble(_end),Cdouble(step))
end

function print_array(arr::AFArray)
    _error(ccall((:af_print_array,af_lib),af_err,(af_array,),arr.arr))
end

function print_array_gen(exp,arr::AFArray,precision::Integer)
    _error(ccall((:af_print_array_gen,af_lib),af_err,(Cstring,af_array,Cint),exp,arr.arr,Cint(precision)))
end

function save_array(key,arr::AFArray,filename,append::Bool)
    index = RefValue{Cint}(0)
    _error(ccall((:af_save_array,af_lib),af_err,(Ptr{Cint},Cstring,af_array,Cstring,Bool),index,key,arr.arr,filename,append))
    index[]
end

function read_array_index(filename,index::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_read_array_index,af_lib),af_err,(Ptr{af_array},Cstring,UInt32),out,filename,UInt32(index)))
    AFArray!(out[])
end

function read_array_key(filename,key)
    out = RefValue{af_array}(0)
    _error(ccall((:af_read_array_key,af_lib),af_err,(Ptr{af_array},Cstring,Cstring),out,filename,key))
    AFArray!(out[])
end

function read_array_key_check(filename,key)
    index = RefValue{Cint}(0)
    _error(ccall((:af_read_array_key_check,af_lib),af_err,(Ptr{Cint},Cstring,Cstring),index,filename,key))
    index[]
end

function array_to_string(exp,arr::AFArray,precision::Integer,transpose::Bool)
    output = RefValue{Cstring}()
    _error(ccall((:af_array_to_string,af_lib),af_err,(Ptr{Cstring},Cstring,af_array,Cint,Bool),output,exp,arr.arr,Cint(precision),transpose))
    output[]
end

function afversion()
    major = RefValue{Cint}(0)
    minor = RefValue{Cint}(0)
    patch = RefValue{Cint}(0)
    _error(ccall((:af_get_version,af_lib),af_err,(Ptr{Cint},Ptr{Cint},Ptr{Cint}),major,minor,patch))
    (major[],minor[],patch[])
end

function get_revision()
    ccall((:af_get_revision,af_lib),Cstring,())
end

function index(_in::AFArray{T,N},ndims::Integer,index) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_index,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq}),out,_in.arr,UInt32(ndims),index))
    AFArray{T,N}(out[])
end

function lookup(_in::AFArray,indices::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_lookup,af_lib),af_err,(Ptr{af_array},af_array,af_array,UInt32),out,_in.arr,indices.arr,UInt32(dim - 1)))
    AFArray!(out[])
end

function assign_seq(lhs::AFArray,ndims::Integer,indices,rhs::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_assign_seq,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq},af_array),out,lhs.arr,UInt32(ndims),indices,rhs.arr))
    AFArray!(out[])
end

function index_gen(_in::AFArray{T,N},ndims::dim_t,indices) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_index_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t}),out,_in.arr,ndims,indices))
    AFArray{T,N}(out[])
end

function create_indexers()
    indexers = RefValue{Ptr{af_index_t}}(0)
    _error(ccall((:af_create_indexers,af_lib),af_err,(Ptr{Ptr{af_index_t}},),indexers))
    indexers[]
end

function create_handle(ndims::Integer,dims,_type::Type)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_create_handle,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),arr,UInt32(ndims),dims,af_type(_type)))
    AFArray!(arr[])
end

function copy(_in::AFArray{T,N}) where {T,N}
    arr = RefValue{af_array}(0)
    _error(ccall((:af_copy_array,af_lib),af_err,(Ptr{af_array},af_array),arr,_in.arr))
    AFArray{T,N}(arr[])
end

function write_array(arr::AFArray,data,bytes::Csize_t,src::af_source)
    _error(ccall((:af_write_array,af_lib),af_err,(af_array,Ptr{Cvoid},Csize_t,af_source),arr.arr,data,bytes,src))
end

function get_data_ptr(data,arr::AFArray)
    _error(ccall((:af_get_data_ptr,af_lib),af_err,(Ptr{Cvoid},af_array),data,arr.arr))
end

function set_manual_eval_flag(flag::Bool)
    _error(ccall((:af_set_manual_eval_flag,af_lib),af_err,(Bool,),flag))
end

function get_manual_eval_flag()
    flag = RefValue{Bool}(0)
    _error(ccall((:af_get_manual_eval_flag,af_lib),af_err,(Ptr{Bool},),flag))
    flag[]
end

function get_elements(arr::AFArray)
    elems = RefValue{dim_t}(0)
    _error(ccall((:af_get_elements,af_lib),af_err,(Ptr{dim_t},af_array),elems,arr.arr))
    elems[]
end

function get_dims(arr::AFArray)
    d0 = RefValue{dim_t}(0)
    d1 = RefValue{dim_t}(0)
    d2 = RefValue{dim_t}(0)
    d3 = RefValue{dim_t}(0)
    _error(ccall((:af_get_dims,af_lib),af_err,(Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},af_array),d0,d1,d2,d3,arr.arr))
    (d0[],d1[],d2[],d3[])
end

function is_empty(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_empty,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_scalar(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_scalar,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_row(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_row,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_column(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_column,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_vector(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_vector,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_complex(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_complex,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_real(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_real,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_double(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_double,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_single(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_single,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_realfloating(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_realfloating,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_floating(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_floating,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_integer(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_integer,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function is_bool(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_bool,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function issparse(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_sparse,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

function set_backend(bknd::af_backend)
    _error(ccall((:af_set_backend,af_lib),af_err,(af_backend,),bknd))
end

function get_backend_count()
    num_backends = RefValue{UInt32}(0)
    _error(ccall((:af_get_backend_count,af_lib),af_err,(Ptr{UInt32},),num_backends))
    num_backends[]
end

function get_available_backends()
    backends = RefValue{Cint}(0)
    _error(ccall((:af_get_available_backends,af_lib),af_err,(Ptr{Cint},),backends))
    backends[]
end

function get_backend_id(_in::AFArray)
    backend = RefValue{af_backend}(0)
    _error(ccall((:af_get_backend_id,af_lib),af_err,(Ptr{af_backend},af_array),backend,_in.arr))
    backend[]
end

function get_active_backend()
    backend = RefValue{af_backend}(0)
    _error(ccall((:af_get_active_backend,af_lib),af_err,(Ptr{af_backend},),backend))
    backend[]
end

function get_device_id(_in::AFArray)
    device = RefValue{Cint}(0)
    _error(ccall((:af_get_device_id,af_lib),af_err,(Ptr{Cint},af_array),device,_in.arr))
    device[]
end

function matmul(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},optLhs::af_mat_prop,optRhs::af_mat_prop) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs.arr,rhs.arr,optLhs,optRhs))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function dot(lhs::AFArray{T1,N1},rhs::AFArray{T2,N2},optLhs::af_mat_prop,optRhs::af_mat_prop) where {T1,N1,T2,N2}
    out = RefValue{af_array}(0)
    _error(ccall((:af_dot,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs.arr,rhs.arr,optLhs,optRhs))
    AFArray{typed(T1,T2),batched(N1,N2)}(out[])
end

function dot_all(lhs::AFArray,rhs::AFArray,optLhs::af_mat_prop,optRhs::af_mat_prop)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_dot_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array,af_mat_prop,af_mat_prop),real,imag,lhs.arr,rhs.arr,optLhs,optRhs))
    (real[],imag[])
end

function transpose_inplace(_in::AFArray,conjugate::Bool)
    _error(ccall((:af_transpose_inplace,af_lib),af_err,(af_array,Bool),_in.arr,conjugate))
end

function range(ndims::Integer,dims,seq_dim::Integer,_type::Type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_range,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},Cint,af_dtype),out,UInt32(ndims),dims,Cint(seq_dim),af_type(_type)))
    AFArray!(out[])
end

function iota(ndims::Integer,dims,t_ndims::Integer,tdims,_type::Type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_iota,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,UInt32(t_ndims),tdims,af_type(_type)))
    AFArray!(out[])
end

function identity(ndims::Integer,dims,_type::Type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_identity,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,af_type(_type)))
    AFArray!(out[])
end

function diagm(_in::AFArray,num::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diag_create,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(num)))
    AFArray!(out[])
end

function diag(_in::AFArray,num::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diag_extract,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(num)))
    AFArray!(out[])
end

function tile(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_tile,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in.arr,UInt32(x),UInt32(y),UInt32(z),UInt32(w)))
    AFArray{T,N}(out[])
end

function reorder(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_reorder,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in.arr,UInt32(x),UInt32(y),UInt32(z),UInt32(w)))
    AFArray{T,N}(out[])
end

function shift(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_shift,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint,Cint,Cint),out,_in.arr,Cint(x),Cint(y),Cint(z),Cint(w)))
    AFArray{T,N}(out[])
end

function flip(_in::AFArray{T,N},dim::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_flip,af_lib),af_err,(Ptr{af_array},af_array,UInt32),out,_in.arr,UInt32(dim - 1)))
    AFArray{T,N}(out[])
end

function lower(_in::AFArray{T,N},is_unit_diag::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_lower,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_unit_diag))
    AFArray{T,N}(out[])
end

function upper(_in::AFArray{T,N},is_unit_diag::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_upper,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_unit_diag))
    AFArray{T,N}(out[])
end

function replace(a::AFArray,cond::AFArray,b::AFArray)
    _error(ccall((:af_replace,af_lib),af_err,(af_array,af_array,af_array),a.arr,cond.arr,b.arr))
end

function replace(a::AFArray,cond::AFArray,b::Real)
    _error(ccall((:af_replace_scalar,af_lib),af_err,(af_array,af_array,Cdouble),a.arr,cond.arr,Cdouble(b)))
end

function afinfo()
    _error(ccall((:af_info,af_lib),af_err,()))
end

function afinit()
    _error(ccall((:af_init,af_lib),af_err,()))
end

function get_device_count()
    num_of_devices = RefValue{Cint}(0)
    _error(ccall((:af_get_device_count,af_lib),af_err,(Ptr{Cint},),num_of_devices))
    num_of_devices[]
end

function get_dbl_support(device::Integer)
    available = RefValue{Bool}(0)
    _error(ccall((:af_get_dbl_support,af_lib),af_err,(Ptr{Bool},Cint),available,Cint(device)))
    available[]
end

function set_device(device::Integer)
    _error(ccall((:af_set_device,af_lib),af_err,(Cint,),Cint(device)))
end

function get_device()
    device = RefValue{Cint}(0)
    _error(ccall((:af_get_device,af_lib),af_err,(Ptr{Cint},),device))
    device[]
end

function sync(device::Integer)
    _error(ccall((:af_sync,af_lib),af_err,(Cint,),Cint(device)))
end

function alloc_device(bytes::dim_t)
    ptr = RefValue{Ptr{Cvoid}}(0)
    _error(ccall((:af_alloc_device,af_lib),af_err,(Ptr{Ptr{Cvoid}},dim_t),ptr,bytes))
    ptr[]
end

function free_device(ptr)
    _error(ccall((:af_free_device,af_lib),af_err,(Ptr{Cvoid},),ptr))
end

function device_array(data,ndims::Integer,dims,_type::Type)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_device_array,af_lib),af_err,(Ptr{af_array},Ptr{Cvoid},UInt32,Ptr{dim_t},af_dtype),arr,data,UInt32(ndims),dims,af_type(_type)))
    AFArray!(arr[])
end

function print_mem_info(msg,device_id::Integer)
    _error(ccall((:af_print_mem_info,af_lib),af_err,(Cstring,Cint),msg,Cint(device_id)))
end

function set_mem_step_size(step_bytes::Csize_t)
    _error(ccall((:af_set_mem_step_size,af_lib),af_err,(Csize_t,),step_bytes))
end

function get_mem_step_size()
    step_bytes = RefValue{Csize_t}(0)
    _error(ccall((:af_get_mem_step_size,af_lib),af_err,(Ptr{Csize_t},),step_bytes))
    step_bytes[]
end

function lock_device_ptr(arr::AFArray)
    _error(ccall((:af_lock_device_ptr,af_lib),af_err,(af_array,),arr.arr))
end

function unlock_device_ptr(arr::AFArray)
    _error(ccall((:af_unlock_device_ptr,af_lib),af_err,(af_array,),arr.arr))
end

function lock_array(arr::AFArray)
    _error(ccall((:af_lock_array,af_lib),af_err,(af_array,),arr.arr))
end

function unlock_array(arr::AFArray)
    _error(ccall((:af_unlock_array,af_lib),af_err,(af_array,),arr.arr))
end

function is_locked_array(arr::AFArray)
    res = RefValue{Bool}(0)
    _error(ccall((:af_is_locked_array,af_lib),af_err,(Ptr{Bool},af_array),res,arr.arr))
    res[]
end

function get_device_ptr(arr::AFArray)
    ptr = RefValue{Ptr{Cvoid}}(0)
    _error(ccall((:af_get_device_ptr,af_lib),af_err,(Ptr{Ptr{Cvoid}},af_array),ptr,arr.arr))
    ptr[]
end

function create_features(num::dim_t)
    feat = RefValue{af_features}(0)
    _error(ccall((:af_create_features,af_lib),af_err,(Ptr{af_features},dim_t),feat,num))
    feat[]
end

function retain_features(feat::af_features)
    out = RefValue{af_features}(0)
    _error(ccall((:af_retain_features,af_lib),af_err,(Ptr{af_features},af_features),out,feat))
    out[]
end

function get_features_num(feat::af_features)
    num = RefValue{dim_t}(0)
    _error(ccall((:af_get_features_num,af_lib),af_err,(Ptr{dim_t},af_features),num,feat))
    num[]
end

function get_features_xpos(feat::af_features)
    out = RefValue{af_array}(0)
    _error(ccall((:af_get_features_xpos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    AFArray!(out[])
end

function get_features_ypos(feat::af_features)
    out = RefValue{af_array}(0)
    _error(ccall((:af_get_features_ypos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    AFArray!(out[])
end

function get_features_score(feat::af_features)
    score = RefValue{af_array}(0)
    _error(ccall((:af_get_features_score,af_lib),af_err,(Ptr{af_array},af_features),score,feat))
    AFArray!(score[])
end

function get_features_orientation(feat::af_features)
    orientation = RefValue{af_array}(0)
    _error(ccall((:af_get_features_orientation,af_lib),af_err,(Ptr{af_array},af_features),orientation,feat))
    AFArray!(orientation[])
end

function get_features_size(feat::af_features)
    size = RefValue{af_array}(0)
    _error(ccall((:af_get_features_size,af_lib),af_err,(Ptr{af_array},af_features),size,feat))
    AFArray!(size[])
end

function release_features(feat::af_features)
    _error(ccall((:af_release_features,af_lib),af_err,(af_features,),feat))
end

function create_window(width::Integer,height::Integer,title)
    out = RefValue{af_window}(0)
    _error(ccall((:af_create_window,af_lib),af_err,(Ptr{af_window},Cint,Cint,Cstring),out,Cint(width),Cint(height),title))
    out[]
end

function set_position(wind::af_window,x::Integer,y::Integer)
    _error(ccall((:af_set_position,af_lib),af_err,(af_window,UInt32,UInt32),wind,UInt32(x),UInt32(y)))
end

function set_title(wind::af_window,title)
    _error(ccall((:af_set_title,af_lib),af_err,(af_window,Cstring),wind,title))
end

function set_size(wind::af_window,w::Integer,h::Integer)
    _error(ccall((:af_set_size,af_lib),af_err,(af_window,UInt32,UInt32),wind,UInt32(w),UInt32(h)))
end

function draw_image(wind::af_window,_in::AFArray,props)
    _error(ccall((:af_draw_image,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,_in.arr,props))
end

function draw_plot(wind::af_window,X::AFArray,Y::AFArray,props)
    _error(ccall((:af_draw_plot,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,props))
end

function draw_plot3(wind::af_window,P::AFArray,props)
    _error(ccall((:af_draw_plot3,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P.arr,props))
end

function draw_plot_nd(wind::af_window,P::AFArray,props)
    _error(ccall((:af_draw_plot_nd,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P.arr,props))
end

function draw_plot_2d(wind::af_window,X::AFArray,Y::AFArray,props)
    _error(ccall((:af_draw_plot_2d,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,props))
end

function draw_plot_3d(wind::af_window,X::AFArray,Y::AFArray,Z::AFArray,props)
    _error(ccall((:af_draw_plot_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,Z.arr,props))
end

function draw_scatter(wind::af_window,X::AFArray,Y::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,marker,props))
end

function draw_scatter3(wind::af_window,P::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter3,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P.arr,marker,props))
end

function draw_scatter_nd(wind::af_window,P::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_nd,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P.arr,marker,props))
end

function draw_scatter_2d(wind::af_window,X::AFArray,Y::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_2d,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,marker,props))
end

function draw_scatter_3d(wind::af_window,X::AFArray,Y::AFArray,Z::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,Z.arr,marker,props))
end

function draw_hist(wind::af_window,X::AFArray,minval::Real,maxval::Real,props)
    _error(ccall((:af_draw_hist,af_lib),af_err,(af_window,af_array,Cdouble,Cdouble,Ptr{af_cell}),wind,X.arr,Cdouble(minval),Cdouble(maxval),props))
end

function draw_surface(wind::af_window,xVals::AFArray,yVals::AFArray,S::AFArray,props)
    _error(ccall((:af_draw_surface,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,xVals.arr,yVals.arr,S.arr,props))
end

function draw_vector_field_nd(wind::af_window,points::AFArray,directions::AFArray,props)
    _error(ccall((:af_draw_vector_field_nd,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,points.arr,directions.arr,props))
end

function draw_vector_field_3d(wind::af_window,xPoints::AFArray,yPoints::AFArray,zPoints::AFArray,xDirs::AFArray,yDirs::AFArray,zDirs::AFArray,props)
    _error(ccall((:af_draw_vector_field_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints.arr,yPoints.arr,zPoints.arr,xDirs.arr,yDirs.arr,zDirs.arr,props))
end

function draw_vector_field_2d(wind::af_window,xPoints::AFArray,yPoints::AFArray,xDirs::AFArray,yDirs::AFArray,props)
    _error(ccall((:af_draw_vector_field_2d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints.arr,yPoints.arr,xDirs.arr,yDirs.arr,props))
end

function set_axes_limits_compute(wind::af_window,x::AFArray,y::AFArray,z::AFArray,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_compute,af_lib),af_err,(af_window,af_array,af_array,af_array,Bool,Ptr{af_cell}),wind,x.arr,y.arr,z.arr,exact,props))
end

function set_axes_limits_2d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_2d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,exact,props))
end

function set_axes_limits_3d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,zmin::Cfloat,zmax::Cfloat,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_3d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,zmin,zmax,exact,props))
end

function set_axes_titles(wind::af_window,xtitle,ytitle,ztitle,props)
    _error(ccall((:af_set_axes_titles,af_lib),af_err,(af_window,Cstring,Cstring,Cstring,Ptr{af_cell}),wind,xtitle,ytitle,ztitle,props))
end

function show_window(wind::af_window)
    _error(ccall((:af_show,af_lib),af_err,(af_window,),wind))
end

function is_window_closed(wind::af_window)
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_window_closed,af_lib),af_err,(Ptr{Bool},af_window),out,wind))
    out[]
end

function set_visibility(wind::af_window,is_visible::Bool)
    _error(ccall((:af_set_visibility,af_lib),af_err,(af_window,Bool),wind,is_visible))
end

function destroy_window(wind::af_window)
    _error(ccall((:af_destroy_window,af_lib),af_err,(af_window,),wind))
end

function gradient(_in::AFArray)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    _error(ccall((:af_gradient,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array),dx,dy,_in.arr))
    (AFArray!(dx[]),AFArray!(dy[]))
end

function load_image(filename,isColor::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image,af_lib),af_err,(Ptr{af_array},Cstring,Bool),out,filename,isColor))
    AFArray!(out[])
end

function save_image(filename,_in::AFArray)
    _error(ccall((:af_save_image,af_lib),af_err,(Cstring,af_array),filename,_in.arr))
end

function load_image_memory(ptr)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image_memory,af_lib),af_err,(Ptr{af_array},Ptr{Cvoid}),out,ptr))
    AFArray!(out[])
end

function save_image_memory(_in::AFArray,format::af_image_format)
    ptr = RefValue{Ptr{Cvoid}}(0)
    _error(ccall((:af_save_image_memory,af_lib),af_err,(Ptr{Ptr{Cvoid}},af_array,af_image_format),ptr,_in.arr,format))
    ptr[]
end

function delete_image_memory(ptr)
    _error(ccall((:af_delete_image_memory,af_lib),af_err,(Ptr{Cvoid},),ptr))
end

function load_image_native(filename)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image_native,af_lib),af_err,(Ptr{af_array},Cstring),out,filename))
    AFArray!(out[])
end

function save_image_native(filename,_in::AFArray)
    _error(ccall((:af_save_image_native,af_lib),af_err,(Cstring,af_array),filename,_in.arr))
end

function is_image_io_available()
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_image_io_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

function resize(_in::AFArray{T,N},odim0::dim_t,odim1::dim_t,method::af_interp_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_resize,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_interp_type),out,_in.arr,odim0,odim1,method))
    AFArray{T,N}(out[])
end

function transform(_in::AFArray,transform::AFArray,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_transform,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t,dim_t,af_interp_type,Bool),out,_in.arr,transform.arr,odim0,odim1,method,inverse))
    AFArray!(out[])
end

function transform_coordinates(tf::AFArray{T,N},d0::Cfloat,d1::Cfloat) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_transform_coordinates,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat),out,tf.arr,d0,d1))
    AFArray{T,N}(out[])
end

function rotate(_in::AFArray{T,N},theta::Cfloat,crop::Bool,method::af_interp_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_rotate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Bool,af_interp_type),out,_in.arr,theta,crop,method))
    AFArray{T,N}(out[])
end

function translate(_in::AFArray{T,N},trans0::Cfloat,trans1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_translate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in.arr,trans0,trans1,odim0,odim1,method))
    AFArray{T,N}(out[])
end

function scale(_in::AFArray{T,N},scale0::Cfloat,scale1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_scale,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in.arr,scale0,scale1,odim0,odim1,method))
    AFArray{T,N}(out[])
end

function skew(_in::AFArray{T,N},skew0::Cfloat,skew1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_skew,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type,Bool),out,_in.arr,skew0,skew1,odim0,odim1,method,inverse))
    AFArray{T,N}(out[])
end

function histogram(_in::AFArray{T,N},nbins::Integer,minval::Real,maxval::Real) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_histogram,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Cdouble,Cdouble),out,_in.arr,UInt32(nbins),Cdouble(minval),Cdouble(maxval)))
    AFArray{T,N}(out[])
end

function dilate(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dilate,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

function dilate3(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dilate3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

function erode(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_erode,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

function erode3(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_erode3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

function bilateral(_in::AFArray{T,N},spatial_sigma::Cfloat,chromatic_sigma::Cfloat,isColor::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_bilateral,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Bool),out,_in.arr,spatial_sigma,chromatic_sigma,isColor))
    AFArray{T,N}(out[])
end

function mean_shift(_in::AFArray{T,N},spatial_sigma::Cfloat,chromatic_sigma::Cfloat,iter::Integer,is_color::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean_shift,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,UInt32,Bool),out,_in.arr,spatial_sigma,chromatic_sigma,UInt32(iter),is_color))
    AFArray{T,N}(out[])
end

function minfilt(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_minfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

function maxfilt(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_maxfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

function regions(_in::AFArray{T,N},connectivity::af_connectivity,ty::Type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_regions,af_lib),af_err,(Ptr{af_array},af_array,af_connectivity,af_dtype),out,_in.arr,connectivity,af_type(ty)))
    AFArray!(out[])
end

function sobel_operator(img::AFArray,ker_size::Integer)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    _error(ccall((:af_sobel_operator,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,UInt32),dx,dy,img.arr,UInt32(ker_size)))
    (AFArray!(dx[]),AFArray!(dy[]))
end

function rgb2gray(_in::AFArray{T,N},rPercent::Cfloat,gPercent::Cfloat,bPercent::Cfloat) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2gray,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in.arr,rPercent,gPercent,bPercent))
    AFArray{T,N}(out[])
end

function gray2rgb(_in::AFArray{T,N},rFactor::Cfloat,gFactor::Cfloat,bFactor::Cfloat) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_gray2rgb,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in.arr,rFactor,gFactor,bFactor))
    AFArray{T,N}(out[])
end

function hist_equal(_in::AFArray,hist::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_hist_equal,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,hist.arr))
    AFArray!(out[])
end

function gaussian_kernel(rows::Integer,cols::Integer,sigma_r::Real,sigma_c::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_gaussian_kernel,af_lib),af_err,(Ptr{af_array},Cint,Cint,Cdouble,Cdouble),out,Cint(rows),Cint(cols),Cdouble(sigma_r),Cdouble(sigma_c)))
    AFArray!(out[])
end

function hsv2rgb(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_hsv2rgb,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function rgb2hsv(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2hsv,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function color_space(image::AFArray{T,N},to::af_cspace_t,from::af_cspace_t) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_color_space,af_lib),af_err,(Ptr{af_array},af_array,af_cspace_t,af_cspace_t),out,image.arr,to,from))
    AFArray{T,N}(out[])
end

function unwrap(_in::AFArray{T,N},wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_unwrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in.arr,wx,wy,sx,sy,px,py,is_column))
    AFArray{T,N}(out[])
end

function wrap(_in::AFArray{T,N},ox::dim_t,oy::dim_t,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_wrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in.arr,ox,oy,wx,wy,sx,sy,px,py,is_column))
    AFArray{T,N}(out[])
end

function sat(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sat,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

function ycbcr2rgb(_in::AFArray{T,N},standard::af_ycc_std) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_ycbcr2rgb,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in.arr,standard))
    AFArray{T,N}(out[])
end

function rgb2ycbcr(_in::AFArray{T,N},standard::af_ycc_std) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2ycbcr,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in.arr,standard))
    AFArray{T,N}(out[])
end

function moments(_in::AFArray{T,N},moment::af_moment_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_moments,af_lib),af_err,(Ptr{af_array},af_array,af_moment_type),out,_in.arr,moment))
    AFArray{T,N}(out[])
end

function moments_all(_in::AFArray,moment::af_moment_type)
    out = RefValue{Cdouble}(0)
    _error(ccall((:af_moments_all,af_lib),af_err,(Ptr{Cdouble},af_array,af_moment_type),out,_in.arr,moment))
    out[]
end

function canny(_in::AFArray{T,N},threshold_type::af_canny_threshold,low_threshold_ratio::Cfloat,high_threshold_ratio::Cfloat,sobel_window::Integer,is_fast::Bool) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_canny,af_lib),af_err,(Ptr{af_array},af_array,af_canny_threshold,Cfloat,Cfloat,UInt32,Bool),out,_in.arr,threshold_type,low_threshold_ratio,high_threshold_ratio,UInt32(sobel_window),is_fast))
    AFArray{T,N}(out[])
end

function anisotropic_diffusion(_in::AFArray{T,N},timestep::Cfloat,conductance::Cfloat,iterations::Integer,fftype::af_flux_function,diffusion_kind::af_diffusion_eq) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_anisotropic_diffusion,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,UInt32,af_flux_function,af_diffusion_eq),out,_in.arr,timestep,conductance,UInt32(iterations),fftype,diffusion_kind))
    AFArray{T,N}(out[])
end

function svd_inplace(_in::AFArray)
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    _error(ccall((:af_svd_inplace,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),u,s,vt,_in.arr))
    (AFArray!(u[]),AFArray!(s[]),AFArray!(vt[]))
end

function lu(_in::AFArray)
    lower = RefValue{af_array}(0)
    upper = RefValue{af_array}(0)
    pivot = RefValue{af_array}(0)
    _error(ccall((:af_lu,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),lower,upper,pivot,_in.arr))
    (AFArray!(lower[]),AFArray!(upper[]),AFArray!(pivot[]))
end

function lu_inplace(_in::AFArray{T,N},is_lapack_piv::Bool) where {T,N}
    pivot = RefValue{af_array}(0)
    _error(ccall((:af_lu_inplace,af_lib),af_err,(Ptr{af_array},af_array,Bool),pivot,_in.arr,is_lapack_piv))
    AFArray{T,N}(pivot[])
end

function qr(_in::AFArray)
    q = RefValue{af_array}(0)
    r = RefValue{af_array}(0)
    tau = RefValue{af_array}(0)
    _error(ccall((:af_qr,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),q,r,tau,_in.arr))
    (AFArray!(q[]),AFArray!(r[]),AFArray!(tau[]))
end

function qr_inplace(_in::AFArray{T,N}) where {T,N}
    tau = RefValue{af_array}(0)
    _error(ccall((:af_qr_inplace,af_lib),af_err,(Ptr{af_array},af_array),tau,_in.arr))
    AFArray{T,N}(tau[])
end

function cholesky_inplace(_in::AFArray,is_upper::Bool)
    info = RefValue{Cint}(0)
    _error(ccall((:af_cholesky_inplace,af_lib),af_err,(Ptr{Cint},af_array,Bool),info,_in.arr,is_upper))
    info[]
end

function solve(a::AFArray,b::AFArray,options::af_mat_prop)
    x = RefValue{af_array}(0)
    _error(ccall((:af_solve,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop),x,a.arr,b.arr,options))
    AFArray!(x[])
end

function solve_lu(a::AFArray,piv::AFArray,b::AFArray,options::af_mat_prop)
    x = RefValue{af_array}(0)
    _error(ccall((:af_solve_lu,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_mat_prop),x,a.arr,piv.arr,b.arr,options))
    AFArray!(x[])
end

function inverse(_in::AFArray{T,N},options::af_mat_prop) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_inverse,af_lib),af_err,(Ptr{af_array},af_array,af_mat_prop),out,_in.arr,options))
    AFArray{T,N}(out[])
end

function rank(_in::AFArray,tol::Real)
    rank = RefValue{UInt32}(0)
    _error(ccall((:af_rank,af_lib),af_err,(Ptr{UInt32},af_array,Cdouble),rank,_in.arr,Cdouble(tol)))
    rank[]
end

function det(_in::AFArray)
    det_real = RefValue{Cdouble}(0)
    det_imag = RefValue{Cdouble}(0)
    _error(ccall((:af_det,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),det_real,det_imag,_in.arr))
    (det_real[],det_imag[])
end

function norm(_in::AFArray,_type::af_norm_type,p::Real,q::Real)
    out = RefValue{Cdouble}(0)
    _error(ccall((:af_norm,af_lib),af_err,(Ptr{Cdouble},af_array,af_norm_type,Cdouble,Cdouble),out,_in.arr,_type,Cdouble(p),Cdouble(q)))
    out[]
end

function is_lapack_available()
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_lapack_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

function create_random_engine(rtype::af_random_engine_type,seed::uintl)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_create_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type,uintl),engine,rtype,seed))
    engine[]
end

function retain_random_engine(engine::af_random_engine)
    out = RefValue{af_random_engine}(0)
    _error(ccall((:af_retain_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine),out,engine))
    out[]
end

function random_engine_set_type(rtype::af_random_engine_type)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_random_engine_set_type,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type),engine,rtype))
    engine[]
end

function random_engine_get_type(engine::af_random_engine)
    rtype = RefValue{af_random_engine_type}(0)
    _error(ccall((:af_random_engine_get_type,af_lib),af_err,(Ptr{af_random_engine_type},af_random_engine),rtype,engine))
    rtype[]
end

function random_uniform(ndims::Integer,dims,_type::Type,engine::af_random_engine)
    out = RefValue{af_array}(0)
    _error(ccall((:af_random_uniform,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,UInt32(ndims),dims,af_type(_type),engine))
    AFArray!(out[])
end

function random_normal(ndims::Integer,dims,_type::Type,engine::af_random_engine)
    out = RefValue{af_array}(0)
    _error(ccall((:af_random_normal,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,UInt32(ndims),dims,af_type(_type),engine))
    AFArray!(out[])
end

function random_engine_set_seed(seed::uintl)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_random_engine_set_seed,af_lib),af_err,(Ptr{af_random_engine},uintl),engine,seed))
    engine[]
end

function get_default_random_engine()
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_get_default_random_engine,af_lib),af_err,(Ptr{af_random_engine},),engine))
    engine[]
end

function set_default_random_engine_type(rtype::af_random_engine_type)
    _error(ccall((:af_set_default_random_engine_type,af_lib),af_err,(af_random_engine_type,),rtype))
end

function random_engine_get_seed(engine::af_random_engine)
    seed = RefValue{uintl}(0)
    _error(ccall((:af_random_engine_get_seed,af_lib),af_err,(Ptr{uintl},af_random_engine),seed,engine))
    seed[]
end

function release_random_engine(engine::af_random_engine)
    _error(ccall((:af_release_random_engine,af_lib),af_err,(af_random_engine,),engine))
end

function set_seed(seed::uintl)
    _error(ccall((:af_set_seed,af_lib),af_err,(uintl,),seed))
end

function get_seed()
    seed = RefValue{uintl}(0)
    _error(ccall((:af_get_seed,af_lib),af_err,(Ptr{uintl},),seed))
    seed[]
end

function approx1(_in::AFArray,pos::AFArray,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_approx1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_interp_type,Cfloat),out,_in.arr,pos.arr,method,offGrid))
    AFArray!(out[])
end

function approx2(_in::AFArray,pos0::AFArray,pos1::AFArray,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_approx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_interp_type,Cfloat),out,_in.arr,pos0.arr,pos1.arr,method,offGrid))
    AFArray!(out[])
end

function convolve1(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

function convolve2(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

function convolve3(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

function convolve2_sep(col_filter::AFArray,row_filter::AFArray,signal::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve2_sep,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_conv_mode),out,col_filter.arr,row_filter.arr,signal.arr,mode))
    AFArray!(out[])
end

function medfilt(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

function medfilt1(_in::AFArray{T,N},wind_width::dim_t,edge_pad::af_border_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt1,af_lib),af_err,(Ptr{af_array},af_array,dim_t,af_border_type),out,_in.arr,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

function medfilt2(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt2,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

function create_sparse_array(nRows::dim_t,nCols::dim_t,values::AFArray,rowIdx::AFArray,colIdx::AFArray,stype::af_storage)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,af_array,af_array,af_array,af_storage),out,nRows,nCols,values.arr,rowIdx.arr,colIdx.arr,stype))
    AFArray!(out[])
end

function create_sparse_array_from_ptr(nRows::dim_t,nCols::dim_t,nNZ::dim_t,values,rowIdx,colIdx,_type::Type,stype::af_storage,src::af_source)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array_from_ptr,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,dim_t,Ptr{Cvoid},Ptr{Cint},Ptr{Cint},af_dtype,af_storage,af_source),out,nRows,nCols,nNZ,values,rowIdx,colIdx,af_type(_type),stype,src))
    AFArray!(out[])
end

function create_sparse_array_from_dense(dense::AFArray,stype::af_storage)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array_from_dense,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,dense.arr,stype))
    AFArray!(out[])
end

function sparse_convert_to(_in::AFArray{T,N},destStorage::af_storage) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_convert_to,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,_in.arr,destStorage))
    AFArray{T,N}(out[])
end

function full(sparse::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_to_dense,af_lib),af_err,(Ptr{af_array},af_array),out,sparse.arr))
    AFArray{T,N}(out[])
end

function sparse_get_info(_in::AFArray)
    values = RefValue{af_array}(0)
    rowIdx = RefValue{af_array}(0)
    colIdx = RefValue{af_array}(0)
    stype = RefValue{af_storage}(0)
    _error(ccall((:af_sparse_get_info,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},Ptr{af_storage},af_array),values,rowIdx,colIdx,stype,_in.arr))
    (AFArray!(values[]),AFArray!(rowIdx[]),AFArray!(colIdx[]),stype[])
end

function sparse_get_values(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_values,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,1}(out[])
end

function sparse_get_row_idx(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_row_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Int32,1}(out[])
end

function sparse_get_col_idx(_in::AFArray{T,N}) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_col_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{Int32,1}(out[])
end

function sparse_get_nnz(_in::AFArray)
    out = RefValue{dim_t}(0)
    _error(ccall((:af_sparse_get_nnz,af_lib),af_err,(Ptr{dim_t},af_array),out,_in.arr))
    out[]
end

function sparse_get_storage(_in::AFArray)
    out = RefValue{af_storage}(0)
    _error(ccall((:af_sparse_get_storage,af_lib),af_err,(Ptr{af_storage},af_array),out,_in.arr))
    out[]
end

function cov(X::AFArray,Y::AFArray,isbiased::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_cov,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,X.arr,Y.arr,isbiased))
    AFArray!(out[])
end

function mean_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_mean_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function mean_all_weighted(_in::AFArray,weights::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_mean_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),real,imag,_in.arr,weights.arr))
    (real[],imag[])
end

function var_all(_in::AFArray,isbiased::Bool)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_var_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Bool),realVal,imagVal,_in.arr,isbiased))
    (realVal[],imagVal[])
end

function var_all_weighted(_in::AFArray,weights::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_var_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,_in.arr,weights.arr))
    (realVal[],imagVal[])
end

function stdev_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_stdev_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

function median_all(_in::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_median_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),realVal,imagVal,_in.arr))
    (realVal[],imagVal[])
end

function corrcoef(X::AFArray,Y::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_corrcoef,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,X.arr,Y.arr))
    (realVal[],imagVal[])
end

function topk(_in::AFArray,k::Integer,dim::Integer,order::af_topk_function)
    values = RefValue{af_array}(0)
    indices = RefValue{af_array}(0)
    _error(ccall((:af_topk,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint,Cint,af_topk_function),values,indices,_in.arr,Cint(k),Cint(dim - 1),order))
    (AFArray!(values[]),AFArray!(indices[]))
end

function fast(_in::AFArray,thr::Cfloat,arc_length::Integer,non_max::Bool,feature_ratio::Cfloat,edge::Integer)
    out = RefValue{af_features}(0)
    _error(ccall((:af_fast,af_lib),af_err,(Ptr{af_features},af_array,Cfloat,UInt32,Bool,Cfloat,UInt32),out,_in.arr,thr,UInt32(arc_length),non_max,feature_ratio,UInt32(edge)))
    out[]
end

function harris(_in::AFArray,max_corners::Integer,min_response::Cfloat,sigma::Cfloat,block_size::Integer,k_thr::Cfloat)
    out = RefValue{af_features}(0)
    _error(ccall((:af_harris,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,UInt32,Cfloat),out,_in.arr,UInt32(max_corners),min_response,sigma,UInt32(block_size),k_thr))
    out[]
end

function orb(_in::AFArray,fast_thr::Cfloat,max_feat::Integer,scl_fctr::Cfloat,levels::Integer,blur_img::Bool)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_orb,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,Cfloat,UInt32,Cfloat,UInt32,Bool),feat,desc,_in.arr,fast_thr,UInt32(max_feat),scl_fctr,UInt32(levels),blur_img))
    (feat[],AFArray!(desc[]))
end

function sift(_in::AFArray,n_layers::Integer,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_sift,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in.arr,UInt32(n_layers),contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],AFArray!(desc[]))
end

function gloh(_in::AFArray,n_layers::Integer,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_gloh,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in.arr,UInt32(n_layers),contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],AFArray!(desc[]))
end

function hamming_matcher(query::AFArray,train::AFArray,dist_dim::dim_t,n_dist::Integer)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    _error(ccall((:af_hamming_matcher,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32),idx,dist,query.arr,train.arr,dist_dim,UInt32(n_dist)))
    (AFArray!(idx[]),AFArray!(dist[]))
end

function nearest_neighbour(query::AFArray,train::AFArray,dist_dim::dim_t,n_dist::Integer,dist_type::af_match_type)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    _error(ccall((:af_nearest_neighbour,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32,af_match_type),idx,dist,query.arr,train.arr,dist_dim,UInt32(n_dist),dist_type))
    (AFArray!(idx[]),AFArray!(dist[]))
end

function match_template(search_img::AFArray,template_img::AFArray,m_type::af_match_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_match_template,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_match_type),out,search_img.arr,template_img.arr,m_type))
    AFArray!(out[])
end

function susan(_in::AFArray,radius::Integer,diff_thr::Cfloat,geom_thr::Cfloat,feature_ratio::Cfloat,edge::Integer)
    out = RefValue{af_features}(0)
    _error(ccall((:af_susan,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,Cfloat,UInt32),out,_in.arr,UInt32(radius),diff_thr,geom_thr,feature_ratio,UInt32(edge)))
    out[]
end

function dog(_in::AFArray{T,N},radius1::Integer,radius2::Integer) where {T,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_dog,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint),out,_in.arr,Cint(radius1),Cint(radius2)))
    AFArray{T,N}(out[])
end

function homography(x_src::AFArray,y_src::AFArray,x_dst::AFArray,y_dst::AFArray,htype::af_homography_type,inlier_thr::Cfloat,iterations::Integer,otype::Type)
    H = RefValue{af_array}(0)
    inliers = RefValue{Cint}(0)
    _error(ccall((:af_homography,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,af_array,af_array,af_array,af_homography_type,Cfloat,UInt32,af_dtype),H,inliers,x_src.arr,y_src.arr,x_dst.arr,y_dst.arr,htype,inlier_thr,UInt32(iterations),af_type(otype)))
    (AFArray!(H[]),inliers[])
end
