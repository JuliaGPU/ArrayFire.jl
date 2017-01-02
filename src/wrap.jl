# Julia wrapper for header: /usr/include/arrayfire.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function af_sum(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_sum

function af_sum_nan(_in::af_array,dim::Cint,nanval::Cdouble)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sum_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in,dim,nanval))
    out[]
end

export af_sum_nan

function af_product(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_product,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_product

function af_product_nan(_in::af_array,dim::Cint,nanval::Cdouble)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_product_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in,dim,nanval))
    out[]
end

export af_product_nan

function af_min(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_min,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_min

function af_max(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_max,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_max

function af_all_true(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_all_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_all_true

function af_any_true(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_any_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_any_true

function af_count(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_count,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_count

function af_sum_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_sum_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_sum_all

function af_sum_nan_all(_in::af_array,nanval::Cdouble)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_sum_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in,nanval))
    (real[],imag[])
end

export af_sum_nan_all

function af_product_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_product_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_product_all

function af_product_nan_all(_in::af_array,nanval::Cdouble)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_product_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in,nanval))
    (real[],imag[])
end

export af_product_nan_all

function af_min_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_min_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_min_all

function af_max_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_max_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_max_all

function af_all_true_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_all_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_all_true_all

function af_any_true_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_any_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_any_true_all

function af_count_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_count_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_count_all

function af_imin(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    af_error(ccall((:af_imin,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in,dim))
    (out[],idx[])
end

export af_imin

function af_imax(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    af_error(ccall((:af_imax,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in,dim))
    (out[],idx[])
end

export af_imax

function af_imin_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    af_error(ccall((:af_imin_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in))
    (real[],imag[],idx[])
end

export af_imin_all

function af_imax_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    af_error(ccall((:af_imax_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in))
    (real[],imag[],idx[])
end

export af_imax_all

function af_accum(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_accum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_accum

function af_scan(_in::af_array,dim::Cint,op::af_binary_op,inclusive_scan::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_scan,af_lib),af_err,(Ptr{af_array},af_array,Cint,af_binary_op,Bool),out,_in,dim,op,inclusive_scan))
    out[]
end

export af_scan

function af_scan_by_key(key::af_array,_in::af_array,dim::Cint,op::af_binary_op,inclusive_scan::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_scan_by_key,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cint,af_binary_op,Bool),out,key,_in,dim,op,inclusive_scan))
    out[]
end

export af_scan_by_key

function af_where(_in::af_array)
    idx = RefValue{af_array}(0)
    af_error(ccall((:af_where,af_lib),af_err,(Ptr{af_array},af_array),idx,_in))
    idx[]
end

export af_where

function af_diff1(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_diff1,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_diff1

function af_diff2(_in::af_array,dim::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_diff2,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,dim))
    out[]
end

export af_diff2

function af_sort(_in::af_array,dim::UInt32,isAscending::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sort,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Bool),out,_in,dim,isAscending))
    out[]
end

export af_sort

function af_sort_index(_in::af_array,dim::UInt32,isAscending::Bool)
    out = RefValue{af_array}(0)
    indices = RefValue{af_array}(0)
    af_error(ccall((:af_sort_index,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,UInt32,Bool),out,indices,_in,dim,isAscending))
    (out[],indices[])
end

export af_sort_index

function af_sort_by_key(keys::af_array,values::af_array,dim::UInt32,isAscending::Bool)
    out_keys = RefValue{af_array}(0)
    out_values = RefValue{af_array}(0)
    af_error(ccall((:af_sort_by_key,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,UInt32,Bool),out_keys,out_values,keys,values,dim,isAscending))
    (out_keys[],out_values[])
end

export af_sort_by_key

function af_set_unique(_in::af_array,is_sorted::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_set_unique,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in,is_sorted))
    out[]
end

export af_set_unique

function af_set_union(first::af_array,second::af_array,is_unique::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_set_union,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first,second,is_unique))
    out[]
end

export af_set_union

function af_set_intersect(first::af_array,second::af_array,is_unique::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_set_intersect,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first,second,is_unique))
    out[]
end

export af_set_intersect

function af_add(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_add,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_add

function af_sub(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sub,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_sub

function af_mul(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_mul,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_mul

function af_div(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_div,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_div

function af_lt(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_lt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_lt

function af_gt(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_gt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_gt

function af_le(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_le,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_le

function af_ge(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ge,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_ge

function af_eq(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_eq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_eq

function af_neq(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_neq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_neq

function af_and(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_and,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_and

function af_or(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_or,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_or

function af_not(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_not,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_not

function af_bitand(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bitand,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_bitand

function af_bitor(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bitor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_bitor

function af_bitxor(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bitxor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_bitxor

function af_bitshiftl(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bitshiftl,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_bitshiftl

function af_bitshiftr(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bitshiftr,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_bitshiftr

function af_cast(_in::af_array,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cast,af_lib),af_err,(Ptr{af_array},af_array,af_dtype),out,_in,_type))
    out[]
end

export af_cast

function af_minof(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_minof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_minof

function af_maxof(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_maxof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_maxof

function af_clamp(_in::af_array,lo::af_array,hi::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_clamp,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,Bool),out,_in,lo,hi,batch))
    out[]
end

export af_clamp

function af_rem(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_rem,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_rem

function af_mod(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_mod,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_mod

function af_abs(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_abs,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_abs

function af_arg(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_arg,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_arg

function af_sign(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sign,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sign

function af_round(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_round,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_round

function af_trunc(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_trunc,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_trunc

function af_floor(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_floor,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_floor

function af_ceil(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ceil,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_ceil

function af_hypot(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_hypot,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_hypot

function af_sin(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sin,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sin

function af_cos(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cos,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_cos

function af_tan(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_tan,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_tan

function af_asin(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_asin,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_asin

function af_acos(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_acos,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_acos

function af_atan(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_atan,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_atan

function af_atan2(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_atan2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_atan2

function af_cplx2(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cplx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_cplx2

function af_cplx(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cplx,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_cplx

function af_real(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_real,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_real

function af_imag(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_imag,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_imag

function af_conjg(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_conjg,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_conjg

function af_sinh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sinh

function af_cosh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_cosh

function af_tanh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_tanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_tanh

function af_asinh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_asinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_asinh

function af_acosh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_acosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_acosh

function af_atanh(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_atanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_atanh

function af_root(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_root,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_root

function af_pow(lhs::af_array,rhs::af_array,batch::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_pow,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs,rhs,batch))
    out[]
end

export af_pow

function af_pow2(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_pow2,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_pow2

function af_exp(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_exp,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_exp

function af_sigmoid(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sigmoid,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sigmoid

function af_expm1(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_expm1,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_expm1

function af_erf(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_erf,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_erf

function af_erfc(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_erfc,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_erfc

function af_log(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_log,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_log

function af_log1p(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_log1p,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_log1p

function af_log10(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_log10,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_log10

function af_log2(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_log2,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_log2

function af_sqrt(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sqrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sqrt

function af_cbrt(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cbrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_cbrt

function af_factorial(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_factorial,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_factorial

function af_tgamma(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_tgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_tgamma

function af_lgamma(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_lgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_lgamma

function af_iszero(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_iszero,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_iszero

function af_isinf(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_isinf,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_isinf

function af_isnan(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_isnan,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_isnan

function af_make_seq(_begin::Cdouble,_end::Cdouble,step::Cdouble)
    ccall((:af_make_seq,af_lib),af_seq,(Cdouble,Cdouble,Cdouble),_begin,_end,step)
end

export af_make_seq

function af_print_array(arr::af_array)
    af_error(ccall((:af_print_array,af_lib),af_err,(af_array,),arr))
end

export af_print_array

function af_print_array_gen(exp,arr::af_array,precision::Cint)
    af_error(ccall((:af_print_array_gen,af_lib),af_err,(Cstring,af_array,Cint),exp,arr,precision))
end

export af_print_array_gen

function af_save_array(key,arr::af_array,filename,append::Bool)
    index = RefValue{Cint}(0)
    af_error(ccall((:af_save_array,af_lib),af_err,(Ptr{Cint},Cstring,af_array,Cstring,Bool),index,key,arr,filename,append))
    index[]
end

export af_save_array

function af_read_array_index(filename,index::UInt32)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_read_array_index,af_lib),af_err,(Ptr{af_array},Cstring,UInt32),out,filename,index))
    out[]
end

export af_read_array_index

function af_read_array_key(filename,key)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_read_array_key,af_lib),af_err,(Ptr{af_array},Cstring,Cstring),out,filename,key))
    out[]
end

export af_read_array_key

function af_read_array_key_check(filename,key)
    index = RefValue{Cint}(0)
    af_error(ccall((:af_read_array_key_check,af_lib),af_err,(Ptr{Cint},Cstring,Cstring),index,filename,key))
    index[]
end

export af_read_array_key_check

function af_array_to_string(exp,arr::af_array,precision::Cint,transpose::Bool)
    output = RefValue{Cstring}(0)
    af_error(ccall((:af_array_to_string,af_lib),af_err,(Ptr{Cstring},Cstring,af_array,Cint,Bool),output,exp,arr,precision,transpose))
    output[]
end

export af_array_to_string

function af_example_function(_in::af_array,param::af_someenum_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_example_function,af_lib),af_err,(Ptr{af_array},af_array,af_someenum_t),out,_in,param))
    out[]
end

export af_example_function

function af_get_version()
    major = RefValue{Cint}(0)
    minor = RefValue{Cint}(0)
    patch = RefValue{Cint}(0)
    af_error(ccall((:af_get_version,af_lib),af_err,(Ptr{Cint},Ptr{Cint},Ptr{Cint}),major,minor,patch))
    (major[],minor[],patch[])
end

export af_get_version

function af_get_revision()
    ccall((:af_get_revision,af_lib),Cstring,())
end

export af_get_revision

function af_get_size_of(_type::af_dtype)
    size = RefValue{Csize_t}(0)
    af_error(ccall((:af_get_size_of,af_lib),af_err,(Ptr{Csize_t},af_dtype),size,_type))
    size[]
end

export af_get_size_of

function af_index(_in::af_array,ndims::UInt32,index)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_index,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq}),out,_in,ndims,index))
    out[]
end

export af_index

function af_lookup(_in::af_array,indices::af_array,dim::UInt32)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_lookup,af_lib),af_err,(Ptr{af_array},af_array,af_array,UInt32),out,_in,indices,dim))
    out[]
end

export af_lookup

function af_assign_seq(lhs::af_array,ndims::UInt32,indices,rhs::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_assign_seq,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq},af_array),out,lhs,ndims,indices,rhs))
    out[]
end

export af_assign_seq

function af_index_gen(_in::af_array,ndims::dim_t,indices)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_index_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t}),out,_in,ndims,indices))
    out[]
end

export af_index_gen

function af_assign_gen(lhs::af_array,ndims::dim_t,indices,rhs::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_assign_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t},af_array),out,lhs,ndims,indices,rhs))
    out[]
end

export af_assign_gen

function af_create_indexers()
    indexers = RefValue{Ptr{af_index_t}}(0)
    af_error(ccall((:af_create_indexers,af_lib),af_err,(Ptr{Ptr{af_index_t}},),indexers))
    indexers[]
end

export af_create_indexers

function af_set_array_indexer(idx::af_array,dim::dim_t)
    indexer = RefValue{af_index_t}(0)
    af_error(ccall((:af_set_array_indexer,af_lib),af_err,(Ptr{af_index_t},af_array,dim_t),indexer,idx,dim))
    indexer[]
end

export af_set_array_indexer

function af_set_seq_indexer(dim::dim_t,is_batch::Bool)
    indexer = RefValue{af_index_t}(0)
    idx = RefValue{af_seq}(0)
    af_error(ccall((:af_set_seq_indexer,af_lib),af_err,(Ptr{af_index_t},Ptr{af_seq},dim_t,Bool),indexer,idx,dim,is_batch))
    (indexer[],idx[])
end

export af_set_seq_indexer

function af_set_seq_param_indexer(_begin::Cdouble,_end::Cdouble,step::Cdouble,dim::dim_t,is_batch::Bool)
    indexer = RefValue{af_index_t}(0)
    af_error(ccall((:af_set_seq_param_indexer,af_lib),af_err,(Ptr{af_index_t},Cdouble,Cdouble,Cdouble,dim_t,Bool),indexer,_begin,_end,step,dim,is_batch))
    indexer[]
end

export af_set_seq_param_indexer

function af_release_indexers()
    indexers = RefValue{af_index_t}(0)
    af_error(ccall((:af_release_indexers,af_lib),af_err,(Ptr{af_index_t},),indexers))
    indexers[]
end

export af_release_indexers

function af_create_array(ndims::UInt32,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    data = RefValue{Void}(0)
    af_error(ccall((:af_create_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),arr,data,ndims,dims,_type))
    (arr[],data[])
end

export af_create_array

function af_create_handle(ndims::UInt32,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_create_handle,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),arr,ndims,dims,_type))
    arr[]
end

export af_create_handle

function af_copy_array(_in::af_array)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_copy_array,af_lib),af_err,(Ptr{af_array},af_array),arr,_in))
    arr[]
end

export af_copy_array

function af_write_array(arr::af_array,data,bytes::Csize_t,src::af_source)
    af_error(ccall((:af_write_array,af_lib),af_err,(af_array,Ptr{Void},Csize_t,af_source),arr,data,bytes,src))
end

export af_write_array

function af_get_data_ptr(arr::af_array)
    data = RefValue{Void}(0)
    af_error(ccall((:af_get_data_ptr,af_lib),af_err,(Ptr{Void},af_array),data,arr))
    data[]
end

export af_get_data_ptr

function af_release_array(arr::af_array)
    af_error(ccall((:af_release_array,af_lib),af_err,(af_array,),arr))
end

export af_release_array

function af_retain_array(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_retain_array,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_retain_array

function af_get_data_ref_count(_in::af_array)
    use_count = RefValue{Cint}(0)
    af_error(ccall((:af_get_data_ref_count,af_lib),af_err,(Ptr{Cint},af_array),use_count,_in))
    use_count[]
end

export af_get_data_ref_count

function af_eval(_in::af_array)
    af_error(ccall((:af_eval,af_lib),af_err,(af_array,),_in))
end

export af_eval

function af_eval_multiple(num::Cint,arrays)
    af_error(ccall((:af_eval_multiple,af_lib),af_err,(Cint,Ptr{af_array}),num,arrays))
end

export af_eval_multiple

function af_set_manual_eval_flag(flag::Bool)
    af_error(ccall((:af_set_manual_eval_flag,af_lib),af_err,(Bool,),flag))
end

export af_set_manual_eval_flag

function af_get_manual_eval_flag()
    flag = RefValue{Bool}(0)
    af_error(ccall((:af_get_manual_eval_flag,af_lib),af_err,(Ptr{Bool},),flag))
    flag[]
end

export af_get_manual_eval_flag

function af_get_elements(arr::af_array)
    elems = RefValue{dim_t}(0)
    af_error(ccall((:af_get_elements,af_lib),af_err,(Ptr{dim_t},af_array),elems,arr))
    elems[]
end

export af_get_elements

function af_get_type(arr::af_array)
    _type = RefValue{af_dtype}(0)
    af_error(ccall((:af_get_type,af_lib),af_err,(Ptr{af_dtype},af_array),_type,arr))
    _type[]
end

export af_get_type

function af_get_dims(arr::af_array)
    d0 = RefValue{dim_t}(0)
    d1 = RefValue{dim_t}(0)
    d2 = RefValue{dim_t}(0)
    d3 = RefValue{dim_t}(0)
    af_error(ccall((:af_get_dims,af_lib),af_err,(Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},af_array),d0,d1,d2,d3,arr))
    (d0[],d1[],d2[],d3[])
end

export af_get_dims

function af_get_numdims(arr::af_array)
    result = RefValue{UInt32}(0)
    af_error(ccall((:af_get_numdims,af_lib),af_err,(Ptr{UInt32},af_array),result,arr))
    result[]
end

export af_get_numdims

function af_is_empty(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_empty,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_empty

function af_is_scalar(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_scalar,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_scalar

function af_is_row(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_row,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_row

function af_is_column(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_column,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_column

function af_is_vector(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_vector,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_vector

function af_is_complex(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_complex,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_complex

function af_is_real(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_real,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_real

function af_is_double(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_double,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_double

function af_is_single(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_single,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_single

function af_is_realfloating(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_realfloating,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_realfloating

function af_is_floating(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_floating,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_floating

function af_is_integer(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_integer,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_integer

function af_is_bool(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_bool,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_bool

function af_is_sparse(arr::af_array)
    result = RefValue{Bool}(0)
    af_error(ccall((:af_is_sparse,af_lib),af_err,(Ptr{Bool},af_array),result,arr))
    result[]
end

export af_is_sparse

function af_set_backend(bknd::af_backend)
    af_error(ccall((:af_set_backend,af_lib),af_err,(af_backend,),bknd))
end

export af_set_backend

function af_get_backend_count()
    num_backends = RefValue{UInt32}(0)
    af_error(ccall((:af_get_backend_count,af_lib),af_err,(Ptr{UInt32},),num_backends))
    num_backends[]
end

export af_get_backend_count

function af_get_available_backends()
    backends = RefValue{Cint}(0)
    af_error(ccall((:af_get_available_backends,af_lib),af_err,(Ptr{Cint},),backends))
    backends[]
end

export af_get_available_backends

function af_get_backend_id(_in::af_array)
    backend = RefValue{af_backend}(0)
    af_error(ccall((:af_get_backend_id,af_lib),af_err,(Ptr{af_backend},af_array),backend,_in))
    backend[]
end

export af_get_backend_id

function af_get_active_backend()
    backend = RefValue{af_backend}(0)
    af_error(ccall((:af_get_active_backend,af_lib),af_err,(Ptr{af_backend},),backend))
    backend[]
end

export af_get_active_backend

function af_get_device_id(_in::af_array)
    device = RefValue{Cint}(0)
    af_error(ccall((:af_get_device_id,af_lib),af_err,(Ptr{Cint},af_array),device,_in))
    device[]
end

export af_get_device_id

function af_matmul(lhs::af_array,rhs::af_array,optLhs::af_mat_prop,optRhs::af_mat_prop)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs,rhs,optLhs,optRhs))
    out[]
end

export af_matmul

function af_dot(lhs::af_array,rhs::af_array,optLhs::af_mat_prop,optRhs::af_mat_prop)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_dot,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs,rhs,optLhs,optRhs))
    out[]
end

export af_dot

function af_dot_all(lhs::af_array,rhs::af_array,optLhs::af_mat_prop,optRhs::af_mat_prop)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_dot_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array,af_mat_prop,af_mat_prop),real,imag,lhs,rhs,optLhs,optRhs))
    (real[],imag[])
end

export af_dot_all

function af_transpose(_in::af_array,conjugate::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_transpose,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in,conjugate))
    out[]
end

export af_transpose

function af_transpose_inplace(_in::af_array,conjugate::Bool)
    af_error(ccall((:af_transpose_inplace,af_lib),af_err,(af_array,Bool),_in,conjugate))
end

export af_transpose_inplace

function af_constant(val::Cdouble,ndims::UInt32,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_constant,af_lib),af_err,(Ptr{af_array},Cdouble,UInt32,Ptr{dim_t},af_dtype),arr,val,ndims,dims,_type))
    arr[]
end

export af_constant

function af_constant_complex(real::Cdouble,imag::Cdouble,ndims::UInt32,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_constant_complex,af_lib),af_err,(Ptr{af_array},Cdouble,Cdouble,UInt32,Ptr{dim_t},af_dtype),arr,real,imag,ndims,dims,_type))
    arr[]
end

export af_constant_complex

function af_constant_long(val::intl,ndims::UInt32,dims)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_constant_long,af_lib),af_err,(Ptr{af_array},intl,UInt32,Ptr{dim_t}),arr,val,ndims,dims))
    arr[]
end

export af_constant_long

function af_constant_ulong(val::uintl,ndims::UInt32,dims)
    arr = RefValue{af_array}(0)
    af_error(ccall((:af_constant_ulong,af_lib),af_err,(Ptr{af_array},uintl,UInt32,Ptr{dim_t}),arr,val,ndims,dims))
    arr[]
end

export af_constant_ulong

function af_range(ndims::UInt32,dims,seq_dim::Cint,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_range,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},Cint,af_dtype),out,ndims,dims,seq_dim,_type))
    out[]
end

export af_range

function af_iota(ndims::UInt32,dims,t_ndims::UInt32,tdims,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_iota,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},UInt32,Ptr{dim_t},af_dtype),out,ndims,dims,t_ndims,tdims,_type))
    out[]
end

export af_iota

function af_identity(ndims::UInt32,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_identity,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,ndims,dims,_type))
    out[]
end

export af_identity

function af_diag_create(_in::af_array,num::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_diag_create,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,num))
    out[]
end

export af_diag_create

function af_diag_extract(_in::af_array,num::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_diag_extract,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in,num))
    out[]
end

export af_diag_extract

function af_join(dim::Cint,first::af_array,second::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_join,af_lib),af_err,(Ptr{af_array},Cint,af_array,af_array),out,dim,first,second))
    out[]
end

export af_join

function af_join_many(dim::Cint,n_arrays::UInt32,inputs)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_join_many,af_lib),af_err,(Ptr{af_array},Cint,UInt32,Ptr{af_array}),out,dim,n_arrays,inputs))
    out[]
end

export af_join_many

function af_tile(_in::af_array,x::UInt32,y::UInt32,z::UInt32,w::UInt32)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_tile,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in,x,y,z,w))
    out[]
end

export af_tile

function af_reorder(_in::af_array,x::UInt32,y::UInt32,z::UInt32,w::UInt32)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_reorder,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in,x,y,z,w))
    out[]
end

export af_reorder

function af_shift(_in::af_array,x::Cint,y::Cint,z::Cint,w::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_shift,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint,Cint,Cint),out,_in,x,y,z,w))
    out[]
end

export af_shift

function af_moddims(_in::af_array,ndims::UInt32,dims)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_moddims,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{dim_t}),out,_in,ndims,dims))
    out[]
end

export af_moddims

function af_flat(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_flat,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_flat

function af_flip(_in::af_array,dim::UInt32)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_flip,af_lib),af_err,(Ptr{af_array},af_array,UInt32),out,_in,dim))
    out[]
end

export af_flip

function af_lower(_in::af_array,is_unit_diag::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_lower,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in,is_unit_diag))
    out[]
end

export af_lower

function af_upper(_in::af_array,is_unit_diag::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_upper,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in,is_unit_diag))
    out[]
end

export af_upper

function af_select(cond::af_array,a::af_array,b::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_select,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),out,cond,a,b))
    out[]
end

export af_select

function af_select_scalar_r(cond::af_array,a::af_array,b::Cdouble)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_select_scalar_r,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cdouble),out,cond,a,b))
    out[]
end

export af_select_scalar_r

function af_select_scalar_l(cond::af_array,a::Cdouble,b::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_select_scalar_l,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,af_array),out,cond,a,b))
    out[]
end

export af_select_scalar_l

function af_replace(a::af_array,cond::af_array,b::af_array)
    af_error(ccall((:af_replace,af_lib),af_err,(af_array,af_array,af_array),a,cond,b))
end

export af_replace

function af_replace_scalar(a::af_array,cond::af_array,b::Cdouble)
    af_error(ccall((:af_replace_scalar,af_lib),af_err,(af_array,af_array,Cdouble),a,cond,b))
end

export af_replace_scalar

function af_info()
    af_error(ccall((:af_info,af_lib),af_err,()))
end

export af_info

function af_init()
    af_error(ccall((:af_init,af_lib),af_err,()))
end

export af_init

function af_info_string(verbose::Bool)
    str = RefValue{Cstring}(0)
    af_error(ccall((:af_info_string,af_lib),af_err,(Ptr{Cstring},Bool),str,verbose))
    str[]
end

export af_info_string

function af_device_info(d_name,d_platform,d_toolkit,d_compute)
    af_error(ccall((:af_device_info,af_lib),af_err,(Cstring,Cstring,Cstring,Cstring),d_name,d_platform,d_toolkit,d_compute))
end

export af_device_info

function af_get_device_count()
    num_of_devices = RefValue{Cint}(0)
    af_error(ccall((:af_get_device_count,af_lib),af_err,(Ptr{Cint},),num_of_devices))
    num_of_devices[]
end

export af_get_device_count

function af_get_dbl_support(device::Cint)
    available = RefValue{Bool}(0)
    af_error(ccall((:af_get_dbl_support,af_lib),af_err,(Ptr{Bool},Cint),available,device))
    available[]
end

export af_get_dbl_support

function af_set_device(device::Cint)
    af_error(ccall((:af_set_device,af_lib),af_err,(Cint,),device))
end

export af_set_device

function af_get_device()
    device = RefValue{Cint}(0)
    af_error(ccall((:af_get_device,af_lib),af_err,(Ptr{Cint},),device))
    device[]
end

export af_get_device

function af_sync(device::Cint)
    af_error(ccall((:af_sync,af_lib),af_err,(Cint,),device))
end

export af_sync

function af_alloc_device(bytes::dim_t)
    ptr = RefValue{Ptr{Void}}(0)
    af_error(ccall((:af_alloc_device,af_lib),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes))
    ptr[]
end

export af_alloc_device

function af_free_device()
    ptr = RefValue{Void}(0)
    af_error(ccall((:af_free_device,af_lib),af_err,(Ptr{Void},),ptr))
    ptr[]
end

export af_free_device

function af_alloc_pinned(bytes::dim_t)
    ptr = RefValue{Ptr{Void}}(0)
    af_error(ccall((:af_alloc_pinned,af_lib),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes))
    ptr[]
end

export af_alloc_pinned

function af_free_pinned()
    ptr = RefValue{Void}(0)
    af_error(ccall((:af_free_pinned,af_lib),af_err,(Ptr{Void},),ptr))
    ptr[]
end

export af_free_pinned

function af_alloc_host(bytes::dim_t)
    ptr = RefValue{Ptr{Void}}(0)
    af_error(ccall((:af_alloc_host,af_lib),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes))
    ptr[]
end

export af_alloc_host

function af_free_host()
    ptr = RefValue{Void}(0)
    af_error(ccall((:af_free_host,af_lib),af_err,(Ptr{Void},),ptr))
    ptr[]
end

export af_free_host

function af_device_array(ndims::UInt32,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    data = RefValue{Void}(0)
    af_error(ccall((:af_device_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),arr,data,ndims,dims,_type))
    (arr[],data[])
end

export af_device_array

function af_device_mem_info()
    alloc_bytes = RefValue{Csize_t}(0)
    alloc_buffers = RefValue{Csize_t}(0)
    lock_bytes = RefValue{Csize_t}(0)
    lock_buffers = RefValue{Csize_t}(0)
    af_error(ccall((:af_device_mem_info,af_lib),af_err,(Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t}),alloc_bytes,alloc_buffers,lock_bytes,lock_buffers))
    (alloc_bytes[],alloc_buffers[],lock_bytes[],lock_buffers[])
end

export af_device_mem_info

function af_print_mem_info(msg,device_id::Cint)
    af_error(ccall((:af_print_mem_info,af_lib),af_err,(Cstring,Cint),msg,device_id))
end

export af_print_mem_info

function af_device_gc()
    af_error(ccall((:af_device_gc,af_lib),af_err,()))
end

export af_device_gc

function af_set_mem_step_size(step_bytes::Csize_t)
    af_error(ccall((:af_set_mem_step_size,af_lib),af_err,(Csize_t,),step_bytes))
end

export af_set_mem_step_size

function af_get_mem_step_size()
    step_bytes = RefValue{Csize_t}(0)
    af_error(ccall((:af_get_mem_step_size,af_lib),af_err,(Ptr{Csize_t},),step_bytes))
    step_bytes[]
end

export af_get_mem_step_size

function af_lock_device_ptr(arr::af_array)
    af_error(ccall((:af_lock_device_ptr,af_lib),af_err,(af_array,),arr))
end

export af_lock_device_ptr

function af_unlock_device_ptr(arr::af_array)
    af_error(ccall((:af_unlock_device_ptr,af_lib),af_err,(af_array,),arr))
end

export af_unlock_device_ptr

function af_lock_array(arr::af_array)
    af_error(ccall((:af_lock_array,af_lib),af_err,(af_array,),arr))
end

export af_lock_array

function af_unlock_array(arr::af_array)
    af_error(ccall((:af_unlock_array,af_lib),af_err,(af_array,),arr))
end

export af_unlock_array

function af_is_locked_array(arr::af_array)
    res = RefValue{Bool}(0)
    af_error(ccall((:af_is_locked_array,af_lib),af_err,(Ptr{Bool},af_array),res,arr))
    res[]
end

export af_is_locked_array

function af_get_device_ptr(arr::af_array)
    ptr = RefValue{Ptr{Void}}(0)
    af_error(ccall((:af_get_device_ptr,af_lib),af_err,(Ptr{Ptr{Void}},af_array),ptr,arr))
    ptr[]
end

export af_get_device_ptr

function af_get_last_error()
    msg = RefValue{Cstring}(0)
    len = RefValue{dim_t}(0)
    ccall((:af_get_last_error,af_lib),Void,(Ptr{Cstring},Ptr{dim_t}),msg,len)
    (msg[],len[])
end

export af_get_last_error

function af_err_to_string(err::af_err)
    ccall((:af_err_to_string,af_lib),Cstring,(af_err,),err)
end

export af_err_to_string

function af_create_features(num::dim_t)
    feat = RefValue{af_features}(0)
    af_error(ccall((:af_create_features,af_lib),af_err,(Ptr{af_features},dim_t),feat,num))
    feat[]
end

export af_create_features

function af_retain_features(feat::af_features)
    out = RefValue{af_features}(0)
    af_error(ccall((:af_retain_features,af_lib),af_err,(Ptr{af_features},af_features),out,feat))
    out[]
end

export af_retain_features

function af_get_features_num(feat::af_features)
    num = RefValue{dim_t}(0)
    af_error(ccall((:af_get_features_num,af_lib),af_err,(Ptr{dim_t},af_features),num,feat))
    num[]
end

export af_get_features_num

function af_get_features_xpos(feat::af_features)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_get_features_xpos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    out[]
end

export af_get_features_xpos

function af_get_features_ypos(feat::af_features)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_get_features_ypos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    out[]
end

export af_get_features_ypos

function af_get_features_score(feat::af_features)
    score = RefValue{af_array}(0)
    af_error(ccall((:af_get_features_score,af_lib),af_err,(Ptr{af_array},af_features),score,feat))
    score[]
end

export af_get_features_score

function af_get_features_orientation(feat::af_features)
    orientation = RefValue{af_array}(0)
    af_error(ccall((:af_get_features_orientation,af_lib),af_err,(Ptr{af_array},af_features),orientation,feat))
    orientation[]
end

export af_get_features_orientation

function af_get_features_size(feat::af_features)
    size = RefValue{af_array}(0)
    af_error(ccall((:af_get_features_size,af_lib),af_err,(Ptr{af_array},af_features),size,feat))
    size[]
end

export af_get_features_size

function af_release_features(feat::af_features)
    af_error(ccall((:af_release_features,af_lib),af_err,(af_features,),feat))
end

export af_release_features

function af_create_window(width::Cint,height::Cint,title)
    out = RefValue{af_window}(0)
    af_error(ccall((:af_create_window,af_lib),af_err,(Ptr{af_window},Cint,Cint,Cstring),out,width,height,title))
    out[]
end

export af_create_window

function af_set_position(wind::af_window,x::UInt32,y::UInt32)
    af_error(ccall((:af_set_position,af_lib),af_err,(af_window,UInt32,UInt32),wind,x,y))
end

export af_set_position

function af_set_title(wind::af_window,title)
    af_error(ccall((:af_set_title,af_lib),af_err,(af_window,Cstring),wind,title))
end

export af_set_title

function af_set_size(wind::af_window,w::UInt32,h::UInt32)
    af_error(ccall((:af_set_size,af_lib),af_err,(af_window,UInt32,UInt32),wind,w,h))
end

export af_set_size

function af_draw_image(wind::af_window,_in::af_array,props)
    af_error(ccall((:af_draw_image,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,_in,props))
end

export af_draw_image

function af_draw_plot(wind::af_window,X::af_array,Y::af_array,props)
    af_error(ccall((:af_draw_plot,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X,Y,props))
end

export af_draw_plot

function af_draw_plot3(wind::af_window,P::af_array,props)
    af_error(ccall((:af_draw_plot3,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P,props))
end

export af_draw_plot3

function af_draw_plot_nd(wind::af_window,P::af_array,props)
    af_error(ccall((:af_draw_plot_nd,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P,props))
end

export af_draw_plot_nd

function af_draw_plot_2d(wind::af_window,X::af_array,Y::af_array,props)
    af_error(ccall((:af_draw_plot_2d,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X,Y,props))
end

export af_draw_plot_2d

function af_draw_plot_3d(wind::af_window,X::af_array,Y::af_array,Z::af_array,props)
    af_error(ccall((:af_draw_plot_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,X,Y,Z,props))
end

export af_draw_plot_3d

function af_draw_scatter(wind::af_window,X::af_array,Y::af_array,marker::af_marker_type,props)
    af_error(ccall((:af_draw_scatter,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X,Y,marker,props))
end

export af_draw_scatter

function af_draw_scatter3(wind::af_window,P::af_array,marker::af_marker_type,props)
    af_error(ccall((:af_draw_scatter3,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P,marker,props))
end

export af_draw_scatter3

function af_draw_scatter_nd(wind::af_window,P::af_array,marker::af_marker_type,props)
    af_error(ccall((:af_draw_scatter_nd,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P,marker,props))
end

export af_draw_scatter_nd

function af_draw_scatter_2d(wind::af_window,X::af_array,Y::af_array,marker::af_marker_type,props)
    af_error(ccall((:af_draw_scatter_2d,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X,Y,marker,props))
end

export af_draw_scatter_2d

function af_draw_scatter_3d(wind::af_window,X::af_array,Y::af_array,Z::af_array,marker::af_marker_type,props)
    af_error(ccall((:af_draw_scatter_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X,Y,Z,marker,props))
end

export af_draw_scatter_3d

function af_draw_hist(wind::af_window,X::af_array,minval::Cdouble,maxval::Cdouble,props)
    af_error(ccall((:af_draw_hist,af_lib),af_err,(af_window,af_array,Cdouble,Cdouble,Ptr{af_cell}),wind,X,minval,maxval,props))
end

export af_draw_hist

function af_draw_surface(wind::af_window,xVals::af_array,yVals::af_array,S::af_array,props)
    af_error(ccall((:af_draw_surface,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,xVals,yVals,S,props))
end

export af_draw_surface

function af_draw_vector_field_nd(wind::af_window,points::af_array,directions::af_array,props)
    af_error(ccall((:af_draw_vector_field_nd,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,points,directions,props))
end

export af_draw_vector_field_nd

function af_draw_vector_field_3d(wind::af_window,xPoints::af_array,yPoints::af_array,zPoints::af_array,xDirs::af_array,yDirs::af_array,zDirs::af_array,props)
    af_error(ccall((:af_draw_vector_field_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints,yPoints,zPoints,xDirs,yDirs,zDirs,props))
end

export af_draw_vector_field_3d

function af_draw_vector_field_2d(wind::af_window,xPoints::af_array,yPoints::af_array,xDirs::af_array,yDirs::af_array,props)
    af_error(ccall((:af_draw_vector_field_2d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints,yPoints,xDirs,yDirs,props))
end

export af_draw_vector_field_2d

function af_grid(wind::af_window,rows::Cint,cols::Cint)
    af_error(ccall((:af_grid,af_lib),af_err,(af_window,Cint,Cint),wind,rows,cols))
end

export af_grid

function af_set_axes_limits_compute(wind::af_window,x::af_array,y::af_array,z::af_array,exact::Bool,props)
    af_error(ccall((:af_set_axes_limits_compute,af_lib),af_err,(af_window,af_array,af_array,af_array,Bool,Ptr{af_cell}),wind,x,y,z,exact,props))
end

export af_set_axes_limits_compute

function af_set_axes_limits_2d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,exact::Bool,props)
    af_error(ccall((:af_set_axes_limits_2d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,exact,props))
end

export af_set_axes_limits_2d

function af_set_axes_limits_3d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,zmin::Cfloat,zmax::Cfloat,exact::Bool,props)
    af_error(ccall((:af_set_axes_limits_3d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,zmin,zmax,exact,props))
end

export af_set_axes_limits_3d

function af_set_axes_titles(wind::af_window,xtitle,ytitle,ztitle,props)
    af_error(ccall((:af_set_axes_titles,af_lib),af_err,(af_window,Cstring,Cstring,Cstring,Ptr{af_cell}),wind,xtitle,ytitle,ztitle,props))
end

export af_set_axes_titles

function af_show(wind::af_window)
    af_error(ccall((:af_show,af_lib),af_err,(af_window,),wind))
end

export af_show

function af_is_window_closed(wind::af_window)
    out = RefValue{Bool}(0)
    af_error(ccall((:af_is_window_closed,af_lib),af_err,(Ptr{Bool},af_window),out,wind))
    out[]
end

export af_is_window_closed

function af_set_visibility(wind::af_window,is_visible::Bool)
    af_error(ccall((:af_set_visibility,af_lib),af_err,(af_window,Bool),wind,is_visible))
end

export af_set_visibility

function af_destroy_window(wind::af_window)
    af_error(ccall((:af_destroy_window,af_lib),af_err,(af_window,),wind))
end

export af_destroy_window

function af_gradient(_in::af_array)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    af_error(ccall((:af_gradient,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array),dx,dy,_in))
    (dx[],dy[])
end

export af_gradient

function af_load_image(filename,isColor::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_load_image,af_lib),af_err,(Ptr{af_array},Cstring,Bool),out,filename,isColor))
    out[]
end

export af_load_image

function af_save_image(filename,_in::af_array)
    af_error(ccall((:af_save_image,af_lib),af_err,(Cstring,af_array),filename,_in))
end

export af_save_image

function af_load_image_memory()
    out = RefValue{af_array}(0)
    ptr = RefValue{Void}(0)
    af_error(ccall((:af_load_image_memory,af_lib),af_err,(Ptr{af_array},Ptr{Void}),out,ptr))
    (out[],ptr[])
end

export af_load_image_memory

function af_save_image_memory(_in::af_array,format::af_image_format)
    ptr = RefValue{Ptr{Void}}(0)
    af_error(ccall((:af_save_image_memory,af_lib),af_err,(Ptr{Ptr{Void}},af_array,af_image_format),ptr,_in,format))
    ptr[]
end

export af_save_image_memory

function af_delete_image_memory()
    ptr = RefValue{Void}(0)
    af_error(ccall((:af_delete_image_memory,af_lib),af_err,(Ptr{Void},),ptr))
    ptr[]
end

export af_delete_image_memory

function af_load_image_native(filename)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_load_image_native,af_lib),af_err,(Ptr{af_array},Cstring),out,filename))
    out[]
end

export af_load_image_native

function af_save_image_native(filename,_in::af_array)
    af_error(ccall((:af_save_image_native,af_lib),af_err,(Cstring,af_array),filename,_in))
end

export af_save_image_native

function af_is_image_io_available()
    out = RefValue{Bool}(0)
    af_error(ccall((:af_is_image_io_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

export af_is_image_io_available

function af_resize(_in::af_array,odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_resize,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_interp_type),out,_in,odim0,odim1,method))
    out[]
end

export af_resize

function af_transform(_in::af_array,transform::af_array,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_transform,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t,dim_t,af_interp_type,Bool),out,_in,transform,odim0,odim1,method,inverse))
    out[]
end

export af_transform

function af_transform_coordinates(tf::af_array,d0::Cfloat,d1::Cfloat)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_transform_coordinates,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat),out,tf,d0,d1))
    out[]
end

export af_transform_coordinates

function af_rotate(_in::af_array,theta::Cfloat,crop::Bool,method::af_interp_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_rotate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Bool,af_interp_type),out,_in,theta,crop,method))
    out[]
end

export af_rotate

function af_translate(_in::af_array,trans0::Cfloat,trans1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_translate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in,trans0,trans1,odim0,odim1,method))
    out[]
end

export af_translate

function af_scale(_in::af_array,scale0::Cfloat,scale1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_scale,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in,scale0,scale1,odim0,odim1,method))
    out[]
end

export af_scale

function af_skew(_in::af_array,skew0::Cfloat,skew1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_skew,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type,Bool),out,_in,skew0,skew1,odim0,odim1,method,inverse))
    out[]
end

export af_skew

function af_histogram(_in::af_array,nbins::UInt32,minval::Cdouble,maxval::Cdouble)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_histogram,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Cdouble,Cdouble),out,_in,nbins,minval,maxval))
    out[]
end

export af_histogram

function af_dilate(_in::af_array,mask::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_dilate,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in,mask))
    out[]
end

export af_dilate

function af_dilate3(_in::af_array,mask::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_dilate3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in,mask))
    out[]
end

export af_dilate3

function af_erode(_in::af_array,mask::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_erode,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in,mask))
    out[]
end

export af_erode

function af_erode3(_in::af_array,mask::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_erode3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in,mask))
    out[]
end

export af_erode3

function af_bilateral(_in::af_array,spatial_sigma::Cfloat,chromatic_sigma::Cfloat,isColor::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_bilateral,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Bool),out,_in,spatial_sigma,chromatic_sigma,isColor))
    out[]
end

export af_bilateral

function af_mean_shift(_in::af_array,spatial_sigma::Cfloat,chromatic_sigma::Cfloat,iter::UInt32,is_color::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_mean_shift,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,UInt32,Bool),out,_in,spatial_sigma,chromatic_sigma,iter,is_color))
    out[]
end

export af_mean_shift

function af_minfilt(_in::af_array,wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_minfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in,wind_length,wind_width,edge_pad))
    out[]
end

export af_minfilt

function af_maxfilt(_in::af_array,wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_maxfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in,wind_length,wind_width,edge_pad))
    out[]
end

export af_maxfilt

function af_regions(_in::af_array,connectivity::af_connectivity,ty::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_regions,af_lib),af_err,(Ptr{af_array},af_array,af_connectivity,af_dtype),out,_in,connectivity,ty))
    out[]
end

export af_regions

function af_sobel_operator(img::af_array,ker_size::UInt32)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    af_error(ccall((:af_sobel_operator,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,UInt32),dx,dy,img,ker_size))
    (dx[],dy[])
end

export af_sobel_operator

function af_rgb2gray(_in::af_array,rPercent::Cfloat,gPercent::Cfloat,bPercent::Cfloat)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_rgb2gray,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in,rPercent,gPercent,bPercent))
    out[]
end

export af_rgb2gray

function af_gray2rgb(_in::af_array,rFactor::Cfloat,gFactor::Cfloat,bFactor::Cfloat)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_gray2rgb,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in,rFactor,gFactor,bFactor))
    out[]
end

export af_gray2rgb

function af_hist_equal(_in::af_array,hist::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_hist_equal,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in,hist))
    out[]
end

export af_hist_equal

function af_gaussian_kernel(rows::Cint,cols::Cint,sigma_r::Cdouble,sigma_c::Cdouble)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_gaussian_kernel,af_lib),af_err,(Ptr{af_array},Cint,Cint,Cdouble,Cdouble),out,rows,cols,sigma_r,sigma_c))
    out[]
end

export af_gaussian_kernel

function af_hsv2rgb(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_hsv2rgb,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_hsv2rgb

function af_rgb2hsv(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_rgb2hsv,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_rgb2hsv

function af_color_space(image::af_array,to::af_cspace_t,from::af_cspace_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_color_space,af_lib),af_err,(Ptr{af_array},af_array,af_cspace_t,af_cspace_t),out,image,to,from))
    out[]
end

export af_color_space

function af_unwrap(_in::af_array,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_unwrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in,wx,wy,sx,sy,px,py,is_column))
    out[]
end

export af_unwrap

function af_wrap(_in::af_array,ox::dim_t,oy::dim_t,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_wrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in,ox,oy,wx,wy,sx,sy,px,py,is_column))
    out[]
end

export af_wrap

function af_sat(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sat,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sat

function af_ycbcr2rgb(_in::af_array,standard::af_ycc_std)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ycbcr2rgb,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in,standard))
    out[]
end

export af_ycbcr2rgb

function af_rgb2ycbcr(_in::af_array,standard::af_ycc_std)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_rgb2ycbcr,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in,standard))
    out[]
end

export af_rgb2ycbcr

function af_moments(_in::af_array,moment::af_moment_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_moments,af_lib),af_err,(Ptr{af_array},af_array,af_moment_type),out,_in,moment))
    out[]
end

export af_moments

function af_moments_all(_in::af_array,moment::af_moment_type)
    out = RefValue{Cdouble}(0)
    af_error(ccall((:af_moments_all,af_lib),af_err,(Ptr{Cdouble},af_array,af_moment_type),out,_in,moment))
    out[]
end

export af_moments_all

function af_svd(_in::af_array)
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    af_error(ccall((:af_svd,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),u,s,vt,_in))
    (u[],s[],vt[])
end

export af_svd

function af_svd_inplace(_in::af_array)
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    af_error(ccall((:af_svd_inplace,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),u,s,vt,_in))
    (u[],s[],vt[])
end

export af_svd_inplace

function af_lu(_in::af_array)
    lower = RefValue{af_array}(0)
    upper = RefValue{af_array}(0)
    pivot = RefValue{af_array}(0)
    af_error(ccall((:af_lu,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),lower,upper,pivot,_in))
    (lower[],upper[],pivot[])
end

export af_lu

function af_lu_inplace(_in::af_array,is_lapack_piv::Bool)
    pivot = RefValue{af_array}(0)
    af_error(ccall((:af_lu_inplace,af_lib),af_err,(Ptr{af_array},af_array,Bool),pivot,_in,is_lapack_piv))
    pivot[]
end

export af_lu_inplace

function af_qr(_in::af_array)
    q = RefValue{af_array}(0)
    r = RefValue{af_array}(0)
    tau = RefValue{af_array}(0)
    af_error(ccall((:af_qr,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),q,r,tau,_in))
    (q[],r[],tau[])
end

export af_qr

function af_qr_inplace(_in::af_array)
    tau = RefValue{af_array}(0)
    af_error(ccall((:af_qr_inplace,af_lib),af_err,(Ptr{af_array},af_array),tau,_in))
    tau[]
end

export af_qr_inplace

function af_cholesky(_in::af_array,is_upper::Bool)
    out = RefValue{af_array}(0)
    info = RefValue{Cint}(0)
    af_error(ccall((:af_cholesky,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,Bool),out,info,_in,is_upper))
    (out[],info[])
end

export af_cholesky

function af_cholesky_inplace(_in::af_array,is_upper::Bool)
    info = RefValue{Cint}(0)
    af_error(ccall((:af_cholesky_inplace,af_lib),af_err,(Ptr{Cint},af_array,Bool),info,_in,is_upper))
    info[]
end

export af_cholesky_inplace

function af_solve(a::af_array,b::af_array,options::af_mat_prop)
    x = RefValue{af_array}(0)
    af_error(ccall((:af_solve,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop),x,a,b,options))
    x[]
end

export af_solve

function af_solve_lu(a::af_array,piv::af_array,b::af_array,options::af_mat_prop)
    x = RefValue{af_array}(0)
    af_error(ccall((:af_solve_lu,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_mat_prop),x,a,piv,b,options))
    x[]
end

export af_solve_lu

function af_inverse(_in::af_array,options::af_mat_prop)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_inverse,af_lib),af_err,(Ptr{af_array},af_array,af_mat_prop),out,_in,options))
    out[]
end

export af_inverse

function af_rank(_in::af_array,tol::Cdouble)
    rank = RefValue{UInt32}(0)
    af_error(ccall((:af_rank,af_lib),af_err,(Ptr{UInt32},af_array,Cdouble),rank,_in,tol))
    rank[]
end

export af_rank

function af_det(_in::af_array)
    det_real = RefValue{Cdouble}(0)
    det_imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_det,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),det_real,det_imag,_in))
    (det_real[],det_imag[])
end

export af_det

function af_norm(_in::af_array,_type::af_norm_type,p::Cdouble,q::Cdouble)
    out = RefValue{Cdouble}(0)
    af_error(ccall((:af_norm,af_lib),af_err,(Ptr{Cdouble},af_array,af_norm_type,Cdouble,Cdouble),out,_in,_type,p,q))
    out[]
end

export af_norm

function af_is_lapack_available()
    out = RefValue{Bool}(0)
    af_error(ccall((:af_is_lapack_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

export af_is_lapack_available

function af_create_random_engine(rtype::af_random_engine_type,seed::uintl)
    engine = RefValue{af_random_engine}(0)
    af_error(ccall((:af_create_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type,uintl),engine,rtype,seed))
    engine[]
end

export af_create_random_engine

function af_retain_random_engine(engine::af_random_engine)
    out = RefValue{af_random_engine}(0)
    af_error(ccall((:af_retain_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine),out,engine))
    out[]
end

export af_retain_random_engine

function af_random_engine_set_type(rtype::af_random_engine_type)
    engine = RefValue{af_random_engine}(0)
    af_error(ccall((:af_random_engine_set_type,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type),engine,rtype))
    engine[]
end

export af_random_engine_set_type

function af_random_engine_get_type(engine::af_random_engine)
    rtype = RefValue{af_random_engine_type}(0)
    af_error(ccall((:af_random_engine_get_type,af_lib),af_err,(Ptr{af_random_engine_type},af_random_engine),rtype,engine))
    rtype[]
end

export af_random_engine_get_type

function af_random_uniform(ndims::UInt32,dims,_type::af_dtype,engine::af_random_engine)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_random_uniform,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,ndims,dims,_type,engine))
    out[]
end

export af_random_uniform

function af_random_normal(ndims::UInt32,dims,_type::af_dtype,engine::af_random_engine)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_random_normal,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,ndims,dims,_type,engine))
    out[]
end

export af_random_normal

function af_random_engine_set_seed(seed::uintl)
    engine = RefValue{af_random_engine}(0)
    af_error(ccall((:af_random_engine_set_seed,af_lib),af_err,(Ptr{af_random_engine},uintl),engine,seed))
    engine[]
end

export af_random_engine_set_seed

function af_get_default_random_engine()
    engine = RefValue{af_random_engine}(0)
    af_error(ccall((:af_get_default_random_engine,af_lib),af_err,(Ptr{af_random_engine},),engine))
    engine[]
end

export af_get_default_random_engine

function af_set_default_random_engine_type(rtype::af_random_engine_type)
    af_error(ccall((:af_set_default_random_engine_type,af_lib),af_err,(af_random_engine_type,),rtype))
end

export af_set_default_random_engine_type

function af_random_engine_get_seed(engine::af_random_engine)
    seed = RefValue{uintl}(0)
    af_error(ccall((:af_random_engine_get_seed,af_lib),af_err,(Ptr{uintl},af_random_engine),seed,engine))
    seed[]
end

export af_random_engine_get_seed

function af_release_random_engine(engine::af_random_engine)
    af_error(ccall((:af_release_random_engine,af_lib),af_err,(af_random_engine,),engine))
end

export af_release_random_engine

function af_randu(ndims::UInt32,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_randu,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,ndims,dims,_type))
    out[]
end

export af_randu

function af_randn(ndims::UInt32,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_randn,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,ndims,dims,_type))
    out[]
end

export af_randn

function af_set_seed(seed::uintl)
    af_error(ccall((:af_set_seed,af_lib),af_err,(uintl,),seed))
end

export af_set_seed

function af_get_seed()
    seed = RefValue{uintl}(0)
    af_error(ccall((:af_get_seed,af_lib),af_err,(Ptr{uintl},),seed))
    seed[]
end

export af_get_seed

function af_approx1(_in::af_array,pos::af_array,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_approx1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_interp_type,Cfloat),out,_in,pos,method,offGrid))
    out[]
end

export af_approx1

function af_approx2(_in::af_array,pos0::af_array,pos1::af_array,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_approx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_interp_type,Cfloat),out,_in,pos0,pos1,method,offGrid))
    out[]
end

export af_approx2

function af_fft(_in::af_array,norm_factor::Cdouble,odim0::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in,norm_factor,odim0))
    out[]
end

export af_fft

function af_fft_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_fft_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_fft_inplace

function af_fft2(_in::af_array,norm_factor::Cdouble,odim0::dim_t,odim1::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in,norm_factor,odim0,odim1))
    out[]
end

export af_fft2

function af_fft2_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_fft2_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_fft2_inplace

function af_fft3(_in::af_array,norm_factor::Cdouble,odim0::dim_t,odim1::dim_t,odim2::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in,norm_factor,odim0,odim1,odim2))
    out[]
end

export af_fft3

function af_fft3_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_fft3_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_fft3_inplace

function af_ifft(_in::af_array,norm_factor::Cdouble,odim0::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ifft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in,norm_factor,odim0))
    out[]
end

export af_ifft

function af_ifft_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_ifft_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_ifft_inplace

function af_ifft2(_in::af_array,norm_factor::Cdouble,odim0::dim_t,odim1::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ifft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in,norm_factor,odim0,odim1))
    out[]
end

export af_ifft2

function af_ifft2_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_ifft2_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_ifft2_inplace

function af_ifft3(_in::af_array,norm_factor::Cdouble,odim0::dim_t,odim1::dim_t,odim2::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_ifft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in,norm_factor,odim0,odim1,odim2))
    out[]
end

export af_ifft3

function af_ifft3_inplace(_in::af_array,norm_factor::Cdouble)
    af_error(ccall((:af_ifft3_inplace,af_lib),af_err,(af_array,Cdouble),_in,norm_factor))
end

export af_ifft3_inplace

function af_fft_r2c(_in::af_array,norm_factor::Cdouble,pad0::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in,norm_factor,pad0))
    out[]
end

export af_fft_r2c

function af_fft2_r2c(_in::af_array,norm_factor::Cdouble,pad0::dim_t,pad1::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft2_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in,norm_factor,pad0,pad1))
    out[]
end

export af_fft2_r2c

function af_fft3_r2c(_in::af_array,norm_factor::Cdouble,pad0::dim_t,pad1::dim_t,pad2::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft3_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in,norm_factor,pad0,pad1,pad2))
    out[]
end

export af_fft3_r2c

function af_fft_c2r(_in::af_array,norm_factor::Cdouble,is_odd::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in,norm_factor,is_odd))
    out[]
end

export af_fft_c2r

function af_fft2_c2r(_in::af_array,norm_factor::Cdouble,is_odd::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft2_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in,norm_factor,is_odd))
    out[]
end

export af_fft2_c2r

function af_fft3_c2r(_in::af_array,norm_factor::Cdouble,is_odd::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft3_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in,norm_factor,is_odd))
    out[]
end

export af_fft3_c2r

function af_convolve1(signal::af_array,filter::af_array,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal,filter,mode,domain))
    out[]
end

export af_convolve1

function af_convolve2(signal::af_array,filter::af_array,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal,filter,mode,domain))
    out[]
end

export af_convolve2

function af_convolve3(signal::af_array,filter::af_array,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal,filter,mode,domain))
    out[]
end

export af_convolve3

function af_convolve2_sep(col_filter::af_array,row_filter::af_array,signal::af_array,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_convolve2_sep,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_conv_mode),out,col_filter,row_filter,signal,mode))
    out[]
end

export af_convolve2_sep

function af_fft_convolve1(signal::af_array,filter::af_array,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal,filter,mode))
    out[]
end

export af_fft_convolve1

function af_fft_convolve2(signal::af_array,filter::af_array,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal,filter,mode))
    out[]
end

export af_fft_convolve2

function af_fft_convolve3(signal::af_array,filter::af_array,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_fft_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal,filter,mode))
    out[]
end

export af_fft_convolve3

function af_fir(b::af_array,x::af_array)
    y = RefValue{af_array}(0)
    af_error(ccall((:af_fir,af_lib),af_err,(Ptr{af_array},af_array,af_array),y,b,x))
    y[]
end

export af_fir

function af_iir(b::af_array,a::af_array,x::af_array)
    y = RefValue{af_array}(0)
    af_error(ccall((:af_iir,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),y,b,a,x))
    y[]
end

export af_iir

function af_medfilt(_in::af_array,wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_medfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in,wind_length,wind_width,edge_pad))
    out[]
end

export af_medfilt

function af_medfilt1(_in::af_array,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_medfilt1,af_lib),af_err,(Ptr{af_array},af_array,dim_t,af_border_type),out,_in,wind_width,edge_pad))
    out[]
end

export af_medfilt1

function af_medfilt2(_in::af_array,wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_medfilt2,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in,wind_length,wind_width,edge_pad))
    out[]
end

export af_medfilt2

function af_set_fft_plan_cache_size(cache_size::Csize_t)
    af_error(ccall((:af_set_fft_plan_cache_size,af_lib),af_err,(Csize_t,),cache_size))
end

export af_set_fft_plan_cache_size

function af_create_sparse_array(nRows::dim_t,nCols::dim_t,values::af_array,rowIdx::af_array,colIdx::af_array,stype::af_storage)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_create_sparse_array,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,af_array,af_array,af_array,af_storage),out,nRows,nCols,values,rowIdx,colIdx,stype))
    out[]
end

export af_create_sparse_array

function af_create_sparse_array_from_ptr(nRows::dim_t,nCols::dim_t,nNZ::dim_t,values,rowIdx,colIdx,_type::af_dtype,stype::af_storage,src::af_source)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_create_sparse_array_from_ptr,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,dim_t,Ptr{Void},Ptr{Cint},Ptr{Cint},af_dtype,af_storage,af_source),out,nRows,nCols,nNZ,values,rowIdx,colIdx,_type,stype,src))
    out[]
end

export af_create_sparse_array_from_ptr

function af_create_sparse_array_from_dense(dense::af_array,stype::af_storage)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_create_sparse_array_from_dense,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,dense,stype))
    out[]
end

export af_create_sparse_array_from_dense

function af_sparse_convert_to(_in::af_array,destStorage::af_storage)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sparse_convert_to,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,_in,destStorage))
    out[]
end

export af_sparse_convert_to

function af_sparse_to_dense(sparse::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sparse_to_dense,af_lib),af_err,(Ptr{af_array},af_array),out,sparse))
    out[]
end

export af_sparse_to_dense

function af_sparse_get_info(_in::af_array)
    values = RefValue{af_array}(0)
    rowIdx = RefValue{af_array}(0)
    colIdx = RefValue{af_array}(0)
    stype = RefValue{af_storage}(0)
    af_error(ccall((:af_sparse_get_info,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},Ptr{af_storage},af_array),values,rowIdx,colIdx,stype,_in))
    (values[],rowIdx[],colIdx[],stype[])
end

export af_sparse_get_info

function af_sparse_get_values(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sparse_get_values,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sparse_get_values

function af_sparse_get_row_idx(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sparse_get_row_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sparse_get_row_idx

function af_sparse_get_col_idx(_in::af_array)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_sparse_get_col_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in))
    out[]
end

export af_sparse_get_col_idx

function af_sparse_get_nnz(_in::af_array)
    out = RefValue{dim_t}(0)
    af_error(ccall((:af_sparse_get_nnz,af_lib),af_err,(Ptr{dim_t},af_array),out,_in))
    out[]
end

export af_sparse_get_nnz

function af_sparse_get_storage(_in::af_array)
    out = RefValue{af_storage}(0)
    af_error(ccall((:af_sparse_get_storage,af_lib),af_err,(Ptr{af_storage},af_array),out,_in))
    out[]
end

export af_sparse_get_storage

function af_mean(_in::af_array,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_mean,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in,dim))
    out[]
end

export af_mean

function af_mean_weighted(_in::af_array,weights::af_array,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_mean_weighted,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t),out,_in,weights,dim))
    out[]
end

export af_mean_weighted

function af_var(_in::af_array,isbiased::Bool,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_var,af_lib),af_err,(Ptr{af_array},af_array,Bool,dim_t),out,_in,isbiased,dim))
    out[]
end

export af_var

function af_var_weighted(_in::af_array,weights::af_array,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_var_weighted,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t),out,_in,weights,dim))
    out[]
end

export af_var_weighted

function af_stdev(_in::af_array,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_stdev,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in,dim))
    out[]
end

export af_stdev

function af_cov(X::af_array,Y::af_array,isbiased::Bool)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_cov,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,X,Y,isbiased))
    out[]
end

export af_cov

function af_median(_in::af_array,dim::dim_t)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_median,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in,dim))
    out[]
end

export af_median

function af_mean_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_mean_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_mean_all

function af_mean_all_weighted(_in::af_array,weights::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_mean_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),real,imag,_in,weights))
    (real[],imag[])
end

export af_mean_all_weighted

function af_var_all(_in::af_array,isbiased::Bool)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    af_error(ccall((:af_var_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Bool),realVal,imagVal,_in,isbiased))
    (realVal[],imagVal[])
end

export af_var_all

function af_var_all_weighted(_in::af_array,weights::af_array)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    af_error(ccall((:af_var_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,_in,weights))
    (realVal[],imagVal[])
end

export af_var_all_weighted

function af_stdev_all(_in::af_array)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    af_error(ccall((:af_stdev_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in))
    (real[],imag[])
end

export af_stdev_all

function af_median_all(_in::af_array)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    af_error(ccall((:af_median_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),realVal,imagVal,_in))
    (realVal[],imagVal[])
end

export af_median_all

function af_corrcoef(X::af_array,Y::af_array)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    af_error(ccall((:af_corrcoef,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,X,Y))
    (realVal[],imagVal[])
end

export af_corrcoef

function af_fast(_in::af_array,thr::Cfloat,arc_length::UInt32,non_max::Bool,feature_ratio::Cfloat,edge::UInt32)
    out = RefValue{af_features}(0)
    af_error(ccall((:af_fast,af_lib),af_err,(Ptr{af_features},af_array,Cfloat,UInt32,Bool,Cfloat,UInt32),out,_in,thr,arc_length,non_max,feature_ratio,edge))
    out[]
end

export af_fast

function af_harris(_in::af_array,max_corners::UInt32,min_response::Cfloat,sigma::Cfloat,block_size::UInt32,k_thr::Cfloat)
    out = RefValue{af_features}(0)
    af_error(ccall((:af_harris,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,UInt32,Cfloat),out,_in,max_corners,min_response,sigma,block_size,k_thr))
    out[]
end

export af_harris

function af_orb(_in::af_array,fast_thr::Cfloat,max_feat::UInt32,scl_fctr::Cfloat,levels::UInt32,blur_img::Bool)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    af_error(ccall((:af_orb,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,Cfloat,UInt32,Cfloat,UInt32,Bool),feat,desc,_in,fast_thr,max_feat,scl_fctr,levels,blur_img))
    (feat[],desc[])
end

export af_orb

function af_sift(_in::af_array,n_layers::UInt32,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    af_error(ccall((:af_sift,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in,n_layers,contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],desc[])
end

export af_sift

function af_gloh(_in::af_array,n_layers::UInt32,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    af_error(ccall((:af_gloh,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in,n_layers,contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],desc[])
end

export af_gloh

function af_hamming_matcher(query::af_array,train::af_array,dist_dim::dim_t,n_dist::UInt32)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    af_error(ccall((:af_hamming_matcher,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32),idx,dist,query,train,dist_dim,n_dist))
    (idx[],dist[])
end

export af_hamming_matcher

function af_nearest_neighbour(query::af_array,train::af_array,dist_dim::dim_t,n_dist::UInt32,dist_type::af_match_type)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    af_error(ccall((:af_nearest_neighbour,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32,af_match_type),idx,dist,query,train,dist_dim,n_dist,dist_type))
    (idx[],dist[])
end

export af_nearest_neighbour

function af_match_template(search_img::af_array,template_img::af_array,m_type::af_match_type)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_match_template,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_match_type),out,search_img,template_img,m_type))
    out[]
end

export af_match_template

function af_susan(_in::af_array,radius::UInt32,diff_thr::Cfloat,geom_thr::Cfloat,feature_ratio::Cfloat,edge::UInt32)
    out = RefValue{af_features}(0)
    af_error(ccall((:af_susan,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,Cfloat,UInt32),out,_in,radius,diff_thr,geom_thr,feature_ratio,edge))
    out[]
end

export af_susan

function af_dog(_in::af_array,radius1::Cint,radius2::Cint)
    out = RefValue{af_array}(0)
    af_error(ccall((:af_dog,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint),out,_in,radius1,radius2))
    out[]
end

export af_dog

function af_homography(x_src::af_array,y_src::af_array,x_dst::af_array,y_dst::af_array,htype::af_homography_type,inlier_thr::Cfloat,iterations::UInt32,otype::af_dtype)
    H = RefValue{af_array}(0)
    inliers = RefValue{Cint}(0)
    af_error(ccall((:af_homography,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,af_array,af_array,af_array,af_homography_type,Cfloat,UInt32,af_dtype),H,inliers,x_src,y_src,x_dst,y_dst,htype,inlier_thr,iterations,otype))
    (H[],inliers[])
end

export af_homography
