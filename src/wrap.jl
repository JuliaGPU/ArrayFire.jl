# Julia wrapper for header: /usr/include/arrayfire.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function sum{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export sum

function sum_nan{T,N}(_in::AFArray{T,N},dim::Integer,nanval::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sum_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in.arr,Cint(dim),Cdouble(nanval)))
    AFArray{T,N}(out[])
end

export sum_nan

function prod{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_product,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export prod

function product_nan{T,N}(_in::AFArray{T,N},dim::Integer,nanval::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_product_nan,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cdouble),out,_in.arr,Cint(dim),Cdouble(nanval)))
    AFArray{T,N}(out[])
end

export product_nan

function minimum{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_min,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export minimum

function maximum{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_max,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export maximum

function all_true(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_all_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray!(out[])
end

export all_true

function any_true(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_any_true,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray!(out[])
end

export any_true

function count{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_count,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export count

function sum_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_sum_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export sum_all

function sum_nan_all(_in::AFArray,nanval::Real)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_sum_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in.arr,Cdouble(nanval)))
    (real[],imag[])
end

export sum_nan_all

function product_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_product_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export product_all

function product_nan_all(_in::AFArray,nanval::Real)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_product_nan_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Cdouble),real,imag,_in.arr,Cdouble(nanval)))
    (real[],imag[])
end

export product_nan_all

function min_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_min_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export min_all

function max_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_max_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export max_all

function all_true_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_all_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export all_true_all

function any_true_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_any_true_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export any_true_all

function count_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_count_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export count_all

function imin(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    _error(ccall((:af_imin,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in.arr,Cint(dim)))
    (AFArray!(out[]),AFArray!(idx[]))
end

export imin

function imax(_in::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    idx = RefValue{af_array}(0)
    _error(ccall((:af_imax,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,Cint),out,idx,_in.arr,Cint(dim)))
    (AFArray!(out[]),AFArray!(idx[]))
end

export imax

function imin_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    _error(ccall((:af_imin_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in.arr))
    (real[],imag[],idx[])
end

export imin_all

function imax_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    idx = RefValue{UInt32}(0)
    _error(ccall((:af_imax_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},Ptr{UInt32},af_array),real,imag,idx,_in.arr))
    (real[],imag[],idx[])
end

export imax_all

function accum{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_accum,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export accum

function scan{T,N}(_in::AFArray{T,N},dim::Integer,op::af_binary_op,inclusive_scan::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_scan,af_lib),af_err,(Ptr{af_array},af_array,Cint,af_binary_op,Bool),out,_in.arr,Cint(dim),op,inclusive_scan))
    AFArray{T,N}(out[])
end

export scan

function scan_by_key(key::AFArray,_in::AFArray,dim::Integer,op::af_binary_op,inclusive_scan::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_scan_by_key,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cint,af_binary_op,Bool),out,key.arr,_in.arr,Cint(dim),op,inclusive_scan))
    AFArray!(out[])
end

export scan_by_key

function where{T,N}(_in::AFArray{T,N})
    idx = RefValue{af_array}(0)
    _error(ccall((:af_where,af_lib),af_err,(Ptr{af_array},af_array),idx,_in.arr))
    AFArray{T,N}(idx[])
end

export where

function diff1{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diff1,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export diff1

function diff2{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diff2,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(dim)))
    AFArray{T,N}(out[])
end

export diff2

function sort{T,N}(_in::AFArray{T,N},dim::Integer,isAscending::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sort,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Bool),out,_in.arr,UInt32(dim),isAscending))
    AFArray{T,N}(out[])
end

export sort

function sort_index(_in::AFArray,dim::Integer,isAscending::Bool)
    out = RefValue{af_array}(0)
    indices = RefValue{af_array}(0)
    _error(ccall((:af_sort_index,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,UInt32,Bool),out,indices,_in.arr,UInt32(dim),isAscending))
    (AFArray!(out[]),AFArray!(indices[]))
end

export sort_index

function sort_by_key(keys::AFArray,values::AFArray,dim::Integer,isAscending::Bool)
    out_keys = RefValue{af_array}(0)
    out_values = RefValue{af_array}(0)
    _error(ccall((:af_sort_by_key,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,UInt32,Bool),out_keys,out_values,keys.arr,values.arr,UInt32(dim),isAscending))
    (AFArray!(out_keys[]),AFArray!(out_values[]))
end

export sort_by_key

function set_unique{T,N}(_in::AFArray{T,N},is_sorted::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_unique,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_sorted))
    AFArray{T,N}(out[])
end

export set_unique

function set_union(first::AFArray,second::AFArray,is_unique::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_union,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first.arr,second.arr,is_unique))
    AFArray!(out[])
end

export set_union

function set_intersect(first::AFArray,second::AFArray,is_unique::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_set_intersect,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,first.arr,second.arr,is_unique))
    AFArray!(out[])
end

export set_intersect

function add(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_add,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export add

function sub(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sub,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export sub

function mul(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mul,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export mul

function div(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_div,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export div

function lt(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_lt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export lt

function gt(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_gt,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export gt

function le(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_le,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export le

function ge(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_ge,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export ge

function eq(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_eq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export eq

function neq(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_neq,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export neq

function and(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_and,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export and

function or(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_or,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export or

function not{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_not,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export not

function bitand(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitand,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export bitand

function bitor(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export bitor

function bitxor(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitxor,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export bitxor

function bitshiftl(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitshiftl,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export bitshiftl

function bitshiftr(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bitshiftr,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export bitshiftr

function cast{T,N}(_in::AFArray{T,N},_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_cast,af_lib),af_err,(Ptr{af_array},af_array,af_dtype),out,_in.arr,_type))
    AFArray{T,N}(out[])
end

export cast

function minof(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_minof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export minof

function maxof(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_maxof,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export maxof

function clamp(_in::AFArray,lo::AFArray,hi::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_clamp,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,Bool),out,_in.arr,lo.arr,hi.arr,batch))
    AFArray!(out[])
end

export clamp

function rem(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_rem,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export rem

function mod(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mod,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export mod

function abs{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_abs,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export abs

function arg{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_arg,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export arg

function signbit{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sign,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export signbit

function round{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_round,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export round

function trunc{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_trunc,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export trunc

function floor{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_floor,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export floor

function ceil{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_ceil,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export ceil

function hypot(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_hypot,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export hypot

function sin{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sin,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sin

function cos{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cos,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export cos

function tan{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_tan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export tan

function asin{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_asin,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export asin

function acos{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_acos,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export acos

function atan{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_atan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export atan

function atan2(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_atan2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export atan2

function cplx2(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_cplx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export cplx2

function cplx{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cplx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export cplx

function real{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_real,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export real

function imag{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_imag,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export imag

function conjg{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_conjg,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export conjg

function sinh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sinh

function cosh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export cosh

function tanh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_tanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export tanh

function asinh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_asinh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export asinh

function acosh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_acosh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export acosh

function atanh{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_atanh,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export atanh

function root(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_root,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export root

function pow(lhs::AFArray,rhs::AFArray,batch::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_pow,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,lhs.arr,rhs.arr,batch))
    AFArray!(out[])
end

export pow

function pow2{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_pow2,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export pow2

function exp{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_exp,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export exp

function sigmoid{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sigmoid,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sigmoid

function expm1{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_expm1,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export expm1

function erf{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_erf,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export erf

function erfc{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_erfc,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export erfc

function log{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_log,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export log

function log1p{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_log1p,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export log1p

function log10{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_log10,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export log10

function log2{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_log2,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export log2

function sqrt{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sqrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sqrt

function cbrt{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_cbrt,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export cbrt

function factorial{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_factorial,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export factorial

function tgamma{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_tgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export tgamma

function lgamma{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_lgamma,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export lgamma

function iszero{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_iszero,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export iszero

function isinf{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_isinf,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export isinf

function isnan{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_isnan,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export isnan

function make_seq(_begin::Real,_end::Real,step::Real)
    ccall((:af_make_seq,af_lib),af_seq,(Cdouble,Cdouble,Cdouble),Cdouble(_begin),Cdouble(_end),Cdouble(step))
end

export make_seq

function print_array(arr::AFArray)
    _error(ccall((:af_print_array,af_lib),af_err,(af_array,),arr.arr))
end

export print_array

function print_array_gen(exp,arr::AFArray,precision::Integer)
    _error(ccall((:af_print_array_gen,af_lib),af_err,(Cstring,af_array,Cint),exp,arr.arr,Cint(precision)))
end

export print_array_gen

function save_array(key,arr::AFArray,filename,append::Bool)
    index = RefValue{Cint}(0)
    _error(ccall((:af_save_array,af_lib),af_err,(Ptr{Cint},Cstring,af_array,Cstring,Bool),index,key,arr.arr,filename,append))
    index[]
end

export save_array

function read_array_index(filename,index::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_read_array_index,af_lib),af_err,(Ptr{af_array},Cstring,UInt32),out,filename,UInt32(index)))
    AFArray!(out[])
end

export read_array_index

function read_array_key(filename,key)
    out = RefValue{af_array}(0)
    _error(ccall((:af_read_array_key,af_lib),af_err,(Ptr{af_array},Cstring,Cstring),out,filename,key))
    AFArray!(out[])
end

export read_array_key

function read_array_key_check(filename,key)
    index = RefValue{Cint}(0)
    _error(ccall((:af_read_array_key_check,af_lib),af_err,(Ptr{Cint},Cstring,Cstring),index,filename,key))
    index[]
end

export read_array_key_check

function array_to_string(output,exp,arr::AFArray,precision::Integer,transpose::Bool)
    _error(ccall((:af_array_to_string,af_lib),af_err,(Ptr{Cstring},Cstring,af_array,Cint,Bool),output,exp,arr.arr,Cint(precision),transpose))
end

export array_to_string

function afversion()
    major = RefValue{Cint}(0)
    minor = RefValue{Cint}(0)
    patch = RefValue{Cint}(0)
    _error(ccall((:af_get_version,af_lib),af_err,(Ptr{Cint},Ptr{Cint},Ptr{Cint}),major,minor,patch))
    (major[],minor[],patch[])
end

export afversion

function get_revision()
    ccall((:af_get_revision,af_lib),Cstring,())
end

export get_revision

function get_size_of(_type::af_dtype)
    size = RefValue{Csize_t}(0)
    _error(ccall((:af_get_size_of,af_lib),af_err,(Ptr{Csize_t},af_dtype),size,_type))
    size[]
end

export get_size_of

function index{T,N}(_in::AFArray{T,N},ndims::Integer,index)
    out = RefValue{af_array}(0)
    _error(ccall((:af_index,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq}),out,_in.arr,UInt32(ndims),index))
    AFArray{T,N}(out[])
end

export index

function lookup(_in::AFArray,indices::AFArray,dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_lookup,af_lib),af_err,(Ptr{af_array},af_array,af_array,UInt32),out,_in.arr,indices.arr,UInt32(dim)))
    AFArray!(out[])
end

export lookup

function assign_seq(lhs::AFArray,ndims::Integer,indices,rhs::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_assign_seq,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{af_seq},af_array),out,lhs.arr,UInt32(ndims),indices,rhs.arr))
    AFArray!(out[])
end

export assign_seq

function index_gen{T,N}(_in::AFArray{T,N},ndims::dim_t,indices)
    out = RefValue{af_array}(0)
    _error(ccall((:af_index_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t}),out,_in.arr,ndims,indices))
    AFArray{T,N}(out[])
end

export index_gen

function assign_gen(lhs::AFArray,ndims::dim_t,indices,rhs::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_assign_gen,af_lib),af_err,(Ptr{af_array},af_array,dim_t,Ptr{af_index_t},af_array),out,lhs.arr,ndims,indices,rhs.arr))
    AFArray!(out[])
end

export assign_gen

function create_indexers()
    indexers = RefValue{Ptr{af_index_t}}(0)
    _error(ccall((:af_create_indexers,af_lib),af_err,(Ptr{Ptr{af_index_t}},),indexers))
    indexers[]
end

export create_indexers

function set_array_indexer(idx::AFArray,dim::dim_t)
    indexer = RefValue{af_index_t}(0)
    _error(ccall((:af_set_array_indexer,af_lib),af_err,(Ptr{af_index_t},af_array,dim_t),indexer,idx.arr,dim))
    indexer[]
end

export set_array_indexer

function set_seq_indexer(dim::dim_t,is_batch::Bool)
    indexer = RefValue{af_index_t}(0)
    idx = RefValue{af_seq}(0)
    _error(ccall((:af_set_seq_indexer,af_lib),af_err,(Ptr{af_index_t},Ptr{af_seq},dim_t,Bool),indexer,idx,dim,is_batch))
    (indexer[],idx[])
end

export set_seq_indexer

function set_seq_param_indexer(_begin::Real,_end::Real,step::Real,dim::dim_t,is_batch::Bool)
    indexer = RefValue{af_index_t}(0)
    _error(ccall((:af_set_seq_param_indexer,af_lib),af_err,(Ptr{af_index_t},Cdouble,Cdouble,Cdouble,dim_t,Bool),indexer,Cdouble(_begin),Cdouble(_end),Cdouble(step),dim,is_batch))
    indexer[]
end

export set_seq_param_indexer

function release_indexers()
    indexers = RefValue{af_index_t}(0)
    _error(ccall((:af_release_indexers,af_lib),af_err,(Ptr{af_index_t},),indexers))
    indexers[]
end

export release_indexers

function create_handle(ndims::Integer,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_create_handle,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),arr,UInt32(ndims),dims,_type))
    AFArray!(arr[])
end

export create_handle

function copy{T,N}(_in::AFArray{T,N})
    arr = RefValue{af_array}(0)
    _error(ccall((:af_copy_array,af_lib),af_err,(Ptr{af_array},af_array),arr,_in.arr))
    AFArray{T,N}(arr[])
end

export copy

function write_array(arr::AFArray,data,bytes::Csize_t,src::af_source)
    _error(ccall((:af_write_array,af_lib),af_err,(af_array,Ptr{Void},Csize_t,af_source),arr.arr,data,bytes,src))
end

export write_array

function get_data_ptr(data,arr::AFArray)
    _error(ccall((:af_get_data_ptr,af_lib),af_err,(Ptr{Void},af_array),data,arr.arr))
end

export get_data_ptr

function release_array(arr::AFArray)
    _error(ccall((:af_release_array,af_lib),af_err,(af_array,),arr.arr))
end

export release_array

function afeval(_in::AFArray)
    _error(ccall((:af_eval,af_lib),af_err,(af_array,),_in.arr))
end

export afeval

function eval_multiple(num::Integer,arrays)
    _error(ccall((:af_eval_multiple,af_lib),af_err,(Cint,Ptr{af_array}),Cint(num),arrays))
end

export eval_multiple

function set_manual_eval_flag(flag::Bool)
    _error(ccall((:af_set_manual_eval_flag,af_lib),af_err,(Bool,),flag))
end

export set_manual_eval_flag

function get_manual_eval_flag()
    flag = RefValue{Bool}(0)
    _error(ccall((:af_get_manual_eval_flag,af_lib),af_err,(Ptr{Bool},),flag))
    flag[]
end

export get_manual_eval_flag

function get_elements(arr::AFArray)
    elems = RefValue{dim_t}(0)
    _error(ccall((:af_get_elements,af_lib),af_err,(Ptr{dim_t},af_array),elems,arr.arr))
    elems[]
end

export get_elements

function get_dims(arr::AFArray)
    d0 = RefValue{dim_t}(0)
    d1 = RefValue{dim_t}(0)
    d2 = RefValue{dim_t}(0)
    d3 = RefValue{dim_t}(0)
    _error(ccall((:af_get_dims,af_lib),af_err,(Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},Ptr{dim_t},af_array),d0,d1,d2,d3,arr.arr))
    (d0[],d1[],d2[],d3[])
end

export get_dims

function is_empty(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_empty,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_empty

function is_scalar(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_scalar,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_scalar

function is_row(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_row,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_row

function is_column(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_column,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_column

function is_vector(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_vector,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_vector

function is_complex(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_complex,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_complex

function is_real(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_real,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_real

function is_double(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_double,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_double

function is_single(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_single,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_single

function is_realfloating(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_realfloating,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_realfloating

function is_floating(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_floating,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_floating

function is_integer(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_integer,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_integer

function is_bool(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_bool,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_bool

function is_sparse(arr::AFArray)
    result = RefValue{Bool}(0)
    _error(ccall((:af_is_sparse,af_lib),af_err,(Ptr{Bool},af_array),result,arr.arr))
    result[]
end

export is_sparse

function set_backend(bknd::af_backend)
    _error(ccall((:af_set_backend,af_lib),af_err,(af_backend,),bknd))
end

export set_backend

function get_backend_count()
    num_backends = RefValue{UInt32}(0)
    _error(ccall((:af_get_backend_count,af_lib),af_err,(Ptr{UInt32},),num_backends))
    num_backends[]
end

export get_backend_count

function get_available_backends()
    backends = RefValue{Cint}(0)
    _error(ccall((:af_get_available_backends,af_lib),af_err,(Ptr{Cint},),backends))
    backends[]
end

export get_available_backends

function get_backend_id(_in::AFArray)
    backend = RefValue{af_backend}(0)
    _error(ccall((:af_get_backend_id,af_lib),af_err,(Ptr{af_backend},af_array),backend,_in.arr))
    backend[]
end

export get_backend_id

function get_active_backend()
    backend = RefValue{af_backend}(0)
    _error(ccall((:af_get_active_backend,af_lib),af_err,(Ptr{af_backend},),backend))
    backend[]
end

export get_active_backend

function get_device_id(_in::AFArray)
    device = RefValue{Cint}(0)
    _error(ccall((:af_get_device_id,af_lib),af_err,(Ptr{Cint},af_array),device,_in.arr))
    device[]
end

export get_device_id

function matmul(lhs::AFArray,rhs::AFArray,optLhs::af_mat_prop,optRhs::af_mat_prop)
    out = RefValue{af_array}(0)
    _error(ccall((:af_matmul,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs.arr,rhs.arr,optLhs,optRhs))
    AFArray!(out[])
end

export matmul

function dot(lhs::AFArray,rhs::AFArray,optLhs::af_mat_prop,optRhs::af_mat_prop)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dot,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop,af_mat_prop),out,lhs.arr,rhs.arr,optLhs,optRhs))
    AFArray!(out[])
end

export dot

function dot_all(lhs::AFArray,rhs::AFArray,optLhs::af_mat_prop,optRhs::af_mat_prop)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_dot_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array,af_mat_prop,af_mat_prop),real,imag,lhs.arr,rhs.arr,optLhs,optRhs))
    (real[],imag[])
end

export dot_all

function transpose{T,N}(_in::AFArray{T,N},conjugate::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_transpose,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,conjugate))
    AFArray{T,N}(out[])
end

export transpose

function transpose_inplace(_in::AFArray,conjugate::Bool)
    _error(ccall((:af_transpose_inplace,af_lib),af_err,(af_array,Bool),_in.arr,conjugate))
end

export transpose_inplace

function constant(val::Real,ndims::Integer,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant,af_lib),af_err,(Ptr{af_array},Cdouble,UInt32,Ptr{dim_t},af_dtype),arr,Cdouble(val),UInt32(ndims),dims,_type))
    AFArray!(arr[])
end

export constant

function constant_complex(real::Real,imag::Real,ndims::Integer,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_complex,af_lib),af_err,(Ptr{af_array},Cdouble,Cdouble,UInt32,Ptr{dim_t},af_dtype),arr,Cdouble(real),Cdouble(imag),UInt32(ndims),dims,_type))
    AFArray!(arr[])
end

export constant_complex

function constant_long(val::intl,ndims::Integer,dims)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_long,af_lib),af_err,(Ptr{af_array},intl,UInt32,Ptr{dim_t}),arr,val,UInt32(ndims),dims))
    AFArray!(arr[])
end

export constant_long

function constant_ulong(val::uintl,ndims::Integer,dims)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_constant_ulong,af_lib),af_err,(Ptr{af_array},uintl,UInt32,Ptr{dim_t}),arr,val,UInt32(ndims),dims))
    AFArray!(arr[])
end

export constant_ulong

function range(ndims::Integer,dims,seq_dim::Integer,_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_range,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},Cint,af_dtype),out,UInt32(ndims),dims,Cint(seq_dim),_type))
    AFArray!(out[])
end

export range

function iota(ndims::Integer,dims,t_ndims::Integer,tdims,_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_iota,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,UInt32(t_ndims),tdims,_type))
    AFArray!(out[])
end

export iota

function identity(ndims::Integer,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_identity,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,_type))
    AFArray!(out[])
end

export identity

function diag_create{T,N}(_in::AFArray{T,N},num::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diag_create,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(num)))
    AFArray{T,N}(out[])
end

export diag_create

function diag_extract{T,N}(_in::AFArray{T,N},num::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_diag_extract,af_lib),af_err,(Ptr{af_array},af_array,Cint),out,_in.arr,Cint(num)))
    AFArray{T,N}(out[])
end

export diag_extract

function join(dim::Integer,first::AFArray,second::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_join,af_lib),af_err,(Ptr{af_array},Cint,af_array,af_array),out,Cint(dim),first.arr,second.arr))
    AFArray!(out[])
end

export join

function join_many(dim::Integer,n_arrays::Integer,inputs)
    out = RefValue{af_array}(0)
    _error(ccall((:af_join_many,af_lib),af_err,(Ptr{af_array},Cint,UInt32,Ptr{af_array}),out,Cint(dim),UInt32(n_arrays),inputs))
    AFArray!(out[])
end

export join_many

function tile{T,N}(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_tile,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in.arr,UInt32(x),UInt32(y),UInt32(z),UInt32(w)))
    AFArray{T,N}(out[])
end

export tile

function reorder{T,N}(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_reorder,af_lib),af_err,(Ptr{af_array},af_array,UInt32,UInt32,UInt32,UInt32),out,_in.arr,UInt32(x),UInt32(y),UInt32(z),UInt32(w)))
    AFArray{T,N}(out[])
end

export reorder

function shift{T,N}(_in::AFArray{T,N},x::Integer,y::Integer,z::Integer,w::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_shift,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint,Cint,Cint),out,_in.arr,Cint(x),Cint(y),Cint(z),Cint(w)))
    AFArray{T,N}(out[])
end

export shift

function moddims{T,N}(_in::AFArray{T,N},ndims::Integer,dims)
    out = RefValue{af_array}(0)
    _error(ccall((:af_moddims,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Ptr{dim_t}),out,_in.arr,UInt32(ndims),dims))
    AFArray{T,N}(out[])
end

export moddims

function flat{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_flat,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export flat

function flip{T,N}(_in::AFArray{T,N},dim::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_flip,af_lib),af_err,(Ptr{af_array},af_array,UInt32),out,_in.arr,UInt32(dim)))
    AFArray{T,N}(out[])
end

export flip

function lower{T,N}(_in::AFArray{T,N},is_unit_diag::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_lower,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_unit_diag))
    AFArray{T,N}(out[])
end

export lower

function upper{T,N}(_in::AFArray{T,N},is_unit_diag::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_upper,af_lib),af_err,(Ptr{af_array},af_array,Bool),out,_in.arr,is_unit_diag))
    AFArray{T,N}(out[])
end

export upper

function select(cond::AFArray,a::AFArray,b::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_select,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),out,cond.arr,a.arr,b.arr))
    AFArray!(out[])
end

export select

function select_scalar_r(cond::AFArray,a::AFArray,b::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_r,af_lib),af_err,(Ptr{af_array},af_array,af_array,Cdouble),out,cond.arr,a.arr,Cdouble(b)))
    AFArray!(out[])
end

export select_scalar_r

function select_scalar_l(cond::AFArray,a::Real,b::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_select_scalar_l,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,af_array),out,cond.arr,Cdouble(a),b.arr))
    AFArray!(out[])
end

export select_scalar_l

function replace(a::AFArray,cond::AFArray,b::AFArray)
    _error(ccall((:af_replace,af_lib),af_err,(af_array,af_array,af_array),a.arr,cond.arr,b.arr))
end

export replace

function replace_scalar(a::AFArray,cond::AFArray,b::Real)
    _error(ccall((:af_replace_scalar,af_lib),af_err,(af_array,af_array,Cdouble),a.arr,cond.arr,Cdouble(b)))
end

export replace_scalar

function afinfo()
    _error(ccall((:af_info,af_lib),af_err,()))
end

export afinfo

function afinit()
    _error(ccall((:af_init,af_lib),af_err,()))
end

export afinit

function get_device_count()
    num_of_devices = RefValue{Cint}(0)
    _error(ccall((:af_get_device_count,af_lib),af_err,(Ptr{Cint},),num_of_devices))
    num_of_devices[]
end

export get_device_count

function get_dbl_support(device::Integer)
    available = RefValue{Bool}(0)
    _error(ccall((:af_get_dbl_support,af_lib),af_err,(Ptr{Bool},Cint),available,Cint(device)))
    available[]
end

export get_dbl_support

function set_device(device::Integer)
    _error(ccall((:af_set_device,af_lib),af_err,(Cint,),Cint(device)))
end

export set_device

function get_device()
    device = RefValue{Cint}(0)
    _error(ccall((:af_get_device,af_lib),af_err,(Ptr{Cint},),device))
    device[]
end

export get_device

function sync(device::Integer)
    _error(ccall((:af_sync,af_lib),af_err,(Cint,),Cint(device)))
end

export sync

function alloc_device(bytes::dim_t)
    ptr = RefValue{Ptr{Void}}(0)
    _error(ccall((:af_alloc_device,af_lib),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes))
    ptr[]
end

export alloc_device

function free_device(ptr)
    _error(ccall((:af_free_device,af_lib),af_err,(Ptr{Void},),ptr))
end

export free_device

function device_array(data,ndims::Integer,dims,_type::af_dtype)
    arr = RefValue{af_array}(0)
    _error(ccall((:af_device_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),arr,data,UInt32(ndims),dims,_type))
    AFArray!(arr[])
end

export device_array

function device_mem_info()
    alloc_bytes = RefValue{Csize_t}(0)
    alloc_buffers = RefValue{Csize_t}(0)
    lock_bytes = RefValue{Csize_t}(0)
    lock_buffers = RefValue{Csize_t}(0)
    _error(ccall((:af_device_mem_info,af_lib),af_err,(Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t}),alloc_bytes,alloc_buffers,lock_bytes,lock_buffers))
    (alloc_bytes[],alloc_buffers[],lock_bytes[],lock_buffers[])
end

export device_mem_info

function print_mem_info(msg,device_id::Integer)
    _error(ccall((:af_print_mem_info,af_lib),af_err,(Cstring,Cint),msg,Cint(device_id)))
end

export print_mem_info

function device_gc()
    _error(ccall((:af_device_gc,af_lib),af_err,()))
end

export device_gc

function set_mem_step_size(step_bytes::Csize_t)
    _error(ccall((:af_set_mem_step_size,af_lib),af_err,(Csize_t,),step_bytes))
end

export set_mem_step_size

function get_mem_step_size()
    step_bytes = RefValue{Csize_t}(0)
    _error(ccall((:af_get_mem_step_size,af_lib),af_err,(Ptr{Csize_t},),step_bytes))
    step_bytes[]
end

export get_mem_step_size

function lock_device_ptr(arr::AFArray)
    _error(ccall((:af_lock_device_ptr,af_lib),af_err,(af_array,),arr.arr))
end

export lock_device_ptr

function unlock_device_ptr(arr::AFArray)
    _error(ccall((:af_unlock_device_ptr,af_lib),af_err,(af_array,),arr.arr))
end

export unlock_device_ptr

function lock_array(arr::AFArray)
    _error(ccall((:af_lock_array,af_lib),af_err,(af_array,),arr.arr))
end

export lock_array

function unlock_array(arr::AFArray)
    _error(ccall((:af_unlock_array,af_lib),af_err,(af_array,),arr.arr))
end

export unlock_array

function is_locked_array(arr::AFArray)
    res = RefValue{Bool}(0)
    _error(ccall((:af_is_locked_array,af_lib),af_err,(Ptr{Bool},af_array),res,arr.arr))
    res[]
end

export is_locked_array

function get_device_ptr(arr::AFArray)
    ptr = RefValue{Ptr{Void}}(0)
    _error(ccall((:af_get_device_ptr,af_lib),af_err,(Ptr{Ptr{Void}},af_array),ptr,arr.arr))
    ptr[]
end

export get_device_ptr

function err_to_string(err::af_err)
    ccall((:af_err_to_string,af_lib),Cstring,(af_err,),err)
end

export err_to_string

function create_features(num::dim_t)
    feat = RefValue{af_features}(0)
    _error(ccall((:af_create_features,af_lib),af_err,(Ptr{af_features},dim_t),feat,num))
    feat[]
end

export create_features

function retain_features(feat::af_features)
    out = RefValue{af_features}(0)
    _error(ccall((:af_retain_features,af_lib),af_err,(Ptr{af_features},af_features),out,feat))
    out[]
end

export retain_features

function get_features_num(feat::af_features)
    num = RefValue{dim_t}(0)
    _error(ccall((:af_get_features_num,af_lib),af_err,(Ptr{dim_t},af_features),num,feat))
    num[]
end

export get_features_num

function get_features_xpos(feat::af_features)
    out = RefValue{af_array}(0)
    _error(ccall((:af_get_features_xpos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    AFArray!(out[])
end

export get_features_xpos

function get_features_ypos(feat::af_features)
    out = RefValue{af_array}(0)
    _error(ccall((:af_get_features_ypos,af_lib),af_err,(Ptr{af_array},af_features),out,feat))
    AFArray!(out[])
end

export get_features_ypos

function get_features_score(feat::af_features)
    score = RefValue{af_array}(0)
    _error(ccall((:af_get_features_score,af_lib),af_err,(Ptr{af_array},af_features),score,feat))
    AFArray!(score[])
end

export get_features_score

function get_features_orientation(feat::af_features)
    orientation = RefValue{af_array}(0)
    _error(ccall((:af_get_features_orientation,af_lib),af_err,(Ptr{af_array},af_features),orientation,feat))
    AFArray!(orientation[])
end

export get_features_orientation

function get_features_size(feat::af_features)
    size = RefValue{af_array}(0)
    _error(ccall((:af_get_features_size,af_lib),af_err,(Ptr{af_array},af_features),size,feat))
    AFArray!(size[])
end

export get_features_size

function release_features(feat::af_features)
    _error(ccall((:af_release_features,af_lib),af_err,(af_features,),feat))
end

export release_features

function create_window(width::Integer,height::Integer,title)
    out = RefValue{af_window}(0)
    _error(ccall((:af_create_window,af_lib),af_err,(Ptr{af_window},Cint,Cint,Cstring),out,Cint(width),Cint(height),title))
    out[]
end

export create_window

function set_position(wind::af_window,x::Integer,y::Integer)
    _error(ccall((:af_set_position,af_lib),af_err,(af_window,UInt32,UInt32),wind,UInt32(x),UInt32(y)))
end

export set_position

function set_title(wind::af_window,title)
    _error(ccall((:af_set_title,af_lib),af_err,(af_window,Cstring),wind,title))
end

export set_title

function set_size(wind::af_window,w::Integer,h::Integer)
    _error(ccall((:af_set_size,af_lib),af_err,(af_window,UInt32,UInt32),wind,UInt32(w),UInt32(h)))
end

export set_size

function draw_image(wind::af_window,_in::AFArray,props)
    _error(ccall((:af_draw_image,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,_in.arr,props))
end

export draw_image

function draw_plot(wind::af_window,X::AFArray,Y::AFArray,props)
    _error(ccall((:af_draw_plot,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,props))
end

export draw_plot

function draw_plot3(wind::af_window,P::AFArray,props)
    _error(ccall((:af_draw_plot3,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P.arr,props))
end

export draw_plot3

function draw_plot_nd(wind::af_window,P::AFArray,props)
    _error(ccall((:af_draw_plot_nd,af_lib),af_err,(af_window,af_array,Ptr{af_cell}),wind,P.arr,props))
end

export draw_plot_nd

function draw_plot_2d(wind::af_window,X::AFArray,Y::AFArray,props)
    _error(ccall((:af_draw_plot_2d,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,props))
end

export draw_plot_2d

function draw_plot_3d(wind::af_window,X::AFArray,Y::AFArray,Z::AFArray,props)
    _error(ccall((:af_draw_plot_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,X.arr,Y.arr,Z.arr,props))
end

export draw_plot_3d

function draw_scatter(wind::af_window,X::AFArray,Y::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,marker,props))
end

export draw_scatter

function draw_scatter3(wind::af_window,P::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter3,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P.arr,marker,props))
end

export draw_scatter3

function draw_scatter_nd(wind::af_window,P::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_nd,af_lib),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P.arr,marker,props))
end

export draw_scatter_nd

function draw_scatter_2d(wind::af_window,X::AFArray,Y::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_2d,af_lib),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,marker,props))
end

export draw_scatter_2d

function draw_scatter_3d(wind::af_window,X::AFArray,Y::AFArray,Z::AFArray,marker::af_marker_type,props)
    _error(ccall((:af_draw_scatter_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X.arr,Y.arr,Z.arr,marker,props))
end

export draw_scatter_3d

function draw_hist(wind::af_window,X::AFArray,minval::Real,maxval::Real,props)
    _error(ccall((:af_draw_hist,af_lib),af_err,(af_window,af_array,Cdouble,Cdouble,Ptr{af_cell}),wind,X.arr,Cdouble(minval),Cdouble(maxval),props))
end

export draw_hist

function draw_surface(wind::af_window,xVals::AFArray,yVals::AFArray,S::AFArray,props)
    _error(ccall((:af_draw_surface,af_lib),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,xVals.arr,yVals.arr,S.arr,props))
end

export draw_surface

function draw_vector_field_nd(wind::af_window,points::AFArray,directions::AFArray,props)
    _error(ccall((:af_draw_vector_field_nd,af_lib),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,points.arr,directions.arr,props))
end

export draw_vector_field_nd

function draw_vector_field_3d(wind::af_window,xPoints::AFArray,yPoints::AFArray,zPoints::AFArray,xDirs::AFArray,yDirs::AFArray,zDirs::AFArray,props)
    _error(ccall((:af_draw_vector_field_3d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints.arr,yPoints.arr,zPoints.arr,xDirs.arr,yDirs.arr,zDirs.arr,props))
end

export draw_vector_field_3d

function draw_vector_field_2d(wind::af_window,xPoints::AFArray,yPoints::AFArray,xDirs::AFArray,yDirs::AFArray,props)
    _error(ccall((:af_draw_vector_field_2d,af_lib),af_err,(af_window,af_array,af_array,af_array,af_array,Ptr{af_cell}),wind,xPoints.arr,yPoints.arr,xDirs.arr,yDirs.arr,props))
end

export draw_vector_field_2d

function grid(wind::af_window,rows::Integer,cols::Integer)
    _error(ccall((:af_grid,af_lib),af_err,(af_window,Cint,Cint),wind,Cint(rows),Cint(cols)))
end

export grid

function set_axes_limits_compute(wind::af_window,x::AFArray,y::AFArray,z::AFArray,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_compute,af_lib),af_err,(af_window,af_array,af_array,af_array,Bool,Ptr{af_cell}),wind,x.arr,y.arr,z.arr,exact,props))
end

export set_axes_limits_compute

function set_axes_limits_2d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_2d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,exact,props))
end

export set_axes_limits_2d

function set_axes_limits_3d(wind::af_window,xmin::Cfloat,xmax::Cfloat,ymin::Cfloat,ymax::Cfloat,zmin::Cfloat,zmax::Cfloat,exact::Bool,props)
    _error(ccall((:af_set_axes_limits_3d,af_lib),af_err,(af_window,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Cfloat,Bool,Ptr{af_cell}),wind,xmin,xmax,ymin,ymax,zmin,zmax,exact,props))
end

export set_axes_limits_3d

function set_axes_titles(wind::af_window,xtitle,ytitle,ztitle,props)
    _error(ccall((:af_set_axes_titles,af_lib),af_err,(af_window,Cstring,Cstring,Cstring,Ptr{af_cell}),wind,xtitle,ytitle,ztitle,props))
end

export set_axes_titles

function show(wind::af_window)
    _error(ccall((:af_show,af_lib),af_err,(af_window,),wind))
end

export show

function is_window_closed(wind::af_window)
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_window_closed,af_lib),af_err,(Ptr{Bool},af_window),out,wind))
    out[]
end

export is_window_closed

function set_visibility(wind::af_window,is_visible::Bool)
    _error(ccall((:af_set_visibility,af_lib),af_err,(af_window,Bool),wind,is_visible))
end

export set_visibility

function destroy_window(wind::af_window)
    _error(ccall((:af_destroy_window,af_lib),af_err,(af_window,),wind))
end

export destroy_window

function gradient(_in::AFArray)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    _error(ccall((:af_gradient,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array),dx,dy,_in.arr))
    (AFArray!(dx[]),AFArray!(dy[]))
end

export gradient

function load_image(filename,isColor::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image,af_lib),af_err,(Ptr{af_array},Cstring,Bool),out,filename,isColor))
    AFArray!(out[])
end

export load_image

function save_image(filename,_in::AFArray)
    _error(ccall((:af_save_image,af_lib),af_err,(Cstring,af_array),filename,_in.arr))
end

export save_image

function load_image_memory(ptr)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image_memory,af_lib),af_err,(Ptr{af_array},Ptr{Void}),out,ptr))
    AFArray!(out[])
end

export load_image_memory

function save_image_memory(_in::AFArray,format::af_image_format)
    ptr = RefValue{Ptr{Void}}(0)
    _error(ccall((:af_save_image_memory,af_lib),af_err,(Ptr{Ptr{Void}},af_array,af_image_format),ptr,_in.arr,format))
    ptr[]
end

export save_image_memory

function delete_image_memory(ptr)
    _error(ccall((:af_delete_image_memory,af_lib),af_err,(Ptr{Void},),ptr))
end

export delete_image_memory

function load_image_native(filename)
    out = RefValue{af_array}(0)
    _error(ccall((:af_load_image_native,af_lib),af_err,(Ptr{af_array},Cstring),out,filename))
    AFArray!(out[])
end

export load_image_native

function save_image_native(filename,_in::AFArray)
    _error(ccall((:af_save_image_native,af_lib),af_err,(Cstring,af_array),filename,_in.arr))
end

export save_image_native

function is_image_io_available()
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_image_io_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

export is_image_io_available

function resize{T,N}(_in::AFArray{T,N},odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_resize,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_interp_type),out,_in.arr,odim0,odim1,method))
    AFArray{T,N}(out[])
end

export resize

function transform(_in::AFArray,transform::AFArray,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_transform,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t,dim_t,af_interp_type,Bool),out,_in.arr,transform.arr,odim0,odim1,method,inverse))
    AFArray!(out[])
end

export transform

function transform_coordinates{T,N}(tf::AFArray{T,N},d0::Cfloat,d1::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_transform_coordinates,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat),out,tf.arr,d0,d1))
    AFArray{T,N}(out[])
end

export transform_coordinates

function rotate{T,N}(_in::AFArray{T,N},theta::Cfloat,crop::Bool,method::af_interp_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_rotate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Bool,af_interp_type),out,_in.arr,theta,crop,method))
    AFArray{T,N}(out[])
end

export rotate

function translate{T,N}(_in::AFArray{T,N},trans0::Cfloat,trans1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_translate,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in.arr,trans0,trans1,odim0,odim1,method))
    AFArray{T,N}(out[])
end

export translate

function scale{T,N}(_in::AFArray{T,N},scale0::Cfloat,scale1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_scale,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type),out,_in.arr,scale0,scale1,odim0,odim1,method))
    AFArray{T,N}(out[])
end

export scale

function skew{T,N}(_in::AFArray{T,N},skew0::Cfloat,skew1::Cfloat,odim0::dim_t,odim1::dim_t,method::af_interp_type,inverse::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_skew,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,dim_t,dim_t,af_interp_type,Bool),out,_in.arr,skew0,skew1,odim0,odim1,method,inverse))
    AFArray{T,N}(out[])
end

export skew

function histogram{T,N}(_in::AFArray{T,N},nbins::Integer,minval::Real,maxval::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_histogram,af_lib),af_err,(Ptr{af_array},af_array,UInt32,Cdouble,Cdouble),out,_in.arr,UInt32(nbins),Cdouble(minval),Cdouble(maxval)))
    AFArray{T,N}(out[])
end

export histogram

function dilate(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dilate,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

export dilate

function dilate3(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dilate3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

export dilate3

function erode(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_erode,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

export erode

function erode3(_in::AFArray,mask::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_erode3,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,mask.arr))
    AFArray!(out[])
end

export erode3

function bilateral{T,N}(_in::AFArray{T,N},spatial_sigma::Cfloat,chromatic_sigma::Cfloat,isColor::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_bilateral,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Bool),out,_in.arr,spatial_sigma,chromatic_sigma,isColor))
    AFArray{T,N}(out[])
end

export bilateral

function mean_shift{T,N}(_in::AFArray{T,N},spatial_sigma::Cfloat,chromatic_sigma::Cfloat,iter::Integer,is_color::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean_shift,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,UInt32,Bool),out,_in.arr,spatial_sigma,chromatic_sigma,UInt32(iter),is_color))
    AFArray{T,N}(out[])
end

export mean_shift

function minfilt{T,N}(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_minfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

export minfilt

function maxfilt{T,N}(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_maxfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

export maxfilt

function regions{T,N}(_in::AFArray{T,N},connectivity::af_connectivity,ty::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_regions,af_lib),af_err,(Ptr{af_array},af_array,af_connectivity,af_dtype),out,_in.arr,connectivity,ty))
    AFArray{T,N}(out[])
end

export regions

function sobel_operator(img::AFArray,ker_size::Integer)
    dx = RefValue{af_array}(0)
    dy = RefValue{af_array}(0)
    _error(ccall((:af_sobel_operator,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,UInt32),dx,dy,img.arr,UInt32(ker_size)))
    (AFArray!(dx[]),AFArray!(dy[]))
end

export sobel_operator

function rgb2gray{T,N}(_in::AFArray{T,N},rPercent::Cfloat,gPercent::Cfloat,bPercent::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2gray,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in.arr,rPercent,gPercent,bPercent))
    AFArray{T,N}(out[])
end

export rgb2gray

function gray2rgb{T,N}(_in::AFArray{T,N},rFactor::Cfloat,gFactor::Cfloat,bFactor::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_gray2rgb,af_lib),af_err,(Ptr{af_array},af_array,Cfloat,Cfloat,Cfloat),out,_in.arr,rFactor,gFactor,bFactor))
    AFArray{T,N}(out[])
end

export gray2rgb

function hist_equal(_in::AFArray,hist::AFArray)
    out = RefValue{af_array}(0)
    _error(ccall((:af_hist_equal,af_lib),af_err,(Ptr{af_array},af_array,af_array),out,_in.arr,hist.arr))
    AFArray!(out[])
end

export hist_equal

function gaussian_kernel(rows::Integer,cols::Integer,sigma_r::Real,sigma_c::Real)
    out = RefValue{af_array}(0)
    _error(ccall((:af_gaussian_kernel,af_lib),af_err,(Ptr{af_array},Cint,Cint,Cdouble,Cdouble),out,Cint(rows),Cint(cols),Cdouble(sigma_r),Cdouble(sigma_c)))
    AFArray!(out[])
end

export gaussian_kernel

function hsv2rgb{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_hsv2rgb,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export hsv2rgb

function rgb2hsv{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2hsv,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export rgb2hsv

function color_space{T,N}(image::AFArray{T,N},to::af_cspace_t,from::af_cspace_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_color_space,af_lib),af_err,(Ptr{af_array},af_array,af_cspace_t,af_cspace_t),out,image.arr,to,from))
    AFArray{T,N}(out[])
end

export color_space

function unwrap{T,N}(_in::AFArray{T,N},wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_unwrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in.arr,wx,wy,sx,sy,px,py,is_column))
    AFArray{T,N}(out[])
end

export unwrap

function wrap{T,N}(_in::AFArray{T,N},ox::dim_t,oy::dim_t,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_wrap,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in.arr,ox,oy,wx,wy,sx,sy,px,py,is_column))
    AFArray{T,N}(out[])
end

export wrap

function sat{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sat,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sat

function ycbcr2rgb{T,N}(_in::AFArray{T,N},standard::af_ycc_std)
    out = RefValue{af_array}(0)
    _error(ccall((:af_ycbcr2rgb,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in.arr,standard))
    AFArray{T,N}(out[])
end

export ycbcr2rgb

function rgb2ycbcr{T,N}(_in::AFArray{T,N},standard::af_ycc_std)
    out = RefValue{af_array}(0)
    _error(ccall((:af_rgb2ycbcr,af_lib),af_err,(Ptr{af_array},af_array,af_ycc_std),out,_in.arr,standard))
    AFArray{T,N}(out[])
end

export rgb2ycbcr

function moments{T,N}(_in::AFArray{T,N},moment::af_moment_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_moments,af_lib),af_err,(Ptr{af_array},af_array,af_moment_type),out,_in.arr,moment))
    AFArray{T,N}(out[])
end

export moments

function moments_all(_in::AFArray,moment::af_moment_type)
    out = RefValue{Cdouble}(0)
    _error(ccall((:af_moments_all,af_lib),af_err,(Ptr{Cdouble},af_array,af_moment_type),out,_in.arr,moment))
    out[]
end

export moments_all

function svd(_in::AFArray)
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    _error(ccall((:af_svd,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),u,s,vt,_in.arr))
    (AFArray!(u[]),AFArray!(s[]),AFArray!(vt[]))
end

export svd

function svd_inplace(_in::AFArray)
    u = RefValue{af_array}(0)
    s = RefValue{af_array}(0)
    vt = RefValue{af_array}(0)
    _error(ccall((:af_svd_inplace,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),u,s,vt,_in.arr))
    (AFArray!(u[]),AFArray!(s[]),AFArray!(vt[]))
end

export svd_inplace

function lu(_in::AFArray)
    lower = RefValue{af_array}(0)
    upper = RefValue{af_array}(0)
    pivot = RefValue{af_array}(0)
    _error(ccall((:af_lu,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),lower,upper,pivot,_in.arr))
    (AFArray!(lower[]),AFArray!(upper[]),AFArray!(pivot[]))
end

export lu

function lu_inplace{T,N}(_in::AFArray{T,N},is_lapack_piv::Bool)
    pivot = RefValue{af_array}(0)
    _error(ccall((:af_lu_inplace,af_lib),af_err,(Ptr{af_array},af_array,Bool),pivot,_in.arr,is_lapack_piv))
    AFArray{T,N}(pivot[])
end

export lu_inplace

function qr(_in::AFArray)
    q = RefValue{af_array}(0)
    r = RefValue{af_array}(0)
    tau = RefValue{af_array}(0)
    _error(ccall((:af_qr,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},af_array),q,r,tau,_in.arr))
    (AFArray!(q[]),AFArray!(r[]),AFArray!(tau[]))
end

export qr

function qr_inplace{T,N}(_in::AFArray{T,N})
    tau = RefValue{af_array}(0)
    _error(ccall((:af_qr_inplace,af_lib),af_err,(Ptr{af_array},af_array),tau,_in.arr))
    AFArray{T,N}(tau[])
end

export qr_inplace

function cholesky(_in::AFArray,is_upper::Bool)
    out = RefValue{af_array}(0)
    info = RefValue{Cint}(0)
    _error(ccall((:af_cholesky,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,Bool),out,info,_in.arr,is_upper))
    (AFArray!(out[]),info[])
end

export cholesky

function cholesky_inplace(_in::AFArray,is_upper::Bool)
    info = RefValue{Cint}(0)
    _error(ccall((:af_cholesky_inplace,af_lib),af_err,(Ptr{Cint},af_array,Bool),info,_in.arr,is_upper))
    info[]
end

export cholesky_inplace

function solve(a::AFArray,b::AFArray,options::af_mat_prop)
    x = RefValue{af_array}(0)
    _error(ccall((:af_solve,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_mat_prop),x,a.arr,b.arr,options))
    AFArray!(x[])
end

export solve

function solve_lu(a::AFArray,piv::AFArray,b::AFArray,options::af_mat_prop)
    x = RefValue{af_array}(0)
    _error(ccall((:af_solve_lu,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_mat_prop),x,a.arr,piv.arr,b.arr,options))
    AFArray!(x[])
end

export solve_lu

function inverse{T,N}(_in::AFArray{T,N},options::af_mat_prop)
    out = RefValue{af_array}(0)
    _error(ccall((:af_inverse,af_lib),af_err,(Ptr{af_array},af_array,af_mat_prop),out,_in.arr,options))
    AFArray{T,N}(out[])
end

export inverse

function rank(_in::AFArray,tol::Real)
    rank = RefValue{UInt32}(0)
    _error(ccall((:af_rank,af_lib),af_err,(Ptr{UInt32},af_array,Cdouble),rank,_in.arr,Cdouble(tol)))
    rank[]
end

export rank

function det(_in::AFArray)
    det_real = RefValue{Cdouble}(0)
    det_imag = RefValue{Cdouble}(0)
    _error(ccall((:af_det,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),det_real,det_imag,_in.arr))
    (det_real[],det_imag[])
end

export det

function norm(_in::AFArray,_type::af_norm_type,p::Real,q::Real)
    out = RefValue{Cdouble}(0)
    _error(ccall((:af_norm,af_lib),af_err,(Ptr{Cdouble},af_array,af_norm_type,Cdouble,Cdouble),out,_in.arr,_type,Cdouble(p),Cdouble(q)))
    out[]
end

export norm

function is_lapack_available()
    out = RefValue{Bool}(0)
    _error(ccall((:af_is_lapack_available,af_lib),af_err,(Ptr{Bool},),out))
    out[]
end

export is_lapack_available

function create_random_engine(rtype::af_random_engine_type,seed::uintl)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_create_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type,uintl),engine,rtype,seed))
    engine[]
end

export create_random_engine

function retain_random_engine(engine::af_random_engine)
    out = RefValue{af_random_engine}(0)
    _error(ccall((:af_retain_random_engine,af_lib),af_err,(Ptr{af_random_engine},af_random_engine),out,engine))
    out[]
end

export retain_random_engine

function random_engine_set_type(rtype::af_random_engine_type)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_random_engine_set_type,af_lib),af_err,(Ptr{af_random_engine},af_random_engine_type),engine,rtype))
    engine[]
end

export random_engine_set_type

function random_engine_get_type(engine::af_random_engine)
    rtype = RefValue{af_random_engine_type}(0)
    _error(ccall((:af_random_engine_get_type,af_lib),af_err,(Ptr{af_random_engine_type},af_random_engine),rtype,engine))
    rtype[]
end

export random_engine_get_type

function random_uniform(ndims::Integer,dims,_type::af_dtype,engine::af_random_engine)
    out = RefValue{af_array}(0)
    _error(ccall((:af_random_uniform,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,UInt32(ndims),dims,_type,engine))
    AFArray!(out[])
end

export random_uniform

function random_normal(ndims::Integer,dims,_type::af_dtype,engine::af_random_engine)
    out = RefValue{af_array}(0)
    _error(ccall((:af_random_normal,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype,af_random_engine),out,UInt32(ndims),dims,_type,engine))
    AFArray!(out[])
end

export random_normal

function random_engine_set_seed(seed::uintl)
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_random_engine_set_seed,af_lib),af_err,(Ptr{af_random_engine},uintl),engine,seed))
    engine[]
end

export random_engine_set_seed

function get_default_random_engine()
    engine = RefValue{af_random_engine}(0)
    _error(ccall((:af_get_default_random_engine,af_lib),af_err,(Ptr{af_random_engine},),engine))
    engine[]
end

export get_default_random_engine

function set_default_random_engine_type(rtype::af_random_engine_type)
    _error(ccall((:af_set_default_random_engine_type,af_lib),af_err,(af_random_engine_type,),rtype))
end

export set_default_random_engine_type

function random_engine_get_seed(engine::af_random_engine)
    seed = RefValue{uintl}(0)
    _error(ccall((:af_random_engine_get_seed,af_lib),af_err,(Ptr{uintl},af_random_engine),seed,engine))
    seed[]
end

export random_engine_get_seed

function release_random_engine(engine::af_random_engine)
    _error(ccall((:af_release_random_engine,af_lib),af_err,(af_random_engine,),engine))
end

export release_random_engine

function randu(ndims::Integer,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_randu,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,_type))
    AFArray!(out[])
end

export randu

function randn(ndims::Integer,dims,_type::af_dtype)
    out = RefValue{af_array}(0)
    _error(ccall((:af_randn,af_lib),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),out,UInt32(ndims),dims,_type))
    AFArray!(out[])
end

export randn

function set_seed(seed::uintl)
    _error(ccall((:af_set_seed,af_lib),af_err,(uintl,),seed))
end

export set_seed

function get_seed()
    seed = RefValue{uintl}(0)
    _error(ccall((:af_get_seed,af_lib),af_err,(Ptr{uintl},),seed))
    seed[]
end

export get_seed

function approx1(_in::AFArray,pos::AFArray,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_approx1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_interp_type,Cfloat),out,_in.arr,pos.arr,method,offGrid))
    AFArray!(out[])
end

export approx1

function approx2(_in::AFArray,pos0::AFArray,pos1::AFArray,method::af_interp_type,offGrid::Cfloat)
    out = RefValue{af_array}(0)
    _error(ccall((:af_approx2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_interp_type,Cfloat),out,_in.arr,pos0.arr,pos1.arr,method,offGrid))
    AFArray!(out[])
end

export approx2

function fft{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

export fft

function fft_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_fft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export fft_inplace

function fft2{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t,odim1::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

export fft2

function fft2_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_fft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export fft2_inplace

function fft3{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t,odim1::dim_t,odim2::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

export fft3

function fft3_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_fft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export fft3_inplace

function ifft{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

export ifft

function ifft_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_ifft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export ifft_inplace

function ifft2{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t,odim1::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

export ifft2

function ifft2_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_ifft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export ifft2_inplace

function ifft3{T,N}(_in::AFArray{T,N},norm_factor::Real,odim0::dim_t,odim1::dim_t,odim2::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

export ifft3

function ifft3_inplace(_in::AFArray,norm_factor::Real)
    _error(ccall((:af_ifft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

export ifft3_inplace

function fft_r2c{T,N}(_in::AFArray{T,N},norm_factor::Real,pad0::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),pad0))
    AFArray{T,N}(out[])
end

export fft_r2c

function fft2_r2c{T,N}(_in::AFArray{T,N},norm_factor::Real,pad0::dim_t,pad1::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1))
    AFArray{T,N}(out[])
end

export fft2_r2c

function fft3_r2c{T,N}(_in::AFArray{T,N},norm_factor::Real,pad0::dim_t,pad1::dim_t,pad2::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1,pad2))
    AFArray{T,N}(out[])
end

export fft3_r2c

function fft_c2r{T,N}(_in::AFArray{T,N},norm_factor::Real,is_odd::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

export fft_c2r

function fft2_c2r{T,N}(_in::AFArray{T,N},norm_factor::Real,is_odd::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

export fft2_c2r

function fft3_c2r{T,N}(_in::AFArray{T,N},norm_factor::Real,is_odd::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

export fft3_c2r

function convolve1(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

export convolve1

function convolve2(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

export convolve2

function convolve3(signal::AFArray,filter::AFArray,mode::af_conv_mode,domain::af_conv_domain)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode,af_conv_domain),out,signal.arr,filter.arr,mode,domain))
    AFArray!(out[])
end

export convolve3

function convolve2_sep(col_filter::AFArray,row_filter::AFArray,signal::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_convolve2_sep,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array,af_conv_mode),out,col_filter.arr,row_filter.arr,signal.arr,mode))
    AFArray!(out[])
end

export convolve2_sep

function fft_convolve1(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

export fft_convolve1

function fft_convolve2(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

export fft_convolve2

function fft_convolve3(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

export fft_convolve3

function fir(b::AFArray,x::AFArray)
    y = RefValue{af_array}(0)
    _error(ccall((:af_fir,af_lib),af_err,(Ptr{af_array},af_array,af_array),y,b.arr,x.arr))
    AFArray!(y[])
end

export fir

function iir(b::AFArray,a::AFArray,x::AFArray)
    y = RefValue{af_array}(0)
    _error(ccall((:af_iir,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_array),y,b.arr,a.arr,x.arr))
    AFArray!(y[])
end

export iir

function medfilt{T,N}(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

export medfilt

function medfilt1{T,N}(_in::AFArray{T,N},wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt1,af_lib),af_err,(Ptr{af_array},af_array,dim_t,af_border_type),out,_in.arr,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

export medfilt1

function medfilt2{T,N}(_in::AFArray{T,N},wind_length::dim_t,wind_width::dim_t,edge_pad::af_border_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_medfilt2,af_lib),af_err,(Ptr{af_array},af_array,dim_t,dim_t,af_border_type),out,_in.arr,wind_length,wind_width,edge_pad))
    AFArray{T,N}(out[])
end

export medfilt2

function set_fft_plan_cache_size(cache_size::Csize_t)
    _error(ccall((:af_set_fft_plan_cache_size,af_lib),af_err,(Csize_t,),cache_size))
end

export set_fft_plan_cache_size

function create_sparse_array(nRows::dim_t,nCols::dim_t,values::AFArray,rowIdx::AFArray,colIdx::AFArray,stype::af_storage)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,af_array,af_array,af_array,af_storage),out,nRows,nCols,values.arr,rowIdx.arr,colIdx.arr,stype))
    AFArray!(out[])
end

export create_sparse_array

function create_sparse_array_from_ptr(nRows::dim_t,nCols::dim_t,nNZ::dim_t,values,rowIdx,colIdx,_type::af_dtype,stype::af_storage,src::af_source)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array_from_ptr,af_lib),af_err,(Ptr{af_array},dim_t,dim_t,dim_t,Ptr{Void},Ptr{Cint},Ptr{Cint},af_dtype,af_storage,af_source),out,nRows,nCols,nNZ,values,rowIdx,colIdx,_type,stype,src))
    AFArray!(out[])
end

export create_sparse_array_from_ptr

function create_sparse_array_from_dense(dense::AFArray,stype::af_storage)
    out = RefValue{af_array}(0)
    _error(ccall((:af_create_sparse_array_from_dense,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,dense.arr,stype))
    AFArray!(out[])
end

export create_sparse_array_from_dense

function sparse_convert_to{T,N}(_in::AFArray{T,N},destStorage::af_storage)
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_convert_to,af_lib),af_err,(Ptr{af_array},af_array,af_storage),out,_in.arr,destStorage))
    AFArray{T,N}(out[])
end

export sparse_convert_to

function sparse_to_dense{T,N}(sparse::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_to_dense,af_lib),af_err,(Ptr{af_array},af_array),out,sparse.arr))
    AFArray{T,N}(out[])
end

export sparse_to_dense

function sparse_get_info(_in::AFArray)
    values = RefValue{af_array}(0)
    rowIdx = RefValue{af_array}(0)
    colIdx = RefValue{af_array}(0)
    stype = RefValue{af_storage}(0)
    _error(ccall((:af_sparse_get_info,af_lib),af_err,(Ptr{af_array},Ptr{af_array},Ptr{af_array},Ptr{af_storage},af_array),values,rowIdx,colIdx,stype,_in.arr))
    (AFArray!(values[]),AFArray!(rowIdx[]),AFArray!(colIdx[]),stype[])
end

export sparse_get_info

function sparse_get_values{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_values,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sparse_get_values

function sparse_get_row_idx{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_row_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sparse_get_row_idx

function sparse_get_col_idx{T,N}(_in::AFArray{T,N})
    out = RefValue{af_array}(0)
    _error(ccall((:af_sparse_get_col_idx,af_lib),af_err,(Ptr{af_array},af_array),out,_in.arr))
    AFArray{T,N}(out[])
end

export sparse_get_col_idx

function sparse_get_nnz(_in::AFArray)
    out = RefValue{dim_t}(0)
    _error(ccall((:af_sparse_get_nnz,af_lib),af_err,(Ptr{dim_t},af_array),out,_in.arr))
    out[]
end

export sparse_get_nnz

function sparse_get_storage(_in::AFArray)
    out = RefValue{af_storage}(0)
    _error(ccall((:af_sparse_get_storage,af_lib),af_err,(Ptr{af_storage},af_array),out,_in.arr))
    out[]
end

export sparse_get_storage

function mean{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in.arr,dim))
    AFArray{T,N}(out[])
end

export mean

function mean_weighted(_in::AFArray,weights::AFArray,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_mean_weighted,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t),out,_in.arr,weights.arr,dim))
    AFArray!(out[])
end

export mean_weighted

function var{T,N}(_in::AFArray{T,N},isbiased::Bool,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_var,af_lib),af_err,(Ptr{af_array},af_array,Bool,dim_t),out,_in.arr,isbiased,dim))
    AFArray{T,N}(out[])
end

export var

function var_weighted(_in::AFArray,weights::AFArray,dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_var_weighted,af_lib),af_err,(Ptr{af_array},af_array,af_array,dim_t),out,_in.arr,weights.arr,dim))
    AFArray!(out[])
end

export var_weighted

function stdev{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_stdev,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in.arr,dim))
    AFArray{T,N}(out[])
end

export stdev

function cov(X::AFArray,Y::AFArray,isbiased::Bool)
    out = RefValue{af_array}(0)
    _error(ccall((:af_cov,af_lib),af_err,(Ptr{af_array},af_array,af_array,Bool),out,X.arr,Y.arr,isbiased))
    AFArray!(out[])
end

export cov

function median{T,N}(_in::AFArray{T,N},dim::dim_t)
    out = RefValue{af_array}(0)
    _error(ccall((:af_median,af_lib),af_err,(Ptr{af_array},af_array,dim_t),out,_in.arr,dim))
    AFArray{T,N}(out[])
end

export median

function mean_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_mean_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export mean_all

function mean_all_weighted(_in::AFArray,weights::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_mean_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),real,imag,_in.arr,weights.arr))
    (real[],imag[])
end

export mean_all_weighted

function var_all(_in::AFArray,isbiased::Bool)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_var_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,Bool),realVal,imagVal,_in.arr,isbiased))
    (realVal[],imagVal[])
end

export var_all

function var_all_weighted(_in::AFArray,weights::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_var_all_weighted,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,_in.arr,weights.arr))
    (realVal[],imagVal[])
end

export var_all_weighted

function stdev_all(_in::AFArray)
    real = RefValue{Cdouble}(0)
    imag = RefValue{Cdouble}(0)
    _error(ccall((:af_stdev_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),real,imag,_in.arr))
    (real[],imag[])
end

export stdev_all

function median_all(_in::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_median_all,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array),realVal,imagVal,_in.arr))
    (realVal[],imagVal[])
end

export median_all

function corrcoef(X::AFArray,Y::AFArray)
    realVal = RefValue{Cdouble}(0)
    imagVal = RefValue{Cdouble}(0)
    _error(ccall((:af_corrcoef,af_lib),af_err,(Ptr{Cdouble},Ptr{Cdouble},af_array,af_array),realVal,imagVal,X.arr,Y.arr))
    (realVal[],imagVal[])
end

export corrcoef

function fast(_in::AFArray,thr::Cfloat,arc_length::Integer,non_max::Bool,feature_ratio::Cfloat,edge::Integer)
    out = RefValue{af_features}(0)
    _error(ccall((:af_fast,af_lib),af_err,(Ptr{af_features},af_array,Cfloat,UInt32,Bool,Cfloat,UInt32),out,_in.arr,thr,UInt32(arc_length),non_max,feature_ratio,UInt32(edge)))
    out[]
end

export fast

function harris(_in::AFArray,max_corners::Integer,min_response::Cfloat,sigma::Cfloat,block_size::Integer,k_thr::Cfloat)
    out = RefValue{af_features}(0)
    _error(ccall((:af_harris,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,UInt32,Cfloat),out,_in.arr,UInt32(max_corners),min_response,sigma,UInt32(block_size),k_thr))
    out[]
end

export harris

function orb(_in::AFArray,fast_thr::Cfloat,max_feat::Integer,scl_fctr::Cfloat,levels::Integer,blur_img::Bool)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_orb,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,Cfloat,UInt32,Cfloat,UInt32,Bool),feat,desc,_in.arr,fast_thr,UInt32(max_feat),scl_fctr,UInt32(levels),blur_img))
    (feat[],AFArray!(desc[]))
end

export orb

function sift(_in::AFArray,n_layers::Integer,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_sift,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in.arr,UInt32(n_layers),contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],AFArray!(desc[]))
end

export sift

function gloh(_in::AFArray,n_layers::Integer,contrast_thr::Cfloat,edge_thr::Cfloat,init_sigma::Cfloat,double_input::Bool,intensity_scale::Cfloat,feature_ratio::Cfloat)
    feat = RefValue{af_features}(0)
    desc = RefValue{af_array}(0)
    _error(ccall((:af_gloh,af_lib),af_err,(Ptr{af_features},Ptr{af_array},af_array,UInt32,Cfloat,Cfloat,Cfloat,Bool,Cfloat,Cfloat),feat,desc,_in.arr,UInt32(n_layers),contrast_thr,edge_thr,init_sigma,double_input,intensity_scale,feature_ratio))
    (feat[],AFArray!(desc[]))
end

export gloh

function hamming_matcher(query::AFArray,train::AFArray,dist_dim::dim_t,n_dist::Integer)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    _error(ccall((:af_hamming_matcher,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32),idx,dist,query.arr,train.arr,dist_dim,UInt32(n_dist)))
    (AFArray!(idx[]),AFArray!(dist[]))
end

export hamming_matcher

function nearest_neighbour(query::AFArray,train::AFArray,dist_dim::dim_t,n_dist::Integer,dist_type::af_match_type)
    idx = RefValue{af_array}(0)
    dist = RefValue{af_array}(0)
    _error(ccall((:af_nearest_neighbour,af_lib),af_err,(Ptr{af_array},Ptr{af_array},af_array,af_array,dim_t,UInt32,af_match_type),idx,dist,query.arr,train.arr,dist_dim,UInt32(n_dist),dist_type))
    (AFArray!(idx[]),AFArray!(dist[]))
end

export nearest_neighbour

function match_template(search_img::AFArray,template_img::AFArray,m_type::af_match_type)
    out = RefValue{af_array}(0)
    _error(ccall((:af_match_template,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_match_type),out,search_img.arr,template_img.arr,m_type))
    AFArray!(out[])
end

export match_template

function susan(_in::AFArray,radius::Integer,diff_thr::Cfloat,geom_thr::Cfloat,feature_ratio::Cfloat,edge::Integer)
    out = RefValue{af_features}(0)
    _error(ccall((:af_susan,af_lib),af_err,(Ptr{af_features},af_array,UInt32,Cfloat,Cfloat,Cfloat,UInt32),out,_in.arr,UInt32(radius),diff_thr,geom_thr,feature_ratio,UInt32(edge)))
    out[]
end

export susan

function dog{T,N}(_in::AFArray{T,N},radius1::Integer,radius2::Integer)
    out = RefValue{af_array}(0)
    _error(ccall((:af_dog,af_lib),af_err,(Ptr{af_array},af_array,Cint,Cint),out,_in.arr,Cint(radius1),Cint(radius2)))
    AFArray{T,N}(out[])
end

export dog

function homography(x_src::AFArray,y_src::AFArray,x_dst::AFArray,y_dst::AFArray,htype::af_homography_type,inlier_thr::Cfloat,iterations::Integer,otype::af_dtype)
    H = RefValue{af_array}(0)
    inliers = RefValue{Cint}(0)
    _error(ccall((:af_homography,af_lib),af_err,(Ptr{af_array},Ptr{Cint},af_array,af_array,af_array,af_array,af_homography_type,Cfloat,UInt32,af_dtype),H,inliers,x_src.arr,y_src.arr,x_dst.arr,y_dst.arr,htype,inlier_thr,UInt32(iterations),otype))
    (AFArray!(H[]),inliers[])
end

export homography
