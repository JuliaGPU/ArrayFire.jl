# Julia wrapper for header: /usr/local/include/arrayfire.h
# Automatically generated using Clang.jl wrap_c, version 0.0.0


function af_sum(out::Base.Ref, _in::AFArray, dim::Integer)
    err= ccall((:af_sum,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_sum_nan(out::Base.Ref, _in::AFArray, dim::Int, nanval::Real)
    err = ccall((:af_sum_nan, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Cdouble),
                out, _in.ptr, dim, nanval)
    err == 0 || throwAFerror(err)
end

function af_product(out::Base.Ref, _in::AFArray, dim::Integer)
    err= ccall((:af_product,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_product_nan(out::Base.Ref, _in::AFArray, dim::Int, nanval::Real)
    err = ccall((:af_product_nan, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Cdouble),
                out, _in.ptr, dim, nanval)
    err == 0 || throwAFerror(err)
end

function af_min(out::Base.Ref, _in::AFArray, dim::Integer)
    err = ccall((:af_min,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_max(out::Base.Ref, _in::AFArray, dim::Integer)
    err = ccall((:af_max,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_all_true(out::Base.Ref, _in::AFArray, dim::Integer)
    err = ccall((:af_all_true,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_any_true(out::Base.Ref, _in::AFArray, dim::Integer)
    err = ccall((:af_any_true,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_count(out::Base.Ref,_in::AFArray,dim::Integer)
    err = ccall((:af_count,af_lib), 
                Cint ,(Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr , dim)
    err == 0 || throwAFerror(err)
end

function af_sum_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err= ccall((:af_sum_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_sum_nan_all(real::Base.Ref, imag::Base.Ref, _in::AFArray, nanval::Real)
    err = ccall((:af_sum_nan_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Cdouble),
                real, imag, _in.ptr, nanval)
    err == 0 || throwAFerror(err)
end

function af_product_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err= ccall((:af_product_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_product_nan_all(real::Base.Ref, imag::Base.Ref, _in::AFArray, nanval::Real)
    err = ccall((:af_product_nan_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Cdouble),
                real, imag, _in.ptr, nanval)
    err == 0 || throwAFerror(err)
end

function af_min_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_min_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_max_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_max_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_all_true_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_all_true_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_any_true_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_any_true_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_count_all(real::Base.Ref, imag::Base.Ref,_in::AFArray)
    err = ccall((:af_count_all, af_lib), 
                Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_imin(out::Base.Ref, idx::Base.Ref, _in::AFArray, dim::Int)
    err = ccall((:af_imin, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint),
                out, idx, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_imax(out::Base.Ref, idx::Base.Ref, _in::AFArray, dim::Int)
    err = ccall((:af_imax, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint),
                out, idx, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_imin_all(real::Base.Ref, imag::Base.Ref, idx::Base.Ref,  _in::AFArray)
    err = ccall((:af_imin_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt32}, Ptr{Void}),
                real, imag, idx, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_imax_all(real::Base.Ref, imag::Base.Ref, idx::Base.Ref, _in::AFArray)
    err = ccall((:af_imax_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{UInt32}, Ptr{Void}),
                real, imag, idx, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_accum(out::Base.Ref, _in::AFArray, dim::Cint)
    err = ccall((:af_accum, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint),
                out, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_where(idx::Base.Ref, _in::AFArray)
    err = ccall((:af_where, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}),
                idx, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_diff1(out::Base.Ref,_in::AFArray, dim::Cint)
    err = ccall((:af_diff1, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint),
                out, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_diff2(out::Base.Ref, _in::AFArray, dim::Cint)
    err = ccall((:af_diff2, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint),
                out, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_sort(out::Base.Ref, _in::AFArray, dim::UInt32, isAscending::Bool)
    err = ccall((:af_sort, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, UInt32, Bool), 
                out, _in.ptr, dim, isAscending)
    err == 0 || throwAFerror(err)
end

function af_sort_index(out::Base.Ref, indices::Base.Ref, _in::AFArray, dim::UInt32, isAscending::Bool)
    err = ccall((:af_sort_index, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32, Bool), 
                out, indices, _in.ptr, dim, isAscending)
    err == 0 || throwAFerror(err)
end

function af_sort_by_key(out_keys::Base.Ref, out_values::Base.Ref, keys::AFArray, 
                        values::AFArray, dim::UInt32, isAscending::Bool)
    err = ccall((:af_sort_by_key, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32, Bool), 
                out_keys , out_values, keys.ptr, values.ptr, dim, isAscending)
    err == 0 || throwAFerror(err)
end

function af_set_unique(out::Base.Ref, _in::AFArray, is_sorted::Bool)
    err = ccall((:af_set_unique, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Bool), 
                out, _in.ptr, is_sorted)
    err == 0 || throwAFerror(err)
end

function af_set_union(out::Base.Ref, first::AFArray, second::AFArray, is_unique::Bool)
    err = ccall((:af_set_union, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool), 
                out, first.ptr, second.ptr, is_unique)
    err == 0 || throwAFerror(err)
end

function af_set_intersect(out::Base.Ref, first::AFArray, second::AFArray, is_unique::Bool)
    err = ccall((:af_set_intersect, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool), 
                out , first.ptr, second.ptr, is_unique)
    err == 0 || throwAFerror(err)
end

function af_add(out::Base.Ref, lhs::AFArray, rhs::AFArray,batch::Bool = true)
    err = ccall((:af_add, af_lib), 
                Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_sub(out::Base.Ref, lhs::AFArray, rhs::AFArray,batch::Bool = true)
    err = ccall((:af_sub, af_lib), 
                Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_mul(out::Base.Ref, lhs::AFArray, rhs::AFArray,batch::Bool = true)
    err = ccall((:af_mul, af_lib), 
                Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_div(out::Base.Ref, lhs::AFArray, rhs::AFArray,batch::Bool = true)
    err = ccall((:af_div, af_lib), 
                Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_lt(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_lt, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_gt(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_gt, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_le(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_le, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_ge(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_ge, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_eq(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_eq, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_neq(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_neq, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_and(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_and, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_or(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_or, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_not(out::Base.Ref, _in::AFArray)
    err = ccall((:af_not, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_bitand(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_bitand, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_bitor(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_bitor, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_bitxor(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_bitxor, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_bitshiftl(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_bitshiftl, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_bitshiftr(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_bitshiftr, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool), 
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_cast(out::Base.Ref, _in::AFArray, _type::DataType)
    err = ccall((:af_cast, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cuint),
                out, _in.ptr, aftype(_type))
    err == 0 || throwAFerror(err)
end

function af_minof(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_minof, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_maxof(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_maxof, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_rem(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_rem, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_mod(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_mod, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_abs(out::Base.Ref, _in::AFArray)
    err = ccall((:af_abs, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_arg(out::Base.Ref, _in::AFArray)
    err = ccall((:af_arg, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_sign(out::Base.Ref, _in::AFArray)
    err = ccall((:af_sign, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_round(out::Base.Ref, _in::AFArray)
    err = ccall((:af_round, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_trunc(out::Base.Ref, _in::AFArray)
    err = ccall((:af_trunc, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_floor(out::Base.Ref, _in::AFArray)
    err = ccall((:af_floor, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_ceil(out::Base.Ref, _in::AFArray)
    err = ccall((:af_ceil, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_hypot(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_hypot, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_sin(out::Base.Ref,_in::AFArray)
    err = ccall((:af_sin, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_cos(out,_in::AFArray)
    err = ccall((:af_cos, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_tan(out,_in::AFArray)
    err = ccall((:af_tan, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_asin(out::Base.Ref, _in::AFArray)
    err = ccall((:af_asin, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_acos(out::Base.Ref, _in::AFArray)
    err = ccall((:af_acos, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_atan(out::Base.Ref, _in::AFArray)
    err = ccall((:af_atan, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_atan2(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_atan2, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_cplx2(out::Base.Ref ,lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_cplx2,af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool), out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end


function af_cplx(out::Base.Ref, _in::AFArray)
    err = ccall((:af_cplx, af_lib), Cint, 
                (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_real(out::Base.Ref, _in::AFArray)
    err = ccall((:af_real, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_imag(out::Base.Ref, _in::AFArray)
    err = ccall((:af_imag, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_conjg(out::Base.Ref, _in::AFArray)
    err = ccall((:af_conjg, af_lib), Cint, 
                (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_sinh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_sinh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_cosh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_cosh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_tanh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_tanh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_asinh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_asinh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_acosh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_acosh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_atanh(out::Base.Ref,_in::AFArray)
    err = ccall((:af_atanh, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_root(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_root, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_pow(out::Base.Ref, lhs::AFArray, rhs::AFArray, batch::Bool)
    err = ccall((:af_pow, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool),
                out, lhs.ptr, rhs.ptr, batch)
    err == 0 || throwAFerror(err)
end

function af_pow2(out,_in::AFArray)
    err = ccall((:af_pow2,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_exp(out,_in::AFArray)
    err = ccall((:af_exp,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_sigmoid(out::Base.Ref, _in::AFArray)
    err = ccall((:af_sigmoid,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_expm1(out,_in::AFArray)
    err = ccall((:af_expm1,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_erf(out::Base.Ref,_in::AFArray)
    err = ccall((:af_erf,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_erfc(out::Base.Ref,_in::AFArray)
    err = ccall((:af_erfc,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_log(out::Base.Ref,_in::AFArray)
    err = ccall((:af_log,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_log1p(out::Base.Ref,_in::AFArray)
    err = ccall((:af_log1p,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_log10(out::Base.Ref,_in::AFArray)
    err = ccall((:af_log10,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_log2(out::Base.Ref,_in::AFArray)
    err = ccall((:af_log2,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_sqrt(out::Base.Ref,_in::AFArray)
    err = ccall((:af_sqrt,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_cbrt(out::Base.Ref,_in::AFArray)
    err = ccall((:af_cbrt,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_factorial(out::Base.Ref,_in::AFArray)
    err = ccall((:af_factorial,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_tgamma(out::Base.Ref,_in::AFArray)
    err = ccall((:af_tgamma,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_lgamma(out::Base.Ref,_in::AFArray)
    err = ccall((:af_lgamma,af_lib), Cint,
            (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_iszero(out::Base.Ref, _in::AFArray)
    err = ccall((:af_iszero, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_isinf(out::Base.Ref, _in::AFArray)
    err = ccall((:af_isinf, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_isnan(out::Base.Ref, _in::AFArray)
    err = ccall((:af_isnan, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_make_seq(_begin::Real, _end::Real, step::Real)
    ccall((:af_make_seq, af_lib), seq, 
                (Cdouble,Cdouble,Cdouble), _begin, _end, step)
end

function af_print_array(arr::AFArray)
    err = ccall((:af_print_array, af_lib), Cint,(Ptr{Void},),arr.ptr)
    err == 0 || throwAFerror(err)
end

#=
function af_print_array_gen(exp,arr::af_array,precision::Cint)
    ccall((:af_print_array_gen,util),af_err,(Cstring,af_array,Cint),exp,arr,precision)
end

function af_save_array(index,key,arr::af_array,filename,append::Bool)
    ccall((:af_save_array,util),af_err,(Ptr{Cint},Cstring,af_array,Cstring,Bool),index,key,arr,filename,append)
end

function af_read_array_index(out,filename,index::UInt32)
    ccall((:af_read_array_index,util),af_err,(Ptr{af_array},Cstring,UInt32),out,filename,index)
end

function af_read_array_key(out,filename,key)
    ccall((:af_read_array_key,util),af_err,(Ptr{af_array},Cstring,Cstring),out,filename,key)
end

function af_read_array_key_check(index,filename,key)
    ccall((:af_read_array_key_check,util),af_err,(Ptr{Cint},Cstring,Cstring),index,filename,key)
end

function af_array_to_string(output,exp,arr::af_array,precision::Cint,transpose::Bool)
    ccall((:af_array_to_string,util),af_err,(Ptr{Cstring},Cstring,af_array,Cint,Bool),output,exp,arr,precision,transpose)
end

function af_example_function(out,_in::af_array,param::af_someenum_t)
    ccall((:af_example_function,util),af_err,(Ptr{af_array},af_array,af_someenum_t),out,_in,param)
end
=#

function af_get_version(major::Base.Ref, minor::Base.Ref, patch::Base.Ref)
    err = ccall((:af_get_version, af_lib), Cint,
                    (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
                    major, minor, patch)
    err == 0 || throwAFerror(err)
end

#=
function af_get_revision()
    ccall((:af_get_revision,util),Cstring,())
end
=#

function af_index(out::Base.Ref, _in::AFArray, ndims::UInt32, index::Union{seq, Vector{seq}})
    err = ccall((:af_index, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, UInt32, Ptr{seq}),
                out, _in.ptr, ndims, pointer(index))
    err == 0 || throwAFerror(err)
end

function af_lookup(out::Base.Ref, _in::AFArray, indices::AFArray, dim::UInt32)
    err = ccall((:af_lookup, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32),
                out, _in.ptr, indices, dim)
    err == 0 || throwAFerror(err)
end

function af_assign_seq(out::Base.Ref, lhs::AFArray, ndims::UInt32, indices::Union{seq, Vector{seq}}, rhs::AFArray)
    err = ccall((:af_assign_seq, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, UInt32, Ptr{seq},Ptr{Void}),
                out, lhs.ptr, ndims, pointer(indices), rhs.ptr)
    err == 0 || throwAFerror(err)
end

function af_index_gen(out::Base.Ref, _in::AFArray, ndims::Int, indices::index)
    err = ccall((:af_index_gen, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Ptr{Void}),
                out, _in.ptr, ndims, indices.ptr)
    err == 0 || throwAFerror(err)
end

function af_assign_gen(out, lhs::AFArray, ndims::Int, indices::index, rhs::AFArray)
    err = ccall((:af_assign_gen, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Ptr{Void}, Ptr{Void}),
                out, lhs.ptr, ndims, indices.ptr, rhs.ptr)
    err == 0 || throwAFerror(err)
end

function af_create_indexers(indexers::Base.Ref)
    err = ccall((:af_create_indexers, af_lib), Cint,
                (Ptr{Void}, ), indexers)
    err == 0 || throwAFerror(err)
end

function af_set_array_indexer(indexer::index, idx::AFArray, dim::Int)
    err = ccall((:af_set_array_indexer, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint), indexer.ptr, idx.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_set_seq_indexer(indexer::index, idx::seq, dim::Int, is_batch::Bool)
    err = ccall((:af_set_seq_indexer, af_lib), Cint, 
                (Ptr{Void}, Ptr{seq}, Cint, Bool),
                indexer.ptr, pointer_from_objref(idx), dim, is_batch)
    err == 0 || throwAFerror(err)
end

function af_set_seq_param_indexer(indexer::Base.Ref, _begin::Real, _end::Real, step::Real, dim::Int, is_batch::Bool)
    err = ccall((:af_set_seq_param_indexer, af_lib), Cint,
                (Ptr{Void}, Cdouble, Cdouble, Cdouble, Cint, Bool),
                indexer, _begin, _end, step, dim, is_batch)
    err == 0 || throwAFerror(err)
end

function af_release_indexers(indexers::index)
    err = ccall((:af_release_indexers, af_lib), Cint,
                (Ptr{Void}, ), indexers.ptr)
    err == 0 || throwAFerror(err)
end

function af_create_array(arr::Base.Ref, data::Ptr{Void}, ndims::UInt32, dims::Vector{Int}, T::DataType)
    err = ccall((:af_create_array, af_lib), Cint,
                    (Ptr{Void}, Ptr{Void}, UInt32, Ptr{Int}, Cuint),
                    arr, data, ndims, dims, aftype(T))
    err == 0 || throwAFerror(err)
end

#=
function af_create_handle(arr,ndims::UInt32,dims,_type::af_dtype)
    ccall((:af_create_handle,array),af_err,(Ptr{af_array},UInt32,Ptr{dim_t},af_dtype),arr,ndims,dims,_type)
end
=#

function af_copy_array(arr::AFArray, _in::Ptr{Void})
    err = ccall((:af_copy_array, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), arr.ptr, _in)
    err == 0 || throwAFerror(err)
end

function af_copy_array(out::Base.Ref, _in::Ptr{Void})
    err = ccall((:af_copy_array, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in)
    err == 0 || throwAFerror(err)
end

#=
function af_write_array(arr::af_array,data,bytes::Cint,src::af_source)
    ccall((:af_write_array,array),af_err,(af_array,Ptr{Void},Cint,af_source),arr,data,bytes,src)
end

function af_get_data_ptr!(ret::Vector{UInt8}, x::AFArray, T::DataType)
    err = ccall((:af_get_data_ptr, af_lib), 
                Cint, (Ptr{T}, Ptr{Ptr{Void}}), 
                pointer(ret), x.ptr)
    err == 0 || throwAFerror(err)
end
=#

function af_release_array(arr::AFArray)
    err = ccall((:af_release_array, af_lib), Cint,
                (Ptr{Void}, ), arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_retain_array(out::Base.Ref, _in::AFArray)
    err = ccall((:af_retain_array, af_lib), Cint,
                (Ptr{Void},Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_get_data_ref_count(use_count::Base.Ref, _in::AFArray)
    err = ccall((:af_get_data_ref_count, af_lib), Cint,
                (Ptr{Cint}, Ptr{Void}), use_count, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_eval(_in::AFArray)
    err = ccall((:af_eval, af_lib), Cint, (Ptr{Void}, ), _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_get_elements(elems::Vector{Int}, arr::AFArray)
    err = ccall((:af_get_elements, af_lib), Cint,
                (Ptr{Int}, Ptr{Void}), elems, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_get_type(_type::Base.Ref, arr::Ptr{Void})
    err = ccall((:af_get_type, af_lib), Cint, 
                (Ptr{Cuint}, Ptr{Void}), _type, arr)
    err == 0 || throwAFerror(err)
end

function af_get_dims!(d1::Base.Ref, d2::Base.Ref, d3::Base.Ref, d4::Base.Ref, arr::AFArray)
    err = ccall((:af_get_dims, af_lib), 
                Cint, 
                (Ptr{Cuint}, Ptr{Cuint}, Ptr{Cuint}, Ptr{Cuint}, Ptr{Void}),
                d1, d2, d3, d4, arr.ptr)
    err == 0 || throwAFerror(err)
end


function af_get_numdims!(n::Base.Ref, ptr::Ptr{Void})
    err = ccall((:af_get_numdims, af_lib),
            Cint, (Ptr{Cuint}, Ptr{Void}),
            n, ptr)
    err == 0 || throwAFerror(err)
end

function af_is_empty(result::Base.Ref, arr::AFArray)
    err = ccall((:af_is_empty, af_lib), Cint, 
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_scalar(result::Base.Ref, arr::AFArray)
    err = ccall((:af_is_scalar, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_row(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_row, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_column(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_column, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_vector(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_vector, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_complex(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_complex, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_real(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_real, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_double(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_double, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_single(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_single, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_realfloating(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_realfloating, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_floating(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_floating, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_integer(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_integer, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_is_bool(result::Base.Ref,arr::AFArray)
    err = ccall((:af_is_bool, af_lib), Cint,
                (Ptr{Bool}, Ptr{Void}), result, arr.ptr)
    err == 0 || throwAFerror(err)
end

function af_set_backend(bknd::Cuint)
    err = ccall((:af_set_backend,af_lib), Cint,
                (Cuint,),bknd)
    err == 0 || throwAFerror(err)
end

function af_get_backend_count(num_backends::Base.Ref)
    err = ccall((:af_get_backend_count, af_lib), Cint,
                (Ptr{UInt32},), num_backends)
    err == 0 || throwAFerror(err)
end

function af_get_available_backends(backends::Base.Ref)
    err = ccall((:af_get_available_backends, af_lib), Cint,
                (Ptr{Cint}, ), backends)
    err == 0 || throwAFerror(err)
end

function af_get_backend_id(backend::Base.Ref,_in::AFArray)
    err = ccall((:af_get_backend_id, af_lib), Cint,
                (Ptr{Cuint}, Ptr{Void}), backend, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_get_device_id(device::Base.Ref,_in::AFArray)
    err = ccall((:af_get_device_id, af_lib), Cint,
                (Ptr{Cint}, Ptr{Void}), device, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_get_active_backend(backend)
    err = ccall((:af_get_active_backend, af_lib), Cint, 
                (Ptr{Cuint},) , backend)
    err == 0 || throwAFerror(err)
end

function af_matmul(out::Base.Ref, lhs::AFArray, rhs::AFArray, optLhs::Int ,optRhs::Int)
    err = ccall((:af_matmul, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, lhs.ptr, rhs.ptr, optLhs, optRhs)
    err == 0 || throwAFerror(err)
end

function af_dot(out::Base.Ref, lhs::AFArray, rhs::AFArray, optLhs::Int, optRhs::Int)
    err = ccall((:af_dot, af_lib), Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, lhs.ptr, rhs.ptr, optLhs, optRhs)
    err == 0 || throwAFerror(err)
end

function af_transpose(out::Base.Ref, _in::AFArray, conjugate::Bool)
    err = ccall((:af_transpose, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Bool), out, _in.ptr, conjugate)
    err == 0 || throwAFerror(err)
end

function af_transpose_inplace(_in::AFArray, conjugate::Bool)
    err = ccall((:af_transpose_inplace, af_lib), Cint,
                (Ptr{Void}, Bool), _in.ptr, conjugate)
    err == 0 || throwAFerror(err)
end

function af_constant!(ptr::Base.Ref, val::Real, n::Int, dims::Vector{Int}, T::DataType)
    err = ccall((:af_constant, af_lib), Cint, 
                (Ptr{Void}, Cdouble, Cint, Ptr{Int}, Cuint),
                ptr, val, n, pointer(dims), aftype(T))
    err == 0 || throwAFerror(err)
end

function af_constant_complex!(ptr::Base.Ref, val::Complex, n::Int, dims::Vector{Int}, T::DataType)
    r = real(val)
    i = imag(val)
    err = ccall((:af_constant_complex, af_lib), Cint, 
                (Ptr{Void}, Cdouble, Cdouble, Cint, Ptr{Int}, Cuint),
                ptr, r, i, n, pointer(dims), aftype(Complex{T}))
    err == 0 || throwAFerror(err)
end

function af_constant_long(arr::Base.Ref, val::Int32, ndims::UInt32, dims::Vector{Int})
    err = ccall((:af_constant_long, af_lib), Cint,
                (Ptr{Void}, Int32, UInt32, Ptr{Int}),
                arr, val, ndims, pointer(dims))
    err == 0 || throwAFerror(err)
end

function af_constant_ulong(arr::Base.Ref, val::UInt32, ndims::UInt32, dims::Vector{Int})
    err = ccall((:af_constant_ulong, af_lib), Cint,
                (Ptr{Void}, UInt32, UInt32, Ptr{Int}),
                arr, val, ndims, pointer(dims))
    err == 0 || throwAFerror(err)
end

function af_range(out::Base.Ref, ndims::UInt32, dims::Vector{Int}, seq_dim::Int, _type::DataType)
    err = ccall((:af_range, af_lib), Cint,
                (Ptr{Void}, UInt32, Ptr{Int}, Cint, Cuint), 
                out, ndims, pointer(dims), seq_dim, aftype(_type))
    err == 0 || throwAFerror(err)
end

function af_iota(out::Base.Ref, ndims::UInt32, dims::Vector{Int}, t_ndims::UInt32, tdims::Vector{Int}, _type::DataType)
    err = ccall((:af_iota, af_lib), Cint, 
                (Ptr{Void}, UInt32, Ptr{Int}, UInt32, Ptr{Int}, Cuint), 
                out, ndims, pointer(dims), t_ndims, pointer(tdims),aftype(_type))
    err == 0 || throwAFerror(err)
end

function af_randu(ptr::Base.Ref, dims::Vector{Int}, T::DataType)
    err = ccall((:af_randu,af_lib), 
                Cint, (Ptr{Ptr{Void}}, Cint, Ptr{Int}, Cuint), 
                ptr , length(dims), pointer(dims), aftype(T))
    err == 0 || throwAFerror(err)
end

function af_randn(ptr::Base.Ref, dims::Vector{Int}, T::DataType)
    err = ccall((:af_randn,af_lib), 
                Cint, (Ptr{Ptr{Void}}, Cint, Ptr{Int}, Cuint), 
                ptr , length(dims), pointer(dims), aftype(T))
    err == 0 || throwAFerror(err)
end

function af_set_seed(seed::Cuint)
    err = ccall((:af_set_seed, af_lib), Cint, 
                (Cuint,), seed)
    err == 0 || throwAFerror(err)
end

function af_get_seed(seed::Base.Ref)
    err = ccall((:af_get_seed, af_lib), Cint,
                (Ptr{Cuint},), seed)
    err == 0 || throwAFerror(err)
end

function af_identity(out::Base.Ref, ndims::UInt32, dims, _type::DataType)
    err = ccall((:af_identity, af_lib), Cint,
                (Ptr{Void}, UInt32, Ptr{Int}, Cuint),
                out, ndims, pointer(dims), aftype(_type))
    err == 0 || throwAFerror(err)
end

function af_diag_create(out::Base.Ref, _in::AFArray, num::Int)
    err = ccall((:af_diag_create, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr, num)
    err == 0 || throwAFerror(err)
end

function af_diag_extract(out::Base.Ref, _in::AFArray, num::Int)
    err = ccall((:af_diag_extract, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr, num)
    err == 0 || throwAFerror(err)
end

function af_join(out::Base.Ref, dim::Int, first::AFArray, second::AFArray)
    err = ccall((:af_join, af_lib), Cint,
                (Ptr{Void}, Cint, Ptr{Void}, Ptr{Void}),
                out, dim, first.ptr, second.ptr)
    err == 0 || throwAFerror(err)
end

function af_join_many(out::Base.Ref, dim::Int, n_arrays::Cuint, inputs::Vector{Ptr{Void}})
    err = ccall((:af_join_many, af_lib), Cint,
                (Ptr{Void}, Cint, UInt32, Ptr{Ptr{Void}}),
                out, dim, n_arrays, pointer(inputs))
    err == 0 || throwAFerror(err)
end

function af_tile(out::Base.Ref, _in::AFArray, x::UInt32, y::UInt32, z::UInt32, w::UInt32)
    err = ccall((:af_tile, af_lib), Cint,
                (Ptr{Void},Ptr{Void},UInt32,UInt32,UInt32,UInt32),
                out, _in.ptr, x, y, z, w)
    err == 0 || throwAFerror(err)
end

function af_reorder(out::Base.Ref, _in::AFArray, x::UInt32, y::UInt32, z::UInt32, w::UInt32)
    err =  ccall((:af_reorder, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, UInt32, UInt32, UInt32, UInt32),
                out, _in.ptr, x, y, z, w)
    err == 0 || throwAFerror(err)
end

function af_shift(out::Base.Ref, _in::AFArray, x::Int, y::Int, z::Int, w::Int)
    err = ccall((:af_shift, af_lib), Cint,
                (Ptr{Void},Ptr{Void},Cint,Cint,Cint,Cint),
                out, _in.ptr, x, y, z, w)
    err == 0 || throwAFerror(err)
end

function af_moddims(out::Base.Ref, _in::AFArray, ndims::UInt32, dims::Vector{Int})
    err = ccall((:af_moddims, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, UInt32, Ptr{Int}),
                out, _in.ptr, ndims, pointer(dims))
    err == 0 || throwAFerror(err)
end

function af_flat(out::Base.Ref, _in::AFArray)
    err = ccall((:af_flat, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_flip(out::Base.Ref, _in::AFArray, dim::UInt32)
    err = ccall((:af_flip, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, UInt32), out, _in.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_lower(out::Base.Ref, _in::AFArray, is_unit_diag::Bool)
    err = ccall((:af_lower, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Bool),
                out, _in.ptr, is_unit_diag)
    err == 0 || throwAFerror(err)
end 

function af_upper(out::Base.Ref, _in::AFArray, is_unit_diag::Bool)
    err = ccall((:af_upper, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Bool),
                out, _in.ptr, is_unit_diag)
    err == 0 || throwAFerror(err)
end

function af_select(out::Base.Ref, cond::AFArray, a::AFArray, b::AFArray)
    err = ccall((:af_select, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
                out, cond.ptr, a.ptr, b.ptr)
    err == 0 || throwAFerror(err)
end

function af_select_scalar_r(out::Base.Ref, cond::AFArray, a::AFArray, b::Real)
    err = ccall((:af_select_scalar_r, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cdouble),
                out, cond.ptr, a.ptr, b)
    err == 0 || throwAFerror(err)
end

function af_select_scalar_l(out::Base.Ref, cond::AFArray, a::Real, b::AFArray)
    err = ccall((:af_select_scalar_l, af_lib), Cint,
                (Ptr{Void},Ptr{Void},Cdouble,Ptr{Void}),
                out, cond.ptr, a, b.ptr)
    err == 0 || throwAFerror(err)
end

function af_replace(a::AFArray, cond::AFArray, b::AFArray)
    err = ccall((:af_replace, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}),
                a.ptr, cond.ptr, b.ptr)
    err == 0 || throwAFerror(err)
end

function af_replace_scalar(a::AFArray, cond::AFArray, b::Real)
    err = ccall((:af_replace_scalar, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble),
                a.ptr, cond.ptr, b.ptr)
    err == 0 || throwAFerror(err)
end

function af_info()
    err = ccall((:af_info, af_lib), Cint, () ) 
    err == 0 || throwAFerror(err)
end

#=
function af_init()
    ccall((:af_init,device),af_err,())
end

function af_info_string(str,verbose::Bool)
    ccall((:af_info_string,device),af_err,(Ptr{Cstring},Bool),str,verbose)
end

function af_device_info(d_name,d_platform,d_toolkit,d_compute)
    ccall((:af_device_info,device),af_err,(Cstring,Cstring,Cstring,Cstring),d_name,d_platform,d_toolkit,d_compute)
end
=#

function af_get_device_count(num_of_devices::Base.Ref)
    err = ccall((:af_get_device_count, af_lib), Cint, (Ptr{Cint},),num_of_devices)
    err == 0 || throwAFerror(err)
end

#=
function af_get_dbl_support(available,device::Cint)
    ccall((:af_get_dbl_support,device),af_err,(Ptr{Bool},Cint),available,device)
end
=#

function af_set_device(device::Cint)
    err = ccall((:af_set_device, af_lib), Cint, (Cint,), device)
    err == 0 || throwAFerror(err)
end

function af_get_device(device::Base.Ref)
    err = ccall((:af_get_device, af_lib), Cint, (Ptr{Cint},), device)
    err == 0 || throwAFerror(err)
end

function af_sync(device::Int)
    err = ccall((:af_sync, af_lib), Cint, (Cint,), device)
    err == 0 || throwAFerror(err)
end

#=
function af_alloc_device(ptr,bytes::dim_t)
    ccall((:af_alloc_device,device),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes)
end
=#

function af_free_device(ptr::Ptr{Void})
    err = ccall((:af_free_device, af_lib), Cint, (Ptr{Void},), ptr)
    err == 0 || throwAFerror(err)
end

#=
function af_alloc_pinned(ptr,bytes::dim_t)
    ccall((:af_alloc_pinned,device),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes)
end

function af_free_pinned(ptr)
    ccall((:af_free_pinned,device),af_err,(Ptr{Void},),ptr)
end

function af_alloc_host(ptr,bytes::dim_t)
    ccall((:af_alloc_host,device),af_err,(Ptr{Ptr{Void}},dim_t),ptr,bytes)
end

function af_free_host(ptr)
    ccall((:af_free_host,device),af_err,(Ptr{Void},),ptr)
end

function af_device_array(arr,data,ndims::UInt32,dims,_type::af_dtype)
    ccall((:af_device_array,device),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),arr,data,ndims,dims,_type)
end
=#

function af_device_mem_info(alloc_bytes::Base.Ref, alloc_buffers::Base.Ref, lock_bytes::Base.Ref, lock_buffers::Base.Ref)
    err = ccall((:af_device_mem_info, af_lib), Cint,
                (Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t},Ptr{Csize_t}),
                alloc_bytes, alloc_buffers, lock_bytes, lock_buffers)
    err == 0 || throwAFerror(err)
end

#=
function af_print_mem_info(msg,device_id::Cint)
    ccall((:af_print_mem_info,device),af_err,(Cstring,Cint),msg,device_id)
end
=#

function af_device_gc()
    err = ccall((:af_device_gc, af_lib), Cint, ())
    err == 0 || throwAFerror(err)
end

#=
function af_set_mem_step_size(step_bytes::Cint)
    ccall((:af_set_mem_step_size,device),af_err,(Cint,),step_bytes)
end

function af_get_mem_step_size(step_bytes)
    ccall((:af_get_mem_step_size,device),af_err,(Ptr{Cint},),step_bytes)
end

function af_lock_device_ptr(arr::af_array)
    ccall((:af_lock_device_ptr,device),af_err,(af_array,),arr)
end

function af_unlock_device_ptr(arr::af_array)
    ccall((:af_unlock_device_ptr,device),af_err,(af_array,),arr)
end

function af_lock_array(arr::af_array)
    ccall((:af_lock_array,device),af_err,(af_array,),arr)
end

function af_unlock_array(arr::af_array)
    ccall((:af_unlock_array,device),af_err,(af_array,),arr)
end
=#

function af_get_device_ptr(ptr::Base.Ref, arr::AFArray)
    err = ccall((:af_get_device_ptr, af_lib),
                    Cint, (Ptr{Void}, Ptr{Void}), ptr, arr.ptr)
    err == 0 || throwAFerror(err)
end

#=
function af_get_last_error(msg,len)
    ccall((:af_get_last_error,exception),Void,(Ptr{Cstring},Ptr{dim_t}),msg,len)
end

function af_err_to_string(err::af_err)
    ccall((:af_err_to_string,exception),Cstring,(af_err,),err)
end

function af_create_features(feat,num::dim_t)
    ccall((:af_create_features,features),af_err,(Ptr{af_features},dim_t),feat,num)
end

function af_retain_features(out,feat::af_features)
    ccall((:af_retain_features,features),af_err,(Ptr{af_features},af_features),out,feat)
end

function af_get_features_num(num,feat::af_features)
    ccall((:af_get_features_num,features),af_err,(Ptr{dim_t},af_features),num,feat)
end

function af_get_features_xpos(out,feat::af_features)
    ccall((:af_get_features_xpos,features),af_err,(Ptr{af_array},af_features),out,feat)
end

function af_get_features_ypos(out,feat::af_features)
    ccall((:af_get_features_ypos,features),af_err,(Ptr{af_array},af_features),out,feat)
end

function af_get_features_score(score,feat::af_features)
    ccall((:af_get_features_score,features),af_err,(Ptr{af_array},af_features),score,feat)
end

function af_get_features_orientation(orientation,feat::af_features)
    ccall((:af_get_features_orientation,features),af_err,(Ptr{af_array},af_features),orientation,feat)
end

function af_get_features_size(size,feat::af_features)
    ccall((:af_get_features_size,features),af_err,(Ptr{af_array},af_features),size,feat)
end

function af_release_features(feat::af_features)
    ccall((:af_release_features,features),af_err,(af_features,),feat)
end

function af_create_window(out,width::Cint,height::Cint,title)
    ccall((:af_create_window,graphics),af_err,(Ptr{af_window},Cint,Cint,Cstring),out,width,height,title)
end

function af_set_position(wind::af_window,x::UInt32,y::UInt32)
    ccall((:af_set_position,graphics),af_err,(af_window,UInt32,UInt32),wind,x,y)
end

function af_set_title(wind::af_window,title)
    ccall((:af_set_title,graphics),af_err,(af_window,Cstring),wind,title)
end

function af_set_size(wind::af_window,w::UInt32,h::UInt32)
    ccall((:af_set_size,graphics),af_err,(af_window,UInt32,UInt32),wind,w,h)
end

function af_draw_image(wind::af_window,_in::af_array,props)
    ccall((:af_draw_image,graphics),af_err,(af_window,af_array,Ptr{af_cell}),wind,_in,props)
end

function af_draw_plot(wind::af_window,X::af_array,Y::af_array,props)
    ccall((:af_draw_plot,graphics),af_err,(af_window,af_array,af_array,Ptr{af_cell}),wind,X,Y,props)
end

function af_draw_scatter(wind::af_window,X::af_array,Y::af_array,marker::af_marker_type,props)
    ccall((:af_draw_scatter,graphics),af_err,(af_window,af_array,af_array,af_marker_type,Ptr{af_cell}),wind,X,Y,marker,props)
end

function af_draw_scatter3(wind::af_window,P::af_array,marker::af_marker_type,props)
    ccall((:af_draw_scatter3,graphics),af_err,(af_window,af_array,af_marker_type,Ptr{af_cell}),wind,P,marker,props)
end

function af_draw_plot3(wind::af_window,P::af_array,props)
    ccall((:af_draw_plot3,graphics),af_err,(af_window,af_array,Ptr{af_cell}),wind,P,props)
end

function af_draw_hist(wind::af_window,X::af_array,minval::Cdouble,maxval::Cdouble,props)
    ccall((:af_draw_hist,graphics),af_err,(af_window,af_array,Cdouble,Cdouble,Ptr{af_cell}),wind,X,minval,maxval,props)
end

function af_draw_surface(wind::af_window,xVals::af_array,yVals::af_array,S::af_array,props)
    ccall((:af_draw_surface,graphics),af_err,(af_window,af_array,af_array,af_array,Ptr{af_cell}),wind,xVals,yVals,S,props)
end

function af_grid(wind::af_window,rows::Cint,cols::Cint)
    ccall((:af_grid,graphics),af_err,(af_window,Cint,Cint),wind,rows,cols)
end

function af_show(wind::af_window)
    ccall((:af_show,graphics),af_err,(af_window,),wind)
end

function af_is_window_closed(out,wind::af_window)
    ccall((:af_is_window_closed,graphics),af_err,(Ptr{Bool},af_window),out,wind)
end

function af_destroy_window(wind::af_window)
    ccall((:af_destroy_window,graphics),af_err,(af_window,),wind)
end

function af_gradient(dx,dy,_in::af_array)
    ccall((:af_gradient,image),af_err,(Ptr{af_array},Ptr{af_array},af_array),dx,dy,_in)
end
=#

function af_load_image(out::Base.Ref,filename::AbstractString, isColor::Bool)
    err = ccall((:af_load_image, af_lib), Cint,
                (Ptr{Void}, Cstring, Bool), out, filename, isColor)
    err == 0 || throwAFerror(err)
end

function af_save_image(filename::AbstractString, _in::AFArray)
    err = ccall((:af_save_image, af_lib), Cint, 
                (Cstring, Ptr{Void}), filename, _in.ptr)
    err == 0 || throwAFerror(err)
end

#=
function af_load_image_memory(out,ptr)
    ccall((:af_load_image_memory,image),af_err,(Ptr{af_array},Ptr{Void}),out,ptr)
end

function af_save_image_memory(ptr,_in::af_array,format::af_image_format)
    ccall((:af_save_image_memory,image),af_err,(Ptr{Ptr{Void}},af_array,af_image_format),ptr,_in,format)
end

function af_delete_image_memory(ptr)
    ccall((:af_delete_image_memory,image),af_err,(Ptr{Void},),ptr)
end

function af_load_image_native(out,filename)
    ccall((:af_load_image_native,image),af_err,(Ptr{af_array},Cstring),out,filename)
end

function af_save_image_native(filename,_in::af_array)
    ccall((:af_save_image_native,image),af_err,(Cstring,af_array),filename,_in)
end
=#

function af_is_image_io_available(out::Base.Ref)
    err = ccall((:af_is_image_io_available, af_lib), Cint,
                (Ptr{Bool}, ), out)
    err == 0 || throwAFerror(err)
end

function af_resize(out::Base.Ref, _in::AFArray, odim0::Int, odim1::Int, method::Int)
    err = ccall((:af_resize, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Cint, Cint),
                out, _in.ptr, odim0, odim1, method)
    err == 0 || throwAFerror(err)
end

function af_transform(out::Base.Ref, _in::AFArray, transform::AFArray, odim0::Int,
                        odim1::Cint, method::Int, inverse::Bool)
    err = ccall((:af_transform, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint, Cint, Bool),
                out, _in.ptr, transform, odim0, odim1, method, inverse)
    err == 0 || throwAFerror(err)
end

function af_transform_coordinates(out::Base.Ref, tf::AFArray, d0::Real, d1::Real)
    err = ccall((:af_transform_coordinates, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat),
                out, tf.ptr, d0, d1)
    err == 0 || throwAFerror(err)
end

function af_rotate(out::Base.Ref, _in::AFArray, theta::Real, crop::Bool, method::Int)
    err = ccall((:af_rotate, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Bool, Cint),
                out, _in.ptr, theta, crop, method)
    err == 0 || throwAFerror(err)
end

function af_translate(out::Base.Ref, _in::AFArray, trans0::Real, trans1::Real, 
                        odim0::Int, odim1::Int, method::Int)
    err = ccall((:af_translate, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Cint, Cint, Cint),
                out, _in.ptr, trans0, trans1, odim0, odim1, method)
    err == 0 || throwAFerror(err)
end

function af_scale(out::Base.Ref, _in::AFArray, scale0::Real, scale1::Real, 
                    odim0::Int, odim1::Int, method::Int)
    err = ccall((:af_scale, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Cint, Cint, Cint), 
                out, _in.ptr, scale0, scale1, odim0, odim1, method)
    err == 0 || throwAFerror(err)
end

function af_skew(out::Base.Ref, _in::AFArray, skew0::Real, skew1::Real, 
                    odim0::Int, odim1::Int, method::Int, inverse::Bool)
    err = ccall((:af_skew, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Cint, Cint, Cint, Bool),
                out, _in.ptr, skew0, skew1, odim0, odim1, method, inverse)
    err == 0 || throwAFerror(err)
end

function af_histogram(out::Base.Ref, _in::AFArray, nbins::UInt32,minval::Cdouble,maxval::Cdouble)
    err = ccall((:af_histogram, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, UInt32, Cdouble, Cdouble),
                out, _in.ptr, nbins, minval, maxval)
    err == 0 || throwAFerror(err)
end

function af_dilate(out::Base.Ref, _in::AFArray, mask::AFArray)
    err = ccall((:af_dilate, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}), out, _in.ptr, mask)
    err == 0 || throwAFerror(err)
end

function af_dilate3(out::Base.Ref, _in::AFArray, mask::AFArray)
    err = ccall((:af_dilate3, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}), out, _in.ptr, mask)
    err == 0 || throwAFerror(err)
end

function af_erode(out::Base.Ref, _in::AFArray, mask::AFArray)
    err = ccall((:af_erode, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}), out, _in.ptr, mask)
    err == 0 || throwAFerror(err)
end

function af_erode3(out::Base.Ref, _in::AFArray, mask::AFArray)
    err = ccall((:af_erode3, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}), out, _in.ptr, mask)
    err == 0 || throwAFerror(err)
end

function af_bilateral(out::Base.Ref, _in::AFArray, spatial_sigma::Cdouble, chromatic_sigma::Cdouble, isColor::Bool)
    err = ccall((:af_bilateral, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Bool),
                out, _in.ptr, spatial_sigma, chromatic_sigma, isColor)
    err == 0 || throwAFerror(err)
end

function af_mean_shift(out::Base.Ref, _in::AFArray, spatial_sigma::Cfloat, chromatic_sigma::Cfloat, iter::Cuint, is_color::Bool)
    err = ccall((:af_mean_shift, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, UInt32, Bool),
                out, _in.ptr, spatial_sigma, chromatic_sigma, iter, is_color)
    err == 0 || throwAFerror(err)
end

function af_medfilt(out::Base.Ref, _in::AFArray, wind_length::Int, wind_width::Int, edge_pad::Int)
    err = ccall((:af_medfilt, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Cint, Cint), 
                out, _in.ptr, wind_length, wind_width, edge_pad)
    err == 0 || throwAFerror(err)
end

function af_minfilt(out::Base.Ref, _in::AFArray, wind_length::Int, wind_width::Int, edge_pad::Int)
    err = ccall((:af_minfilt, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cint, Cint, Cint),
                out, _in.ptr, wind_length, wind_width, edge_pad)
    err == 0 || throwAFerror(err)
end

function af_maxfilt(out::Base.Ref, _in::AFArray, wind_length::Int, wind_width::Int, edge_pad::Int)
    err = ccall((:af_maxfilt, af_lib), Cint, 
                (Ptr{af_array}, Ptr{Void}, Cint, Cint, Cint),
                out, _in.ptr, wind_length, wind_width, edge_pad)
    err == 0 || throwAFerror(err)
end

function af_regions(out::Base.Ref, _in::AFArray, connectivity::Int, ty::DataType)
    err = ccall((:af_regions, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cint, Cuint),
                out, _in.ptr, Cint, aftype(ty))
    err == 0 || throwAFerror(err)
end

function af_sobel_operator(dx::Base.Ref, dy::Base.Ref, img::AFArray, ker_size::UInt32)
    err = ccall((:af_sobel_operator, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32),
                dx, dy, img.ptr, ker_size)
    err == 0 || throwAFerror(err)
end

function af_rgb2gray(out::Base.Ref, _in::AFArray, rPercent::Cfloat, gPercent::Cfloat, bPercent::Cfloat)
    err = ccall((:af_rgb2gray, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Cfloat),
                out, _in.ptr, rPercent, gPercent, bPercent)
    err == 0 || throwAFerror(err)
end

function af_gray2rgb(out::Base.Ref, _in::AFArray, rFactor::Cdouble, gFactor::Cdouble, bFactor::Cdouble)
    err = ccall((:af_gray2rgb, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cfloat, Cfloat, Cfloat),
                out, _in.ptr, rFactor, gFactor, bFactor)
    err == 0 || throwAFerror(err)
end

function af_hist_equal(out::Base.Ref, _in::AFArray, hist::AFArray)
    err = ccall((:af_hist_equal, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}),
                out, _in.ptr, hist.ptr)
    err == 0 || throwAFerror(err)
end

function af_gaussian_kernel(out::Base.Ref, rows::Int, cols::Int, sigma_r::Real, sigma_c::Real)
    err = ccall((:af_gaussian_kernel, af_lib), Cint,
                (Ptr{Void}, Cint, Cint, Cdouble, Cdouble),
                out, rows, cols, sigma_r, sigma_c)
    err == 0 || throwAFerror(err)
end

function af_hsv2rgb(out::Base.Ref, _in::AFArray)
    err = ccall((:af_hsv2rgb, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_rgb2hsv(out::Base.Ref, _in::AFArray)
    err = ccall((:af_rgb2hsv, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}), out, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_color_space(out::Base.Ref, image::AFArray, to::Int, from::Int)
    err = ccall((:af_color_space, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint, Cint), 
                out, image.ptr, to, from)
    err == 0 || throwAFerror(err)
end

#=
function af_unwrap(out,_in::af_array,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    ccall((:af_unwrap,image),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in,wx,wy,sx,sy,px,py,is_column)
end

function af_wrap(out,_in::af_array,ox::dim_t,oy::dim_t,wx::dim_t,wy::dim_t,sx::dim_t,sy::dim_t,px::dim_t,py::dim_t,is_column::Bool)
    ccall((:af_wrap,image),af_err,(Ptr{af_array},af_array,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,dim_t,Bool),out,_in,ox,oy,wx,wy,sx,sy,px,py,is_column)
end

function af_sat(out,_in::af_array)
    ccall((:af_sat,image),af_err,(Ptr{af_array},af_array),out,_in)
end
=#

function af_ycbcr2rgb(out::Base.Ref, _in::AFArray, standard::Int)
    err = ccall((:af_ycbcr2rgb, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr, standard)
    err == 0 || throwAFerror(err)
end

function af_rgb2ycbcr(out::Base.Ref, _in::AFArray, standard::Int)
    err = ccall((:af_rgb2ycbcr, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr, standard)
    err == 0 || throwAFerror(err)
end

function af_svd(u::Base.Ref, s::Base.Ref, vt::Base.Ref, _in::AFArray)
    err = ccall((:af_svd, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}), 
                u, s, vt, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_svd_inplace(u::Base.Ref, s::Base.Ref, vt::Base.Ref, _in::AFArray)
    err = ccall((:af_svd_inplace, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
                u, s, vt, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_lu(lower::Base.Ref, upper::Base.Ref, pivot::Base.Ref, _in::AFArray)
    err = ccall((:af_lu, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
                lower, upper, pivot, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_lu_inplace(pivot::Base.Ref, _in::AFArray, is_lapack_piv::Bool)
    err = ccall((:af_lu_inplace, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Bool),
                pivot, _in.ptr, is_lapack_piv)
    err == 0 || throwAFerror(err)
end

function af_qr(q::Base.Ref, r::Base.Ref, tau::Base.Ref, _in::AFArray)
    err = ccall((:af_qr, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
                q, r, tau, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_qr_inplace(tau::Base.Ref, _in::AFArray)
    err = ccall((:af_qr_inplace, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}), tau, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_cholesky(out::Base.Ref, info::Base.Ref, _in::AFArray, is_upper::Bool)
    err = ccall((:af_cholesky, af_lib), Cint, 
                (Ptr{Void}, Ptr{Cint}, Ptr{Void}, Bool),    
                out, info, _in.ptr, is_upper)
    err == 0 || throwAFerror(err)
end

function af_cholesky_inplace(info::Base.Ref, _in::AFArray, is_upper::Bool)
    err = ccall((:af_cholesky_inplace, af_lib), Cint,
                (Ptr{Cint}, Ptr{Void}, Bool),
                info, _in.ptr, is_upper)
    err == 0 || throwAFerror(err)
end

function af_solve(x::Base.Ref, a::AFArray, b::AFArray, options::Int)
    err = ccall((:af_solve, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint),
                x, a.ptr, b.ptr, options)
    err == 0 || throwAFerror(err)
end

function af_solve_lu(x::Base.Ref, a::AFArray, piv::AFArray, b::AFArray, options::Int)
    err = ccall((:af_solve_lu, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint),
                x, a.ptr, piv.ptr, b.ptr, options)
    err == 0 || throwAFerror(err)
end

function af_inverse(out::Base.Ref, _in::AFArray, options::Int)
    err = ccall((:af_inverse, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cint), out, _in.ptr, options)
    err == 0 || throwAFerror(err)
end

function af_rank(rank::Base.Ref, _in::AFArray, tol::Cdouble)
    err = ccall((:af_rank, af_lib), Cint, 
                (Ptr{UInt32}, Ptr{Void}, Cdouble), rank, _in.ptr, tol)
    err == 0 || throwAFerror(err)
end

function af_det(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_det, af_lib), Cint, 
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_norm(out::Base.Ref, _in::AFArray, _type::Int, p::Real, q::Real)
    err = ccall((:af_norm, af_lib), Cint, 
                (Ptr{Cdouble}, Ptr{Void}, Cint, Cdouble, Cdouble),
                out, _in.ptr, _type, p, q)
    err == 0 || throwAFerror(err)
end

function af_is_lapack_available(out::Base.Ref)
    err = ccall((:af_is_lapack_available, af_lib), Cint,
                (Ptr{Bool},),out)
    err == 0 || throwAFerror(err)
end

function af_approx1(out::Base.Ref, _in::AFArray, pos::AFArray, method::Int, offGrid::Real)
    err = ccall((:af_approx1, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cfloat),
                out, _in.ptr, pos.ptr, method, offGrid)
    err == 0 || throwAFerror(err)
end

function af_approx2(out::Base.Ref, _in::AFArray, pos0::AFArray, pos1::AFArray, method::Int, offGrid::Real)
    err = ccall((:af_approx2, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cfloat),
                out, _in.ptr, pos0.ptr, pos1.ptr, method, offGrid)
    err == 0 || throwAFerror(err)
end

function af_fft(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int)
    err = ccall((:af_fft, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint),
                out, _in.ptr, norm_factor, odim0)
    err == 0 || throwAFerror(err)
end

function af_fft_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_fft_inplace, af_lib), Cint,
                (Ptr{Void}, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_fft2(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int, odim1::Int)
    err = ccall((:af_fft2, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint, Cint),
                out, _in.ptr, norm_factor, odim0, odim1)
    err == 0 || throwAFerror(err)
end

function af_fft2_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_fft2_inplace, af_lib), Cint, 
            (AFArray, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_fft3(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int, odim1::Int, odim2::Int)
    err = ccall((:af_fft3, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint, Cint, Cint),
                out, _in.ptr, norm_factor, odim0, odim1, odim2)
    err == 0 || throwAFerror(err)
end

function af_fft3_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_fft3_inplace, af_lib), Cint,
                (Ptr{Void}, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_ifft(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int)
    err = ccall((:af_ifft, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint),
                out, _in.ptr, norm_factor, odim0)
    err == 0 || throwAFerror(err)
end

function af_ifft_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_ifft_inplace, af_lib), Cint,
                (Ptr{Void}, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_ifft2(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int, odim1::Int)
    err = ccall((:af_ifft2, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint, Cint), 
                out, _in.ptr, norm_factor, odim0, odim1)
    err == 0 || throwAFerror(err)
end

function af_ifft2_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_ifft2_inplace, af_lib), Cint,
                (Ptr{Void}, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_ifft3(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, odim0::Int, odim1::Int, odim2::Int)
    err = ccall((:af_ifft3, af_lib), Cint, 
                (Ptr{Void},Ptr{Void}, Cdouble, Cint, Cint, Cint),
                out, _in.ptr, norm_factor, odim0, odim1, odim2)
    err == 0 || throwAFerror(err)
end

function af_ifft3_inplace(_in::AFArray, norm_factor::Cdouble)
    err = ccall((:af_ifft3_inplace, af_lib), Cint,
                (Ptr{Void}, Cdouble), _in.ptr, norm_factor)
    err == 0 || throwAFerror(err)
end

function af_fft_r2c(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, pad0::Int)
    err = ccall((:af_fft_r2c, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble, dim_t),
                out, _in.ptr, norm_factor, pad0)
    err == 0 || throwAFerror(err)
end

function af_fft2_r2c(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, pad0::Int, pad1::Int)
    err = ccall((:af_fft2_r2c, af_lib), Cint, 
            (Ptr{Void}, Ptr{Void}, Cdouble, Cint, Cint),
            out, _in.ptr, norm_factor, pad0, pad1)
    err == 0 || throwAFerror(err)
end

function af_fft3_r2c(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, pad0::Int, pad1::Int, pad2::Int)
    errr = ccall((:af_fft3_r2c, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Cdouble, Cint, Cint, Cint),
                out, _in.ptr, norm_factor, pad0, pad1, pad2)
    err == 0 || throwAFerror(err)
end

function af_fft_c2r(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, is_odd::Bool)
    err = ccall((:af_fft_c2r, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cdouble, Bool),
                out, _in.ptr, norm_factor, is_odd)
    err == 0 || throwAFerror(err)
end

function af_fft2_c2r(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, is_odd::Bool)
    err = ccall((:af_fft2_c2r, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cdouble, Bool),
                out, _in.ptr, norm_factor, is_odd)
    err == 0 || throwAFerror(err)
end

function af_fft3_c2r(out::Base.Ref, _in::AFArray, norm_factor::Cdouble, is_odd::Bool)
    err = ccall((:af_fft3_c2r, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cdouble, Bool),
                out, _in.ptr, norm_factor, is_odd)
    err == 0 || throwAFerror(err)
end

function af_convolve1(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int, domain::Int)
    err = ccall((:af_convolve1, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, signal.ptr, filter.ptr, mode, domain)
    err == 0 || throwAFerror(err)
end

function af_convolve2(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int, domain::Int)
    err = ccall((:af_convolve2, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, signal.ptr, filter.ptr, mode, domain)
    err == 0 || throwAFerror(err)
end

function af_convolve3(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int, domain::Int)
    err = ccall((:af_convolve3, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, signal.ptr, filter.ptr, mode, domain)
    err == 0 || throwAFerror(err)
end

#=
function af_convolve2_sep(out,col_filter::af_array,row_filter::af_array,signal::af_array,mode::af_conv_mode)
    ccall((:af_convolve2_sep,signal),af_err,(Ptr{af_array},af_array,af_array,af_array,af_conv_mode),out,col_filter,row_filter,signal,mode)
end
=#

function af_fft_convolve1(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int)
    err = ccall((:af_fft_convolve1, af_lib), Cint,
                (Ptr{Void},Ptr{Void},Ptr{Void}, Cint),
                out, signal.ptr, filter.ptr, mode)
    err == 0 || throwAFerror(err)
end

function af_fft_convolve2(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int)
    err = ccall((:af_fft_convolve2, af_lib), Cint,
                (Ptr{Void},Ptr{Void},Ptr{Void}, Cint),
                out, signal.ptr, filter.ptr, mode)
    err == 0 || throwAFerror(err)
end

function af_fft_convolve3(out::Base.Ref, signal::AFArray, filter::AFArray, mode::Int)
    err = ccall((:af_fft_convolve3, af_lib), Cint,
                (Ptr{Void},Ptr{Void},Ptr{Void}, Cint),
                out, signal.ptr, filter.ptr, mode)
    err == 0 || throwAFerror(err)
end

function af_fir(y::Base.Ref, b::AFArray, x::AFArray)
    err = ccall((:af_fir, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}),
                y, b.ptr, x.ptr)
    err == 0 || throwAFerror(err)
end

function af_iir(y::Base.Ref, b::AFArray, a::AFArray, x::AFArray)
    err = ccall((:af_iir,af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}),
                y, b.ptr, a.ptr, x.ptr)
    err == 0 || throwAFerror(err)
end

function af_mean(out::Base.Ref, _in::AFArray, dim::Cuint)
    err = ccall((:af_mean, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cuint), out, _in.ptr ,dim)
    err == 0 || throwAFerror(err)
end

function af_mean_weighted(out::Base.Ref, _in::AFArray,weights::AFArray, dim::Cuint)
    err = ccall((:af_mean_weighted, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cuint), out, _in.ptr, weights.ptr, dim)
    err == 0 || throwAFerror(err)
end

function af_var(out::Base.Ref, _in::AFArray, isbiased::Bool, dim::Cuint)
    err = ccall((:af_var, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Bool, Cuint), out, _in.ptr , isbiased, dim)
    err == 0 || throwAFerror(err)
end

function af_var_weighted(out::Base.Ref, _in::AFArray, weights::AFArray, dim::Cuint)
    err = ccall((:af_var_weighted, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cuint), out, _in.ptr, weights.ptr, dim)
    err == 0 || throwAFerror(err)
end


function af_stdev(out::Base.Ref, _in::AFArray, dim::Cuint)
    err = ccall((:af_stdev, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cuint), out, _in.ptr ,dim)
    err == 0 || throwAFerror(err)
end

function af_cov(out::Base.Ref, X::AFArray, Y::AFArray, isbiased::Bool)
    err = ccall((:af_cov, af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Bool), 
                out, X.ptr, Y.ptr, isbiased)
    err == 0 || throwAFerror(err)
end

function af_median(out::Base.Ref, _in::AFArray, dim::Cuint)
    err = ccall((:af_median, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cuint), out, _in.ptr ,dim)
    err == 0 || throwAFerror(err)
end

function af_mean_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_mean_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_mean_all_weighted(real::Base.Ref, imag::Base.Ref, _in::AFArray, weights::AFArray)
    err = ccall((:af_mean_all_weighted, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Ptr{Void}), real, imag, _in.ptr, weights)
    err == 0 || throwAFerror(err)
end

function af_var_all(real::Base.Ref, imag::Base.Ref, _in::AFArray, isbiased::Bool)
    err = ccall((:af_var_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Bool), real, imag, _in.ptr, isbiased)
    err == 0 || throwAFerror(err)
end

function af_var_all_weighted(real::Base.Ref,imag::Base.Ref,_in::AFArray,weights::AFArray)
    err = ccall((:af_var_all_weighted, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Ptr{Void}), real, imag, _in.ptr, weights)
    err == 0 || throwAFerror(err)
end

function af_stdev_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_stdev_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_median_all(real::Base.Ref, imag::Base.Ref, _in::AFArray)
    err = ccall((:af_median_all, af_lib), Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}), real, imag, _in.ptr)
    err == 0 || throwAFerror(err)
end

function af_corrcoef(real::Base.Ref, imag::Base.Ref, X::AFArray, Y::AFArray)
    err = ccall((:af_corrcoef,af_lib),Cint,
                (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Void}, Ptr{Void}), real, imag, X.ptr, Y.ptr)
    err == 0 || throwAFerror(err)
end

function af_fast(out::Base.Ref, _in::AFArray, thr::Real, arc_length::UInt32, 
                    non_max::Bool, feature_ratio::Real, edge::UInt32)
    err = ccall((:af_fast, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cfloat, UInt32, Bool, Cfloat, UInt32),
                out, _in.ptr, thr, arc_length, non_max, feature_ratio, edge)
    err == 0 || throwAFerror(err)
end

function af_harris(out::Base.Ref, _in::AFArray, max_corners::UInt32, 
                    min_response::Real, sigma::Real, block_size::UInt32, k_thr::Real)
    err = ccall((:af_harris, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, UInt32, Cfloat, Cfloat, UInt32, Cfloat),
                out, _in.ptr, max_corners, min_response, sigma, block_size, k_thr)
    err == 0 || throwAFerror(err)
end

function af_orb(feat::Base.Ref, desc::Base.Ref, _in::AFArray, 
                fast_thr::Real, max_feat::UInt32, scl_fctr::Real, 
                levels::UInt32, blur_img::Bool)
    err = ccall((:af_orb, af_lib), Cint,
            (Ptr{Void}, Ptr{Void}, Ptr{Void}, Cfloat, UInt32, Cfloat, UInt32, Bool),
            feat, desc, _in.ptr, fast_thr, max_feat, scl_fctr, levels, blur_img)
    err == 0 || throwAFerror(err)
end

function af_sift(feat::Base.Ref, desc::Base.Ref, _in::AFArray, n_layers::UInt32, 
                contrast_thr::Real, edge_thr::Real, init_sigma::Real, double_input::Bool, 
                intensity_scale::Real, feature_ratio::Real)
    err = ccall((:af_sift, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32, Cfloat, 
                Cfloat, Cfloat, Bool, Cfloat, Cfloat),
                feat,desc, _in.ptr, n_layers, contrast_thr, edge_thr, 
                init_sigma, double_input, intensity_scale, feature_ratio)
    err == 0 || throwAFerror(err)
end

function af_gloh(feat::Base.Ref, desc::Base.Ref, _in::AFArray, n_layers::UInt32, 
                contrast_thr::Real, edge_thr::Real, init_sigma::Real, double_input::Bool, 
                intensity_scale::Real, feature_ratio::Real)
    err = ccall((:af_gloh, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, UInt32, Cfloat, 
                Cfloat, Cfloat, Bool, Cfloat, Cfloat),
                feat,desc, _in.ptr, n_layers, contrast_thr, edge_thr, 
                init_sigma, double_input, intensity_scale, feature_ratio)
    err == 0 || throwAFerror(err)
end

function af_hamming_matcher(idx::Base.Ref, dist::Base.Ref, query::AFArray, 
                            train::AFArray, dist_dim::Int, n_dist::UInt32)
    err = ccall((:af_hamming_matcher,vision), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, UInt32),
                idx, dist, query.ptr, train.ptr, dist_dim, n_dist)
    err == 0 || throwAFerror(err)
end

function af_nearest_neighbour(idx::Base.Ref, dist::Base.Ref, query::AFArray, 
                                train::AFArray, dist_dim::Int, n_dist::UInt32,dist_type::Int)
    err = ccall((:af_nearest_neighbour,vision), Cint,
                (Ptr{Void}, Ptr{Void}, Ptr{Void}, Ptr{Void}, Cint, UInt32, Cint),
                idx, dist, query.ptr, train.ptr, dist_dim, n_dist, dist_type)
    err == 0 || throwAFerror(err)
end

#=
function af_match_template(out,search_img::af_array,template_img::af_array,m_type::af_match_type)
    ccall((:af_match_template,vision),af_err,(Ptr{af_array},af_array,af_array,af_match_type),out,search_img,template_img,m_type)
end
=#

function af_susan(out::Base.Ref, _in::AFArray, radius::UInt32, 
                    diff_thr::Real, geom_thr::Real, feature_ratio::Real, edge::UInt32)
    err = ccall((:af_susan,af_lib), Cint,
                (Ptr{Void}, Ptr{Void}, UInt32, Cfloat, Cfloat, Cfloat, UInt32),
                out, _in.ptr, radius, diff_thr, geom_thr, feature_ratio, edge)
    err == 0 || throwAFerror(err)
end

function af_dog(out::Base.Ref, _in::AFArray, radius1::Int, radius2::Int)
    err = ccall((:af_dog, af_lib), Cint, 
                (Ptr{Void}, Ptr{Void}, Cint, Cint),
                out, _in.ptr, radius1, radius2)
    err == 0 || throwAFerror(err)
end

#=
function af_homography(H,inliers,x_src::af_array,y_src::af_array,x_dst::af_array,y_dst::af_array,htype::af_homography_type,inlier_thr::Cfloat,iterations::UInt32,otype::af_dtype)
    ccall((:af_homography,vision),af_err,(Ptr{af_array},Ptr{Cint},af_array,af_array,af_array,af_array,af_homography_type,Cfloat,UInt32,af_dtype),H,inliers,x_src,y_src,x_dst,y_dst,htype,inlier_thr,iterations,otype)
end
=#
