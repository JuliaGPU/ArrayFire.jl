
import Base: fft!, ifft!, rfft, irfft
export fft1, fft1!, fft2, fft2!, fft3, fft3!, ifft1, ifft2, ifft3, rfft1, rfft2, rfft3
export irfft1, irfft2, irfft3, ifft1!, ifft2!, ifft3!
export fft_convolve1, fft_convolve2, fft_convolve3
export set_fft_plan_cache_size

function fft(_in::AFArray{T,1},norm_factor=1.0,odim0::dim_t=0) where {T<:Complex}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,1}(out[])
end

function fft(_in::AFArray{T,1},norm_factor=1.0,odim0::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{Complex{T},1}(out[])
end

function fft1(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0) where {T<:Complex,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

function fft1(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{Complex{T},N}(out[])
end

function fft!(_in::AFVector,norm_factor=1.0)
    _error(ccall((:af_fft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function fft1!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function fft(_in::AFArray{T,2},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0) where {T<:Complex}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,2}(out[])
end

function fft(_in::AFArray{T,2},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{Complex{T},2}(out[])
end

function fft2(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0) where {T<:Complex,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

function fft2(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{Complex{T},N}(out[])
end

function fft!(_in::AFMatrix,norm_factor=1.0)
    _error(ccall((:af_fft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function fft2!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function fft(_in::AFArray{T,3},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T<:Complex}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,3}(out[])
end

function fft(_in::AFArray{T,3},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{Complex{T},3}(out[])
end

function fft3(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T<:Complex,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

function fft3(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{Complex{T},N}(out[])
end

function fft!(_in::AFVolume,norm_factor=1.0)
    _error(ccall((:af_fft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function fft3!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft(_in::AFArray{T,1},norm_factor=-1.0,odim0::dim_t=0) where {T}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,1}(out[])
end

function ifft1(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

function ifft!(_in::AFVector,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft1!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft(_in::AFArray{T,2},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0) where {T}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,2}(out[])
end

function ifft2(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

function ifft!(_in::AFMatrix,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft2!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft(_in::AFArray{T,3},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,3}(out[])
end

function ifft3(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

function ifft!(_in::AFVolume,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function ifft3!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
    _in
end

function rfft(_in::AFArray{T,1},norm_factor=1.0,pad0::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),pad0))
    AFArray{Complex{T},1}(out[])
end

function rfft1(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),pad0))
    AFArray{Complex{T},N}(out[])
end

function rfft(_in::AFArray{T,2},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1))
    AFArray{Complex{T},2}(out[])
end

function rfft2(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1))
    AFArray{Complex{T},N}(out[])
end

function rfft(_in::AFArray{T,3},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0,pad2::dim_t=0) where {T<:Real}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1,pad2))
    AFArray{Complex{T},3}(out[])
end

function rfft3(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0,pad2::dim_t=0) where {T<:Real,N}
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1,pad2))
    AFArray{Complex{T},N}(out[])
end

function irfft(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

function irfft2(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

function irfft3(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false) where {T,N}
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

function fft_convolve1(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve1,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

function fft_convolve2(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve2,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

function fft_convolve3(signal::AFArray,filter::AFArray,mode::af_conv_mode)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_convolve3,af_lib),af_err,(Ptr{af_array},af_array,af_array,af_conv_mode),out,signal.arr,filter.arr,mode))
    AFArray!(out[])
end

function set_fft_plan_cache_size(cache_size)
    _error(ccall((:af_set_fft_plan_cache_size,af_lib),af_err,(Csize_t,),Csize_t(cache_size)))
end
