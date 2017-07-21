
import Base: fft!, ifft!, rfft, irfft
export fft, fft2, fft2_c2r, fft2!, ifft, ifft2, ifft3, rfft2, rfft3
export fft3, irfft3, fft3!, irfft2
export fft_convolve1, fft_convolve2, fft_convolve3, fft!
export set_fft_plan_cache_size, ifft!, ifft2!, ifft3!

function fft{T<:Complex,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

function fft{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{Complex{T},N}(out[])
end

function fft!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function fft2{T<:Complex,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

function fft2{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{Complex{T},N}(out[])
end

function fft2!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function fft3{T<:Complex,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

function fft3{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{Complex{T},N}(out[])
end

function fft3!(_in::AFArray,norm_factor=1.0)
    _error(ccall((:af_fft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function ifft{T,N}(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),odim0))
    AFArray{T,N}(out[])
end

function ifft!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function ifft2{T,N}(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft2,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1))
    AFArray{T,N}(out[])
end

function ifft2!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft2_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function ifft3{T,N}(_in::AFArray{T,N},norm_factor=-1.0,odim0::dim_t=0,odim1::dim_t=0,odim2::dim_t=0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_ifft3,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),odim0,odim1,odim2))
    AFArray{T,N}(out[])
end

function ifft3!(_in::AFArray,norm_factor=-1.0)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    _error(ccall((:af_ifft3_inplace,af_lib),af_err,(af_array,Cdouble),_in.arr,Cdouble(norm_factor)))
end

function rfft{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t),out,_in.arr,Cdouble(norm_factor),pad0))
    AFArray{Complex{T},N}(out[])
end

function rfft2{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1))
    AFArray{Complex{T},N}(out[])
end

function rfft3{T<:Real,N}(_in::AFArray{T,N},norm_factor=1.0,pad0::dim_t=0,pad1::dim_t=0,pad2::dim_t=0)
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft3_r2c,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,dim_t,dim_t,dim_t),out,_in.arr,Cdouble(norm_factor),pad0,pad1,pad2))
    AFArray{Complex{T},N}(out[])
end

function irfft{T,N}(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

function irfft2{T,N}(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false)
    if norm_factor < 0
        norm_factor = 1/length(_in)
    end
    out = RefValue{af_array}(0)
    _error(ccall((:af_fft2_c2r,af_lib),af_err,(Ptr{af_array},af_array,Cdouble,Bool),out,_in.arr,Cdouble(norm_factor),is_odd))
    AFArray{T,N}(out[])
end

function irfft3{T,N}(_in::AFArray{Complex{T},N},norm_factor=-1.0,is_odd::Bool=false)
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

function set_fft_plan_cache_size(cache_size::Csize_t)
    _error(ccall((:af_set_fft_plan_cache_size,af_lib),af_err,(Csize_t,),cache_size))
end
