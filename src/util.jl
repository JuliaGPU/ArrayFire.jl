function af_error(err::af_err)
    if err == 0; return; end
    str = af_err_to_string(err)
    throw(ErrorException("ArrayFire Error ($err) : $(unsafe_string(str))"))
end

af_type(::Type{Float32})          = f32
af_type(::Type{Complex{Float32}}) = c32
af_type(::Type{Float64})          = f64
af_type(::Type{Complex{Float64}}) = c64
af_type(::Type{Bool})             = b8
af_type(::Type{Int32})            = s32
af_type(::Type{UInt32})           = u32
af_type(::Type{UInt8})            = u8
af_type(::Type{Int64})            = s64
af_type(::Type{UInt64})           = u64

af_jltype(::Val{f32}) = Float32
af_jltype(::Val{c32}) = Complex{Float32}
af_jltype(::Val{f64}) = Float64
af_jltype(::Val{c64}) = Complex{Float64}
af_jltype(::Val{b8})  = Bool
af_jltype(::Val{s32}) = Int32
af_jltype(::Val{u32}) = UInt32
af_jltype(::Val{u8})  = UInt8
af_jltype(::Val{s64}) = Int64
af_jltype(::Val{u64}) = UInt64

function af_get_numdims!(arr::af_array)
    result = RefValue{UInt32}(0)
    af_error(ccall((:af_get_numdims,af_lib),af_err,(Ptr{UInt32},af_array),result,arr))
    Int(result[])
end

function af_get_type!(arr::af_array)
    _type = RefValue{af_dtype}(0)
    af_error(ccall((:af_get_type,af_lib),af_err,(Ptr{af_dtype},af_array),_type,arr))
    af_jltype(Val{_type[]}())
end

function af_create_array{T,N}(data::AbstractArray{T,N})
    arr = RefValue{af_array}(0)
    sz = size(data)
    af_error(ccall((:af_create_array,af_lib),af_err,(Ptr{af_array},Ptr{Void},UInt32,Ptr{dim_t},af_dtype),
                   arr,Ref(data),UInt32(length(sz)),Ref([sz...]),af_type(T)))
    AFArray{T,N}(arr[])
end

AFArray!(arr::af_array) = AFArray{af_get_type!(arr), af_get_numdims!(arr)}(arr)
