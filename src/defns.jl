const f32 = UInt32(0) # 32-bit floating point values
const c32 = UInt32(1) # 32-bit complex floating point values
const f64 = UInt32(2) # 64-bit complex floating point values
const c64 = UInt32(3) # 64-bit complex floating point values
const b8  = UInt32(4) #  8-bit boolean values
const s32 = UInt32(5) # 32-bit signed integral values
const u32 = UInt32(6) # 32-bit unsigned integral values
const u8  = UInt32(7) #  8-bit unsigned integral values
const s64 = UInt32(8) # 64-bit signed integral values
const u64 = UInt32(9) # 64-bit unsigned integral values

aftype(::Type{Float32})          = f32
aftype(::Type{Complex{Float32}}) = c32
aftype(::Type{Float64})          = f64
aftype(::Type{Complex{Float64}}) = c64
aftype(::Type{Bool})             = b8
aftype(::Type{Int32})            = s32
aftype(::Type{UInt32})           = u32
aftype(::Type{UInt8})            = u8
aftype(::Type{Int64})            = s64
aftype(::Type{UInt64})           = u64

function jltype(af_dtype)
    if af_dtype == f32
        return Float32
    elseif af_dtype == c32
        return Complex{Float32}
    elseif af_dtype == f64
        return Float64
    elseif af_dtype == c64
        return Complex{Float64}
    elseif af_dtype == b8
        return Bool
    elseif af_dtype == s32
        return Int32
    elseif af_dtype == u32
        return UInt32
    elseif af_dtype == u8
        return UInt8
    elseif af_dtype == s64
        return Int64
    elseif af_dtype == u64
        return UInt64
    end
end

function af_promote{T,S}(::Type{T},::Type{S})
    if T == S
        return T
    elseif T == Complex{Float64} || S == Complex{Float64}
        return Complex{Float64}
    elseif T == Complex{Float32} || S == Complex{Float32}
        (T == Float64 || S == Float64) && return Complex{Float64}
        return Complex{Float32}
    elseif T == Float64 || S == Float64
        return Float64
    elseif T == Float32 || S == Float32
        return Float32
    elseif T == UInt64 || S == UInt64
        return UInt64
    elseif T == Int64 || S == Int64
        return Int64
    elseif T == UInt32 || S == UInt32
        return UInt32
    elseif T == Int32 || S == Int32
        return Int32
    elseif T == UInt8 || S == UInt8
        return UInt8
    elseif T == Bool || S == Bool
        return Bool
    else
        return Float32
    end
end


typealias dim Int
#=typealias dim4 Union{Tuple{dim}, Tuple{dim,dim}, Tuple{dim,dim,dim}, Tuple{dim,dim,dim,dim}}

function dims_to_dim4(dims)
	if length(dims) == 1
        (dims[1],)
    elseif length(dims) == 2
        (dims[1], dims[2])
    elseif length(dims) == 3
        (dims[1], dims[2], dims[3])
    elseif length(dims) == 4
        (dims[1], dims[2], dims[3], dims[4])
    else
        throw(ArgumentError("Too many dimensions"))
    end
end=#

