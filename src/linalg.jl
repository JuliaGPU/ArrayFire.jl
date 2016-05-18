### Linear Algebra

import Base: dot, transpose, ctranspose, transpose!, ctranspose!, det, inv,
                norm

# Constants

AF_MAT_NONE       = 0
AF_MAT_TRANS      = 1
AF_MAT_CTRANS     = 2
AF_MAT_CONJ       = 4    
AF_MAT_UPPER      = 32   
AF_MAT_LOWER      = 64   
AF_MAT_DIAG_UNIT  = 128  
AF_MAT_SYM        = 512  
AF_MAT_POSDEF     = 1024 
AF_MAT_ORTHOG     = 2048 
AF_MAT_TRI_DIAG   = 4096 
AF_MAT_BLOCK_DIAG = 8192  

AF_NORM_VECTOR_1 = 0      
AF_NORM_VECTOR_INF = 1   
AF_NORM_VECTOR_2 = 2 
AF_NORM_VECTOR_P = 3    
AF_NORM_MATRIX_1 = 4   
AF_NORM_MATRIX_INF = 5   
AF_NORM_MATRIX_2 = 6 
AF_NORM_MATRIX_L_PQ = 7   
AF_NORM_EUCLID = AF_NORM_VECTOR_2

# Dot

function dot{T,S}(a::AFVector{T}, b::AFVector{S}, lhsprop = AF_MAT_NONE, rhsprop = AF_MAT_NONE)
    out = new_ptr()
    af_dot(out, a, b, lhsprop, rhsprop)
    AFVector{af_promote(T,S)}(out[])
end
    
# Transpose

function transpose{T}(a::AFArray{T})
    out = new_ptr()
    af_transpose(out, a, false)
    AFArray{T}(out[])
end

function ctranspose{T}(a::AFArray{T})
    out = new_ptr()
    af_transpose(out, a, true)
    AFArray{T}(out[])
end

function transpose!{T}(a::AFArray{T})
    af_transpose_inplace(a, false)
    a
end

function ctranspose!{T}(a::AFArray{T})
    af_transpose_inplace(a, true)
    a
end
 
# Matrix Operations

function det{T<:Real}(a::AFArray{T})
    real = Base.Ref{Cdouble}(0)
    imag = Base.Ref{Cdouble}(0)
    af_det(real, imag, a)
    real[]
end

function det{T<:Complex}(a::AFArray{T})
    real = Base.Ref{Cdouble}(0)
    imag = Base.Ref{Cdouble}(0)
    af_det(real, imag, a)
    complex(real[], imag[])
end

function inv(a::AFArray, options::Int = AF_MAT_NONE)
    out = new_ptr()
    af_inverse(out, a, options)
    AFArray{backend_eltype(out[])}(out[])
end

function norm(a::AFArray; options = AF_NORM_EUCLID, p::Real = 1, q::Real = 1)
    out = Base.Ref{Cdouble}(0)
    af_norm(out, a, options, p, q)
    out[]
end
