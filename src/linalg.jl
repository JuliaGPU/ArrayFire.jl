### Linear Algebra

# Import functions from Base

import Base:dot, 
            transpose, 
            ctranspose, 
            transpose!, 
            ctranspose!,
            det, 
            inv,
            norm, 
            rank, 
            *, 
            A_mul_Bt, 
            At_mul_B, 
            At_mul_Bt, 
            Ac_mul_B, 
            A_mul_Bc, 
            Ac_mul_Bc, 
            chol, lu, 
            lufact!, 
            qr, 
            qrfact!, 
            svd,
            svdfact!, 
            \, 
            diag

# Export Functions

export  isLAPACKAvailable, 
        chol!, 
        solveLU, 
        upper, 
        lower

# Export constants

export  AF_MAT_NONE,
        AF_MAT_TRANS,
        AF_MAT_CTRANS,
        AF_MAT_CONJ ,
        AF_MAT_UPPER,
        AF_MAT_LOWER,
        AF_MAT_DIAG_UNIT,
        AF_MAT_SYM,
        AF_MAT_POSDEF,
        AF_MAT_ORTHOG,
        AF_MAT_TRI_DIAG,
        AF_MAT_BLOCK_DIAG,
        AF_NORM_VECTOR_1,
        AF_NORM_VECTOR_INF,
        AF_NORM_VECTOR_2,
        AF_NORM_VECTOR_P,
        AF_NORM_MATRIX_1,
        AF_NORM_MATRIX_INF,
        AF_NORM_MATRIX_2,
        AF_NORM_MATRIX_L_PQ,
        AF_NORM_EUCLID 

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
    Array(AFVector{af_promote(T,S)}(out[]))[1]
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

function rank(a::AFArray; tol::Cdouble = 1e-5)
    out = Base.Ref{Cuint}(0)
    af_rank(out, a, tol)
    Int(out[])
end

# Matrix Multiply

function *{T,S}(a::AFArray{T}, b::AFArray{S}; lhsprop = AF_MAT_NONE, rhsprop = AF_MAT_NONE)
    out = new_ptr()
    af_matmul(out, a, b, lhsprop, rhsprop)
    AFMatrix{af_promote(T,S)}(out[])
end

A_mul_Bt{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, rhsprop = AF_MAT_TRANS)

At_mul_B{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, lhsprop = AF_MAT_TRANS)

At_mul_Bt{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, lhsprop = AF_MAT_TRANS, rhsprop = AF_MAT_TRANS)

A_mul_Bc{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, rhsprop=AF_MAT_CTRANS)

Ac_mul_B{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, lhsprop=AF_MAT_CTRANS)

Ac_mul_Bc{T,S}(a::AFMatrix{T}, b::AFMatrix{S}) = *(a, b, lhsprop=AF_MAT_CTRANS, rhsprop=AF_MAT_CTRANS)

# LAPACK

function isLAPACKAvailable()
    out = Base.Ref{Bool}(0)
    af_is_lapack_available(out)
    out[]
end

# Matrix Factorizations

function _chol{T}(a::AFMatrix{T}, is_upper::Bool)
    out = new_ptr()
    info = Base.Ref{Cint}(0)
    af_cholesky(out, info, a, is_upper)
    info[] > 0 && throw(PosDefException(info))
    AFArray{T}(out[])
end

function _chol!{T}(a::AFMatrix{T}, is_upper::Bool)
    info = Base.Ref{Cint}(0)
    af_cholesky_inplace(info, a, is_upper)
    info[] > 0 && throw(PosDefException(info))
    a
end

chol(a::AFMatrix, ::Type{Val{:U}}) = _chol(a, true)
chol(a::AFMatrix, ::Type{Val{:L}}) = _chol(a, false)
chol(a::AFMatrix) = chol(a,Val{:U})

chol!(a::AFMatrix, ::Type{Val{:U}}) = _chol!(a, true)
chol!(a::AFMatrix, ::Type{Val{:L}}) = _chol!(a, false)
chol!(a::AFMatrix) = chol!(a,Val{:U})

function lu(a::AFMatrix)
    l = new_ptr()
    u = new_ptr()
    p = new_ptr()
    af_lu(l, u, p, a)
    AFArray{backend_eltype(l[])}(l[]), AFArray{backend_eltype(u[])}(u[]), 
    (AFArray{backend_eltype(p[])}(p[]) + 1)
end

function lufact!(a::AFMatrix)
    p = new_ptr()
    af_lu_inplace(p, a, true)
    a, AFArray{backend_eltype(p[])}(p[])
end

function qr(a::AFMatrix)
    q = new_ptr()
    r = new_ptr()
    tau = new_ptr()
    af_qr(q, r, tau, a)
    AFArray{backend_eltype(q[])}(q[]), AFArray{backend_eltype(r[])}(r[]), 
    AFArray{backend_eltype(tau[])}(tau[])
end

function qrfact!(a::AFMatrix)
    tau = new_ptr()
    af_qr_inplace(tau, a)
    a, AFArray{backend_eltype(tau[])}(tau[])
end

function svd(a::AFArray)
    u = new_ptr()
    s = new_ptr()
    vt = new_ptr()
    af_svd(u, s, vt, a)
    AFArray{backend_eltype(u[])}(u[]),
    AFArray{backend_eltype(s[])}(s[]),
    AFArray{backend_eltype(vt[])}(vt[])
end

function svdfact!(a::AFArray)
    u = new_ptr()
    s = new_ptr()
    vt = new_ptr()
    af_svd_inplace(u, s, vt, a)
    AFArray{backend_eltype(u[])}(u[]),
    AFArray{backend_eltype(s[])}(s[]),
    AFArray{backend_eltype(vt[])}(vt[])
end

function \(a::AFArray, b::AFArray; options = AF_MAT_NONE)
    out = new_ptr()
    af_solve(out, a, b, options)
    AFArray{backend_eltype(out[])}(out[])
end

function solveLU(a::AFArray, p::AFArray, b::AFArray; options = AF_MAT_NONE)
    out = new_ptr()
    af_solve_lu(out, a, p, b, options)
    AFArray{backend_eltype(out[])}(out[])
end

function diag{T}(a::AFMatrix{T}, k::Integer = 1)
    out = new_ptr()
    af_diag_extract(out, a, k - 1)
    AFVector{T}(out[])
end

function upper{T}(a::AFMatrix{T}; is_unit_diag = false)
    out = new_ptr()
    af_upper(out, a, is_unit_diag)
    AFMatrix{T}(out[])
end

function lower{T}(a::AFMatrix{T}; is_unit_diag = false)
    out = new_ptr()
    af_lower(out, a, is_unit_diag)
    AFMatrix{T}(out[])
end
