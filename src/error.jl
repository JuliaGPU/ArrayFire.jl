import Base: showerror

const AF_SUCCESS = 0
const AF_ERR_NO_MEM = 1
const AF_ERR_DRIVER = 2
const AF_ERR_RUNTIME = 3
const AF_ERR_INVALID_ARRAY = 4
const AF_ERR_ARG  = 5
const AF_ERR_SIZE = 6
const AF_ERR_TYPE = 7
const AF_ERR_DIFF_TYPE = 8
const AF_ERR_BATCH = 9
const AF_ERR_DEVICE = 10
const AF_ERR_NOT_SUPPORTED = 11
const AF_ERR_NOT_CONFIGURED  = 12
const AF_ERR_NONFREE = 13
const AF_ERR_NO_DBL = 14
const AF_ERR_NO_GFX = 15
const AF_ERR_LOAD_LIB = 16
const AF_ERR_LOAD_SYM  = 17
const AF_ERR_ARR_BKND_MISMATCH = 18
const AF_ERR_INTERNAL = 19
const AF_ERR_UNKNOWN = 20

export  AF_SUCCESS,
        AF_ERR_NO_MEM,
        AF_ERR_DRIVER,
        AF_ERR_RUNTIME,
        AF_ERR_INVALID_ARRAY,
        AF_ERR_ARG ,
        AF_ERR_SIZE,
        AF_ERR_TYPE,
        AF_ERR_DIFF_TYPE,
        AF_ERR_BATCH,
        AF_ERR_DEVICE ,
        AF_ERR_NOT_SUPPORTED ,
        AF_ERR_NOT_CONFIGURED  ,
        AF_ERR_NONFREE ,
        AF_ERR_NO_DBL ,
        AF_ERR_NO_GFX ,
        AF_ERR_LOAD_LIB ,
        AF_ERR_LOAD_SYM  ,
        AF_ERR_ARR_BKND_MISMATCH ,
        AF_ERR_INTERNAL ,
        AF_ERR_UNKNOWN

function throwAFerror(err::Integer)
    str = ccall((:af_err_to_string, "libaf"), 
                Cstring, (Cint, ), err)
    throw("ArrayFire Error ($err) : $(bytestring(str))")
end
