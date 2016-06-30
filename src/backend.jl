# Unified API Functions

export getActiveBackend, getBackendCount, AF_BACKEND_DEFAULT,
        AF_BACKEND_CPU, AF_BACKEND_CUDA, AF_BACKEND_OPENCL, 
        getAvailableBackends, setBackend, getBackendId, sync,
        getActiveBackendId, getDeviceId, getVersionNumber

AF_BACKEND_DEFAULT = UInt32(0)
AF_BACKEND_CPU = UInt32(1)
AF_BACKEND_CUDA = UInt32(2)
AF_BACKEND_OPENCL = UInt32(4)

# Available only v3.3.0 onwards
function getActiveBackend()
    backend = getActiveBackendId()
    if backend == 0
        println("Default Backend")
    elseif backend == 1
        println("CPU Backend")
    elseif backend == 2
        println("CUDA Backend")
    elseif backend == 4
        println("OpenCL Backend")
    end
end 

# Available only v3.3.0 onwards
function getActiveBackendId()
    maj, min, pat = getVersionNumber()
    if maj <= 3 && min < 3
        throw("This function is supported only on arrayfire v3.3.0 and onward, 
                current version is v$(maj).$(min).$(pat). Please install a newer 
                version of arrayfire.")
    backend = Base.Ref{Cuint}(0)
    af_get_active_backend(backend)
    backend = Int(backend[])
end

function getAvailableBackends()
    backends = Base.Ref{Cint}(0)
    af_get_available_backends(backends)
    b = Int(backends[])
    if b == 0
        println("None")
    elseif b == 1
        println("CPU")
    elseif b == 2
        println("CUDA")
    elseif b == 3
        println("CPU and CUDA")
    elseif b == 4
        println("OpenCL")
    elseif b == 5
        println("CPU and OpenCL")
    elseif b == 6
        println("CUDA and OpenCL")
    elseif b == 7
        println("CPU, CUDA and OpenCL")
    end
end

function getBackendCount()
    b = Base.Ref{Cuint}(0)
    af_get_backend_count(b)
    Int(b[])
end

function setBackend(back::Cuint)
    af_set_backend(back)
end

function getBackendId(a::AFArray)
    backend = Base.Ref{Cuint}(0)
    af_get_backend_id(backend, a)
    backend = backend[]
    if backend == 0
        println("Default Backend")
    elseif backend == 1
        println("CPU Backend")
    elseif backend == 2
        println("CUDA Backend")
    elseif backend == 4
        println("OpenCL Backend")
    end
    backend
end

function sync()
    b = getActiveBackendId()
    af_sync(b)
end

function sync(b::Int)
    af_sync(b)
end

# There is a af_get_device_id only in ArrayFire  > v3.3.x, possibly only in v3.4.0 onwards

function getDeviceId(a::AFArray)
    maj, min, pat = getVersionNumber()
    if maj < 3 || (maj == 3 && min < 4)
        throw("ArrayFire Version is v$maj.$min.$pat, which does not support this function. Install the latest version of ArrayFire.")
    end
    device = Base.Ref{Cint}(0)
    af_get_device_id(device, a)
    Int(device[])
end

function getVersionNumber()
    maj = Base.Ref{Cint}(0)
    min = Base.Ref{Cint}(0)
    pat = Base.Ref{Cint}(0)
    af_get_version(maj, min, pat)
    maj[], min[], pat[]
end

