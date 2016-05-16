# Unified API Functions

export getActiveBackend, getBackendCount, AF_BACKEND_DEFAULT,
        AF_BACKEND_CPU, AF_BACKEND_CUDA, AF_BACKEND_OPENCL, 
        getAvailableBackends, setBackend, getBackendId

AF_BACKEND_DEFAULT = UInt32(0)
AF_BACKEND_CPU = UInt32(1)
AF_BACKEND_CUDA = UInt32(2)
AF_BACKEND_OPENCL = UInt32(4)

function getActiveBackend()
    backend = Base.Ref{Cuint}(0)
    af_get_active_backend(backend)
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

# There isn't a af_get_device in the .so

#function getDeviceId(a::AFArray)
#    device = Base.Ref{Cint}(0)
#    af_get_device_id(device, a)
#    Int(device[])
#end
