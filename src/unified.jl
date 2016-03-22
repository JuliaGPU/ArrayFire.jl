#Functions to deal with unified backend

export getActiveBackend, setBackend, getAvailableBackends, getBackendCount, getBackend
export AF_BACKEND_CPU, AF_BACKEND_CUDA, AF_BACKEND_OPENCL

function getActiveBackend()
    backend = af_getActiveBackend()
    if backend == AF_BACKEND_CPU
        println("CPU Backend")
    elseif backend == AF_BACKEND_CUDA
        println("CUDA Backend")
    elseif backend == AF_BACKEND_OPENCL
        println("OPENCL Backend")
    end
end

setBackend(backend::Cxx.CppEnum{:af_backend}) = af_setBackend(backend)

function getAvailableBackends() 
    n = af_getAvailableBackends()
    if n == 0
        println("None")
    elseif n == 1
        println("CPU")
    elseif n == 2
        println("CUDA")
    elseif n == 3
        println("CPU and CUDA")
    elseif n == 4
        println("OpenCL")
    elseif n == 5
        println("CPU and OpenCL")
    elseif n == 6
        println("CUDA and OpenCL")
    elseif n == 7
        println("CPU, CUDA and OpenCL")
    end
end

getBackendCount() = af_getBackendCount()
function getBackend(a::AFAbstractArray) 
    b = af_getBackendId(a)
    if b == AF_BACKEND_CPU
        println("CPU Backend")
    elseif b == AF_BACKEND_OPENCL
        println("OpenCL Backend")
    elseif b == AF_BACKEND_CUDA
        println("CUDA Backend")
    end
end
