# ArrayFire.jl
[![Build Status](https://travis-ci.org/JuliaGPU/ArrayFire.jl.svg)](https://travis-ci.org/JuliaGPU/ArrayFire.jl)
[![codecov](https://codecov.io/gh/JuliaGPU/ArrayFire.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGPU/ArrayFire.jl)


[ArrayFire](http://ArrayFire.com) is a library for GPU and accelerated computing. ArrayFire.jl wraps the ArrayFire library for [Julia](https://JuliaLang.org), and provides a Julia interface.

## Installation

Install ArrayFire library: either [download a binary](http://arrayfire.com/download) from the official site, or you can [build from source](https://github.com/arrayfire/arrayfire).

In Julia 1.0 and up:
```julia
] add ArrayFire
```

## Simple Usage
Congratulations, you've now installed `ArrayFire.jl`! Now what can you do?

Let's say you have a simple Julia array on the CPU:
```julia
a = rand(10, 10)
```

You can transfer this array to the device by calling the `AFArray` constructor on it.
```julia
using ArrayFire  # Don't forget to load the library
ad = AFArray(a)
```

Now let us perform some simple arithmetic on it:
```julia
bd = (ad + 1) / 5
```

Of course, you can do much more than just add and divide numbers. Check the supported functions section for more information.

Now that you're done with all your device computation, you can bring your array back to the CPU (or host):
```julia
b = Array(bd)
```

Here are other examples of simple usage:

```julia
using ArrayFire, LinearAlgebra

# Random number generation
a = rand(AFArray{Float64}, 100, 100)
b = randn(AFArray{Float64}, 100, 100)

# Transfer to device from the CPU
host_to_device = AFArray(rand(100,100))

# Transfer back to CPU
device_to_host = Array(host_to_device)

# Basic arithmetic operations
c = sin(a) + 0.5
d = a * 5

# Logical operations
c = a .> b
any_trues = any(c)

# Reduction operations
total_max = maximum(a)
colwise_min = min(a,2)

# Matrix operations
determinant = det(a)
b_positive = abs(b)
product = a * b
dot_product = a .* b
transposer = a'

# Linear Algebra
lu_fact = lu(a)
cholesky_fact = cholesky(a*a')  # Multiplied to create a positive definite matrix
qr_fact = qr(a)
svd_fact = svd(a)

# FFT
fast_fourier = fft(a)
```

## The Execution Model
`ArrayFire.jl` introduces an `AFArray` type that is a subtype of `AbstractArray`. Operations on `AFArrays` create other `AFArrays`, so data always remains on the device unless it is specifically transferred back. This wrapper provides a simple Julian interface that aims to mimic Base Julia's versatility and ease of use.

**REPL Behaviour**: On the REPL, whenever you create an `AFArray`, the REPL displays the values, just like in Base Julia. This happens because the `showarray` method is overloaded to ensure that every time it is needed to display on the REPL, values are transferred from device to host. This means that every single operation on the REPL involves an implicit memory transfer. This may lead to some slowdown while working interactively depending on the size of the data and memory bandwidth available. You can use a semicolon (`;`) at the end of each statement to disable displaying and avoid that memory transfer. Also, note that in a script, there would be no memory transfer unless a display function is explicitly called (or if you use the `Array` constructor like in the above example).

**Async Behaviour**: `arrayfire` is an asynchronous library. This essentially means that whenever you call a particular function in `ArrayFire.jl`, it would return control to the host almost immediately (which in this case in Julia) and continue executing on the device. This is pretty useful because it would mean that host code that's independent of the device can simply execute while the device computes, resulting in better real world performance.

The library also performs some kernel fusions on elementary arithmetic operations (see the arithmetic section of the Supported Functions). `arrayfire` has an intelligent runtime JIT compliation engine which converts array expressions into the smallest number of OpenCL/CUDA kernels. Kernel fusion not only decreases the number of kernel calls, but also avoids extraneous global memory operations. This asynchronous behaviour ends only when a non-JIT operation is called or an explicit synchronization barrier `sync(array)` is called.

**A note on benchmarking**: In Julia, one would use the `@time` macro to time execution times of functions. However, in this particular case, `@time` would simply time the function call, and the library would execute asynchronously in the background. This would often lead to misleading timings. Therefore, the right way to time individual operations is to run them multiple times, place an explicit synchronization barrier at the end, and take the average of multiple runs.

Also, note that this does not affect how the user writes code. Users can simply write normal Julia code using `ArrayFire.jl` and this asynchronous behaviour is abstracted out. Whenever the data is needed back onto the CPU, an implicit barrier ensures that the computatation is complete, and the values are transferred back.

**Operations between CPU and device arrays**:  Consider the following code. It will return an error:
```julia
a = rand(Float32, 10, 10)
b = AFArray(a)
a - b # Throws Error
```
This is because the two arrays reside in different regions of memory (host and device), and for any coherent operation to be performed, one array would have to be transferred to other region in memory. `ArrayFire.jl` does _not_ do this automatically for performance considerations. Therefore, to make this work, you would have to manually transfer one of the arrays to the other memory. The following operations would work:
```julia
a - Array(b) # Works!
AFArray(a) - b # This works too!
```

**A note on correctness**: Sometimes, `ArrayFire.jl` and Base Julia might return marginally different values from their computation. This is because Julia and `ArrayFire.jl` sometimes use different lower level libraries for BLAS, FFT, etc. For example, Julia uses OpenBLAS for BLAS operations, but `ArrayFire.jl` would use clBLAS for the OpenCL backend and CuBLAS for the CUDA backend, and these libraries might not always the exact same values as OpenBLAS after a certain decimal point. In light of this, users are encouraged to keep testing their codes for correctness.

**A note on performance**: Some operations can be slow due to `Base`'s generic implementations. This is intentional, to enable a "make it work, then make it fast" workflow. When you're ready you can disable slow fallback methods:

```julia
julia> allowslow(AFArray, false)
julia> xs[5]
ERROR: getindex is disabled
```


## Supported Functions

### Creating AFArrays
* `rand`, `randn`, `convert`, `diagm`, `eye`, `range`, `zeros`, `ones`, `trues`, `falses`
* `constant`, `getSeed`, `setSeed`, `iota`

### Arithmetic
* `+`, `-`, `*`, `/`, `^`, `&`, `$`, `|`
* `.+`, `.-`, `.*`, `./`, `.>`, `.>=`, `.<`, `.<=`, `.==`, `.!=, `
* `complex`, `conj`, `real`, `imag`, `max`, `min`, `abs`, `round`, `floor`, `hypot`
* `sigmoid`
* `signbit` (works only in vectorized form on Julia v0.5 - Ref issue #109)

### Linear Algebra
* `cholesky`, `svd`, `lu`, `qr`, `svdfact!`, `lufact!`, `qrfact!`
* `*(matmul)`, `A_mul_Bt`, `At_mul_B`, `At_mul_Bt`, `Ac_mul_B`, `A_mul_Bc`, `Ac_mul_Bc`
* `transpose`, `transpose!`, `ctranspose`, `ctranspose!`
* `det`, `inv`, `rank`, `norm`, `dot`, `diag`, `\`
* `isLAPACKAvailable`, `chol!`, `solveLU`, `upper`, `lower`

### Signal Processing
* `fft`, `ifft`, `fft!`, `ifft!`
* `conv`, `conv2`
* `fftC2R`, `fftR2C`, `conv3`, `convolve`, `fir`, `iir`, `approx1`, `approx2`

### Statistics
* `mean`, `median`, `std`, `var`, `cov`
* `meanWeighted`, `varWeighted`, `corrcoef`

### Vector Algorithms
* `sum`, `min`, `max`, `minimum`, `maximum`, `findmax`, `findmin`
* `countnz`, `any`, `all`, `sort`, `union`, `find`, `cumsum`, `diff`
* `sortIndex`, `sortByKey`, `diff2`, `minidx`, `maxidx`

### Backend Functions
* `get_active_backend`, `get_backend_count`, `get_available_backends`, `set_backend`, `get_backend_id`, `sync`, `get_active_backend_id`

### Device Functions
* `get_device`, `set_device`, `get_device_count`

### Image Processing
* `scale`, `hist`
* `loadImage`, `saveImage`
* `isImageIOAvailable`
* `colorspace`, `gray2rgb`, `rgb2gray`, `rgb2hsv`, `rgb2ycbcr`, `ycbcr2rgb`, `hsv2rgb`
* `regions`, `SAT`
* `bilateral`, `maxfilt`, `meanshift`, `medfilt`, `minfilt`, `sobel`, `histequal`
* `resize`, `rotate`, `skew`, `transform`, `transformCoordinates`, `translate`
* `dilate`, `erode`, `dilate3d`, `erode3d`, `gaussiankernel`

### Computer Vision
* `orb`, `sift`, `gloh`, `diffOfGaussians`, `fast`, `harris`, `susan`, `hammingMatcher`, `nearestNeighbour`, `matchTemplate`

## Performance
ArrayFire was benchmarked on commonly used operations.

<img width="537" alt="general" src="https://cloud.githubusercontent.com/assets/9101377/15921168/36b4f1fc-2e3d-11e6-871a-c8989c5bd279.png">


Another interesting benchmark is [**Non-negative Matrix Factorization**](https://www.wikiwand.com/en/Non-negative_matrix_factorization):

![NMF Benchmark](https://cloud.githubusercontent.com/assets/9101377/15921185/62ad8198-2e3d-11e6-911e-469375a99ecb.png)

CPU: Intel(R) Xeon(R) CPU E5-2670 0 @ 2.60GHz.

GPU: GRID K520, 4096 MB, CUDA Compute 3.0.

ArrayFire v3.4.0

The benchmark scripts are in the benchmark folder, and be run from there by doing by doing:
```julia
include("benchmark.jl")
include("nmf_benchmark.jl")
```

## Backends
There are three backends in `ArrayFire.jl`:
* CUDA Backend
* OpenCL Backend
* CPU Backend

There is yet another backend which essentially allows the user to switch backends at runtime. This is called the unified
backend. `ArrayFire.jl` starts up with the unified backend.

If the backend selected by ArrayFire by default (depends on the available drivers) is not the desired one (depending on the available hardware), you can override the default by setting the environment variable `$JULIA_ARRAYFIRE_BACKEND` before starting Julia (more specifically, before loading the `ArrayFire` module). Possible values for `$JULIA_ARRAYFIRE_BACKEND` are `cpu`, `cuda` and `opencl`.

You may also change the backend at runtime via, e.g., `set_backend(AF_BACKEND_CPU)` (resp. `AF_BACKEND_CUDA` or `AF_BACKEND_OPENCL`).
The unified backend isn't a computational backend by itself but represents an interface to switch
between different backends at runtime. `ArrayFire.jl` starts up with the unified backend, but `get_active_backend()`
will return either a particular default backend, depending on how you have installed the library. For example, if you have
built `ArrayFire.jl` with the CUDA backend, `get_active_backend()` will return `AF_BACKEND_CUDA` backend.


## Troubleshooting
`ArrayFire.jl` isn't working! What do I do?
>Error loading `libaf`

Try adding the path to `libaf` to your `LD_LIBRARY_PATH`.
> `ArrayFire Error (998): Internal Error` whenever you call `rand`

If you're using the CUDA backend, try checking if `libcudart` and `libnvvm` are both in your `LD_LIBRARY_PATH`. This is because `libafcuda` will try to link to these libraries when it loads into Julia. If they're not in your system, install CUDA for your platform.

>`ArrayFire.jl` loads, but `a = rand(AFArray{Float32}, 10)` is stuck.

If you want to use the CUDA backend, check if you have installed CUDA for your platform. If you've installed CUDA, simply downloaded a binary and it still doens't work, try adding `libnvvm`, `libcudart` to your path.

> `ArrayFire.jl` does not work with Atom.

Create a file in your home directory called `.juliarc.jl` and write `ENV["LD_LIBRARY_PATH"] = "/usr/local/lib/"` (or the path to `libaf`) in it. Atom should now be able to load it.

> `ERROR: ArrayFire Error (401) : Double precision not supported for this device`

This error message pops up on devices that do not support double precision: a good example would be the Iris Pro on Macbooks. If you get this message, you should work with single precision. For example, if you're generating random numbers directly on the device, the correct usage in this scenario would be `rand(AFArray{Float32}, 10)` instead of `rand(AFArray{Float64}, 10)`.
