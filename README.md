# ArrayFire.jl

[![Build Status](https://travis-ci.org/gaika/ArrayFire.jl.svg?branch=master)](https://travis-ci.org/gaika/ArrayFire.jl)
[![codecov](https://codecov.io/gh/gaika/ArrayFire.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gaika/ArrayFire.jl)


[ArrayFire](http://arrayfire.com) is a library for GPU and accelerated computing. ArrayFire.jl wraps the ArrayFire library for Julia, and provides a Julian interface.

## Installation

### OSX

If you are on OSX, the easiest way to install arrayfire is by doing
```
brew install arrayfire
```
This would download and install `arrayfire` and link the libraries `libafcpu.so`, and `libafopencl.so` to your `usr/local/lib/` and link `arrayfire.h` to `/usr/local/include`.

Note that this binary contains libraries only for the CPU (`libafcpu`) and OpenCL backends (`libafopencl`). If you want the CUDA backend, you have to [download a different binary](http://arrayfire.com/login/?redirect_to=http%3A%2F%2Farrayfire.com%2Fdownload), or [build the library from source](https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-OSX).

**NOTE**:
* Even if you do download an `arrayfire` binary with the CUDA backend (`libafcuda`), you need to have CUDA installed on your system. If you don't already, [check out these instructions](http://docs.nvidia.com/cuda/cuda-installation-guide-mac-os-x/index.html#axzz4Axqo0CMQ) on how to install it on a Mac.
* You have to build from source for any custom configurations too (such as linking to a different BLAS library).

### Linux
On Linux, you can either [download a binary](http://arrayfire.com/login/?redirect_to=http%3A%2F%2Farrayfire.com%2Fdownload) from the official site, or you can [build from source](https://github.com/arrayfire/arrayfire/wiki/Build-Instructions-for-Linux).

Now that you have `arrayfire` installed, make sure `libaf` in your system path or `LD_LIBRARY_PATH`. Note that `libaf` is the library for the unified backend. For more information on the unified backend, refer to the backends section.

Now, start Julia, and do:
```julia
Pkg.add("ArrayFire")
```
You can also get the latest nightly version of `ArrayFire.jl` by doing:
```julia
Pkg.checkout("ArrayFire")
```
Check if `ArrayFire.jl` works by running the tests:
```julia
Pkg.test("ArrayFire")
```
If you have any issues getting `ArrayFire.jl` to work, please check the Troubleshooting section below. If it still doesn't work, please file an issue.

### Windows
Just download the [installer](http://arrayfire.com/download/) after creating an account.
Follow the installation steps and make sure that you include the library directory into the PATH variable as advised by the installer.
Now try:
```julia
Pkg.add("ArrayFire")
Pkg.test("ArrayFire")
```

Arrayfire requires vcomp120.dll. If you do not have Visual Studio installed, install the [Visual C++ redistributable](https://www.microsoft.com/en-us/download/details.aspx?id=40784).
## Simple Usage
Congratulations, you've now installed `ArrayFire.jl`! Now what can you do?

Let's say you have a simple Julia array on the CPU:
```julia
a = rand(10, 10)
```
You can transfer this array to the device by calling the `AFArray` constructor on it.
```julia
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
using ArrayFire

#Random number generation
a = rand(AFArray{Float64}, 100, 100)
b = randn(AFArray{Float64}, 100, 100)

#Transfer to device from the CPU
host_to_device = AFArray(rand(100,100))

#Transfer back to CPU
device_to_host = Array(host_to_device)

#Basic arithmetic operations
c = sin(a) + 0.5
d = a * 5

#Logical operations
c = a .> b
any_trues = any(c)

#Reduction operations
total_max = maximum(a)
colwise_min = min(a,2)

#Matrix operations
determinant = det(a)
b_positive = abs(b)
product = a * b
dot_product = a .* b
transposer = a'

#Linear Algebra
lu_fact = lu(a)
cholesky_fact = chol(a*a') #Multiplied to create a positive definite matrix
qr_fact = qr(a)
svd_fact = svd(a)

#FFT
fast_fourier = fft(a)

```
## The Execution Model
`ArrayFire.jl` introduces an `AFArray` type that is a subtype of `AbstractArray`. Operations on `AFArrays` create other `AFArrays`, so data always remains on the device unless it is specifically transferred back. This wrapper provides a simple Julian interface that aims to mimic Base Julia's versatility and ease of use.

**Note on REPL Behaviour**: On the REPL, whenever you create an `AFArray`, the REPL displays the values, just like in Base Julia. This happens because the `showarray` method is overloaded to ensure that every time it is needed to display on the REPL, values are transferred from device to host. This means that every single operation on the REPL involves an implicit memory transfer. This may lead to some slowdown while working interactively depending on the size of the data and memory bandwidth available. You can use a semicolon (`;`) at the end of each statement to disable displaying and avoid that memory transfer. Also, note that in a script, there would be no memory transfer unless a display function is explicitly called (or if you use the `Array` constructor like in the above example).

`arrayfire` is an asynchronous library. This essentially means that whenever you call a particular function in `ArrayFire.jl`, it would return control to the host almost immediately (which in this case in Julia) and continue executing on the device. This is pretty useful because it would mean that host code that's independent of the device can simply execute while the device computes, resulting in better real world performance.

The library also performs some kernel fusions on elementary arithmetic operations (see the arithmetic section of the Supported Functions). `arrayfire` has an intelligent runtime JIT compliation engine which converts array expressions into the smallest number of OpenCL/CUDA kernels. Kernel fusion not only decreases the number of kernel calls, but also avoids extraneous global memory operations. This asynchronous behaviour ends only when a non-JIT operation is called or an explicit synchronization barrier (`sync(array)`) is called.

**A note on benchmarking** : In Julia, one would use the `@time` macro to time execution times of functions. However, in this particular case, `@time` would simply time the function call, and the library would execute asynchronously in the background. This would often lead to misleading timings. Therefore, the right way to time individual operations is to run them multiple times, place an explicit synchronization barrier at the end, and take the average of multiple runs.

Also, note that this doesn't affect how the user writes code. Users can simply write normal Julia code using `ArrayFire.jl` and this asynchronous behaviour is abstracted out. Whenever the data is needed back onto the CPU, an implicit barrier ensures that the computatation is complete, and the values are transferred back.

**A note on operations between CPU and device arrays**:  Consider the following code. It will return an error:
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

## Backends
There are three backends in `ArrayFire.jl`:
* CUDA Backend
* OpenCL Backend
* CPU Backend

There is yet another backend which essentially allows the user to switch backends at runtime. This is called the unified backend. `ArrayFire.jl` starts up with the unified backend. You can switch backends by doing:
```julia
setBackend(AF_BACKEND_CPU)
setBackend(AF_BACKEND_OPENCL)
setBackend(AF_BACKEND_CUDA)
```
You can check which backend you're currently on by doing:
```julia
getActiveBackend()
```

**NOTE**: The unified backend isn't a computational backend by itself but represents an interface to switch between different backends at runtime. `ArrayFire.jl` starts up with the unified backend, but`getActiveBackend()` will return either a particular default backend, depending on how you've installed the library. For example, if you've built `ArrayFire.jl` with the CUDA backend, `getActiveBackend()` will return `CUDA` backend.

## Supported Functions

### Creating AFArrays
* `rand, randn, convert, diagm, eye, range, zeros, ones, trues, falses`
* `constant, getSeed, setSeed, iota`

### Arithmetic
* `+, -, *, /, ^, &, $, | `
* `.+, .-, .*, ./, .>, .>=, .<, .<=, .==, .!=, `
* `complex, conj, real, imag, max, min, abs, round, sign, floor, hypot`
* `sigmoid`

### Linear Algebra
* `chol, svd, lu, qr, lufact!, qrfact!, svdfact!`
* `*(matmul), A_mul_Bt, At_mul_B, At_mul_Bt, Ac_mul_B, A_mul_Bc, Ac_mul_Bc`
* `transpose, transpose!, ctranspose, ctranspose!`
* `det, inv, rank, norm, dot, diag, \`
* `isLAPACKAvailable, chol!, solveLU, upper, lower`

### Signal Processing
* `fft, ifft, fft!, ifft!`
* `conv, conv2`
* `fftC2R, fftR2C, conv3, convolve, fir, iir, approx1, approx2`

### Statistics
* `mean, median, std, var, cov`
* `meanWeighted, varWeighted, corrcoef`

### Vector Algorithms
* `sum, min, max, minimum, maximum, findmax, findmin`
* `countnz, any, all, sort, union, find, cumsum, diff`
* `sortIndex, sortByKey, diff2, minidx, maxidx`

### Backend Functions
* `getActiveBackend, getBackendCount, getAvailableBackends, setBackend, getBackendId, sync, getActiveBackendId`

### Device Functions
* `getDevice`, `setDevice`, `getNumDevices`

### Image Processing
* `scale, hist`
* `loadImage, saveImage`
* `isImageIOAvailable`
* `colorspace, gray2rgb, rgb2gray, rgb2hsv, rgb2ycbcr, ycbcr2rgb, hsv2rgb`
* `regions, SAT`
* `bilateral, maxfilt, meanshift, medfilt, minfilt, sobel, histequal`
* `resize, rotate, skew, transform, transformCoordinates, translate`
* `dilate, erode, dilate3d, erode3d, gaussiankernel`

### Computer Vision
* `orb, sift, gloh, diffOfGaussians, fast, harris, susan, hammingMatcher, nearestNeighbour, matchTemplate`

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

## Backend Selection

When using ArrayFire's [Unified Backend](http://arrayfire.org/docs/unifiedbackend.htm), you can select between the CPU, CUDA and OpenCL backend.

If the backend selected by ArrayFire by default (depends on the available drivers) is not the desired one (depending on the available hardware), you can override the default by setting the environment variable `$JULIA_ARRAYFIRE_BACKEND` before starting Julia (more specifically, before loading the `ArrayFire` module). Possible values for `$JULIA_ARRAYFIRE_BACKEND` are `cpu`, `cuda` and `opencl`.

You may also change the backend at runtime via, e.g., `set_backend(AF_BACKEND_CPU)` (resp. `AF_BACKEND_CUDA` or `AF_BACKEND_OPENCL`).

## Troubleshooting
`ArrayFire.jl` isn't working! What do I do?
>Error loading `libaf`

Try adding the path to `libaf` to your `LD_LIBRARY_PATH`.
> `ArrayFire Error (998): Internal Error` whenever you call `rand`

If you're using the CUDA backend, try checking if `libcudart` and `libnvvm` are both in your `LD_LIBRARY_PATH`. This is because `libafcuda` will try to link to these libraries when it loads into Julia. If they're not in your system, install CUDA for your platform.

>`ArrayFire.jl` loads, but `a = rand(AFArray{Float32}, 10)` is stuck.

If you want to use the CUDA backend, check if you have installed CUDA for your platform. If you've installed CUDA, simply downloaded a binary and it still doens't work, try adding `libnvvm`, `libcudart` to your path.

> `ArrayFire.jl` doesn't work with Atom.

Create a file in your home directory called `.juliarc.jl` and write `ENV["LD_LIBRARY_PATH"] = "/usr/local/lib/"` (or the path to `libaf`) in it. Atom should now be able to load it.

> `ERROR: ArrayFire Error (401) : Double precision not supported for this device`

This error message pops up on devices that do not support double precision: a good example would be the Iris Pro on Macbooks. If you get this message, you should work with single precision. For example, if you're generating random numbers directly on the device, the correct usage in this scenario would be `rand(AFArray{Float32}, 10)` instead of `rand(AFArray{Float64}, 10)`.
