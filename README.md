# ArrayFire

[![Build Status](https://travis-ci.org/Keno/ArrayFire.jl.svg?branch=master)](https://travis-ci.org/Keno/ArrayFire.jl)

ArrayFire is a library for GPU and accelerated computing in Julia. It is a wrapper around [arrayfire](https://github.com/arrayfire/arrayfire), a C++ library, using [Cxx.jl](https://github.com/Keno/Cxx.jl).

##Installation 

First make sure you build Julia with special flags in accordance with Cxx.jl. The instructions for doing that are [here](https://github.com/Keno/Cxx.jl.git). Then build Cxx by doing `Pkg.build("Cxx"`).

Then clone the ArrayFire package by doing:
```julia
Pkg.clone(https://github.com/JuliaComputing/ArrayFire.jl.git)
```
Then, you must build arrayfire for your system. Follow the instructions to build the arrayfire library [here](https://github.com/arrayfire/arrayfire/wiki). Make sure you've installed all the dependencies. If you've chosen to build all three backends, you should see `libafcpu.so`, `libafcuda.so` and `libafopencl.so` in `path_to_arrayfire/build/src/backend/` in the `cpu`, `opencl` and `cuda` folders respectively.

Then, once you have built the library, copy all the build files into the `deps` folder in `~/.julia/v0.5/ArrayFire/deps/build/`

If all goes correctly, `using ArrayFire` should work without errors.

## Usage
ArrayFire creates pointers to GPU memory using the `AFArray` type. Operations on `AFArray` types return `AFArray` types, thereby keeping data on the GPU.

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
## Performance 
ArrayFire was benchmarked on commonly used operations.

![Performance Chart](https://cloud.githubusercontent.com/assets/9101377/13040156/8f33516c-d3ce-11e5-9873-766a2ab67781.png)

CPU: Intel(R) Xeon(R) CPU E5-2670 0 @ 2.60GHz.

GPU: GRID K520, 4096 MB, CUDA Compute 3.0.

ArrayFire v3.3.0

Please contribute to the development of this package by filing issues [here](https://github.com/JuliaComputing/ArrayFire.jl/issues). 
