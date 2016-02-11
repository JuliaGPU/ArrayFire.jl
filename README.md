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
