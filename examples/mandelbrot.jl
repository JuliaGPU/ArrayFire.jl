using ArrayFire

maxIterations = 500
gridSize = 2048
xlim = [-0.748766713922161, -0.748766707771757]
ylim = [ 0.123640844894862,  0.123640851045266]

x = linspace( xlim[1], xlim[2], gridSize )
y = linspace( ylim[1], ylim[2], gridSize )

xGrid = [i for i in x, j in y]
yGrid = [j for i in x, j in y]

z0 = xGrid + im*yGrid

function mandelbrotCPU(z0, maxIterations)
    z = copy(z0)
    count = ones( size(z) )

    for n in 1:maxIterations
        z .= z.*z .+ z0
        count .+= abs.( z ).<=2
    end
    count = log.( count )
end

function mandelbrotGPU(z0, maxIterations)
    z = z0
    count = ones(AFArray{Float32}, size(z) )

    for n in 1:maxIterations
        z = z .* z .+ z0
        count = count + (abs(z)<= 2)
    end
    sync(log( count ))
end

# warmup
count = mandelbrotCPU(z0, 1)

cpu_time = @elapsed count = mandelbrotCPU(z0, maxIterations)

count .-= minimum(count)
count ./= maximum(count)
img = AFArray(Array{Float32}(count))

ArrayFire.image(img)

count = mandelbrotGPU(AFArray(z0), 1)
count = mandelbrotGPU(AFArray(z0), maxIterations)
gpu_time = @elapsed count = mandelbrotGPU(AFArray(z0), maxIterations)

ArrayFire.figure(2)
count -= min_all(count)[1]
count /= max_all(count)[1]
img = AFArray{Float32}(count)

ArrayFire.image(img)


@show cpu_time, gpu_time, cpu_time/gpu_time
