
#Image
function loadImage(str::String)
	a = pointer(str)
	return AFArray{Float32}(icxx"loadImage($a);")
end

function saveImage(str::String, a::AFArray{Float32})
	b = pointer(str)
	icxx"saveImage($b, $a);"
end

function rotate(a::AFArray{Float32}, theta)
	icxx"rotate($a, $theta);"
end

function scale(a::AFArray{Float32}, scale1, scale2)
	icxx"scale($a, $scale1, $scale2);"
end

function transform(a::AFArray{Float32}, b::AFArray)
	icxx"transform($a, $b);"
end
