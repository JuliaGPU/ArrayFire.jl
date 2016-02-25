immutable AFWindow
    win::vcpp"af::Window"
    function AFWindow(w::vcpp"af::Window")
        new(w)
    end
end

show(io::IO, w::AFWindow) = print(io, "ArrayFire Window")
show(io::IO, w::vcpp"af::Window") = print(io, "ArrayFire Window")

Cxx.cppconvert(w::AFWindow) = w.win

export AFWindow, window, setTitle, setPlot, setImage, setSurfacePlot

window() = AFWindow(af_window())
window(str::AbstractString) = AFWindow(af_window(str))
window(width::Integer, height::Integer, str::AbstractString) = AFWindow(af_window(width, height, str))

setTitle(w::AFWindow, str::AbstractString) = af_setTitle(w, str)
setImage(w::AFWindow, img::AFAbstractArray, str = "") = af_setImage(w, img/255.0, str)
setPlot(w::AFWindow, x::AFAbstractArray, y::AFAbstractArray, title = "") = af_setPlot(w, x, y, title)
setSurfacePlot(w::AFWindow, x::AFAbstractArray) = af_setSurface(w, x)
show(w::AFWindow) = af_show(w)
