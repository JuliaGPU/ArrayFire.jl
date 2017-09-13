global const cwin = Ref{Int}(0)
global const windows = Dict{Int,af_window}()

# export figure, plot, image

function figure(n)
    if !haskey(windows, n)
        windows[n] = create_window(640, 480, "figure $n")
    else
        show(windows[n])
        if is_window_closed(windows[n])
            windows[n] = create_window(640, 480, "figure $n")
        end
    end
    show(windows[n])
    cwin[] = n
end

function plot(x, y)
    if cwin[] == 0
        figure(1)
    end
    draw_plot(windows[cwin[]], x, y, Ref(af_cell(0, 0, pointer(""), AF_COLORMAP_DEFAULT)))
    show(windows[cwin[]])
end

function image(img)
    if cwin[] == 0
        figure(1)
    end
    draw_image(windows[cwin[]], img, Ref(af_cell(0, 0, pointer(""), AF_COLORMAP_SPECTRUM)))
    show(windows[cwin[]])
end

function grid(wind::af_window,rows::Integer,cols::Integer)
    _error(ccall((:af_grid,af_lib),af_err,(af_window,Cint,Cint),wind,Cint(rows),Cint(cols)))
end
