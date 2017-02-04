gc()
device_gc()

start = device_mem_info()[4]

a = AFArray([1 2 3])
b = AFArray([1 2 3])

@scope function sc1(a, b)
    c = a + b + a .* b
    b = a + c + a .* b
    d = c + a.*c + b
    return c, d
end

@scope function sc2(a, b)
    for k = 1:10
        c, b = sc1(a, b)
        afeval(c)
        afeval(b)
    end
    return b
end

b = sc2(a, b)

afeval(b)

ending = device_mem_info()[4]

@test ending - start == 3
