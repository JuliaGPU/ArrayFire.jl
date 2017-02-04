gc()
device_gc()

start = device_mem_info()[4]

a = AFArray([1 2 3])
b = AFArray([1 2 3])

@afgc function sc1(a, b)
    @afgc c = a + b + a .* b
    @afgc b = a + c + a .* b
    @afgc d = c + a.*c + b
    return c, d
end

@afgc function sc2(a, b)
    for k = 1:10
        @afgc c, b = sc1(a, b)
        afeval(c)
        afeval(b)
    end
    return b
end

b = sc2(a, b)

afeval(b)

ending = device_mem_info()[4]

@test ending - start == 3
