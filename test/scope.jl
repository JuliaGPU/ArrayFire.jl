gc()
device_gc()

start = device_mem_info()[4]

a = AFArray([1 2 3])
b = AFArray([1 2 3])

function sc1(a, b)
    afscope(1) do
        c = a + b + a .* b
        return c
    end
end

function sc2(a, b)
    afscope(1) do
        for k = 1:10
            b = sc1(a, b)
            afeval(b)
        end
        return b
    end
end

b = sc2(a, b)

afeval(b)

ending = device_mem_info()[4]

@test ending - start == 3
