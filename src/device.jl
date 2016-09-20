### Device Functions

export  getDevice,
        setDevice,
        getNumDevices

function getDevice()
    n = Base.Ref{Cint}(0)
    af_get_device(n)
    Int(n[])
end

function setDevice(n::Integer)
    af_set_device(Cint(n))
end

function getNumDevices()
    n = Base.Ref{Cint}(0)
    af_get_device_count(n)
    Int(n[])
end
