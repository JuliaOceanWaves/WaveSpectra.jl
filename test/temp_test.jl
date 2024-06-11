using WaveSpectra

using Plots
using Unitful: Hz, s, m
using DimensionfulAngles: radᵃ as rad, °ᵃ as °


func = x -> (0.1Hz < x < 1.5Hz) ? 1m^2/Hz : 0m^2/Hz

S = OmnidirectionalSpectrum(func, typeof(1.0Hz))

plot(S)
plot(S, 0Hz, 3Hz; legend=false)






# f = [0.1, 0.2, 0.3, 0.4, 0.5]
# s = [0, 1, 1.5, 1, 0.5]

# S1 = OmnidirectionalSpectrum(f, s)
# S2 = OmnidirectionalSpectrum(f*u"Hz", s*u"m^2/Hz")

# func = x -> (0.1<x<1.5) ? 1 : 0
# u
# S3 = OmnidirectionalSpectrum(func; dims_frequency=Unitful.NoDims, dims_value=Unitful.NoDims)
# S4 = OmnidirectionalSpectrum(ufunc)


# m = vcat(s',s',s')
# ss = s*u"m^2/Hz"
# ff = f*u"Hz"
# mm = vcat(ss',ss',ss')
