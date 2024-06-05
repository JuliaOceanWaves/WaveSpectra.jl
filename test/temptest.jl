
using WaveSpectra, Unitful, Plots

f = [0.1, 0.2, 0.3, 0.4, 0.5]
s = [0, 1, 1.5, 1, 0.5]

S1 = OmnidirectionalSpectrum(f, s)
S2 = OmnidirectionalSpectrum(f*u"Hz", s*u"m^2/Hz")

func = x -> (0.1<x<1.5) ? 1 : 0
ufunc = x -> (0.1u"Hz" < x < 1.5u"Hz") ? 1u"m^2/Hz" : 0u"m^2/Hz"
S3 = OmnidirectionalSpectrum(func; dims_frequency=Unitful.NoDims, dims_value=Unitful.NoDims)
S4 = OmnidirectionalSpectrum(ufunc)


m = vcat(s',s',s')
ss = s*u"m^2/Hz"
ff = f*u"Hz"
mm = vcat(ss',ss',ss')
