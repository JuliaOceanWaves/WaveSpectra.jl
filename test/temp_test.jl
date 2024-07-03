using WaveSpectra

using Plots
using Unitful: Hz, s, m
using DimensionfulAngles: radᵃ as rad, °ᵃ as °


func = x -> (0.1Hz < x < 1.5Hz) ? 1m^2/Hz : 0m^2/Hz

S = OmnidirectionalSpectrum(func, typeof(1.0Hz))

plot(S)
plot(S, 0Hz, 3Hz; legend=false)

mse(x,y) = .-(x, y) |> z -> .^(z, 2) |> sum |> z -> /(z, 47) |> sqrt




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

# m²Hz⁻¹ = m^2*Hz^-1
# s1 = [0.0m²Hz⁻¹, 0.0m²Hz⁻¹, 0.0m²Hz⁻¹, 0.0m²Hz⁻¹, 0.0m²Hz⁻¹, 2.01m²Hz⁻¹, 9.35m²Hz⁻¹, 23.69m²Hz⁻¹, 21.27m²Hz⁻¹, 15.1m²Hz⁻¹, 35.6m²Hz⁻¹, 35.71m²Hz⁻¹, 18.71m²Hz⁻¹, 14.37m²Hz⁻¹, 8.61m²Hz⁻¹, 7.16m²Hz⁻¹, 5.44m²Hz⁻¹, 4.08m²Hz⁻¹, 1.58m²Hz⁻¹, 2.49m²Hz⁻¹, 1.76m²Hz⁻¹, 1.34m²Hz⁻¹, 1.38m²Hz⁻¹, 0.82m²Hz⁻¹, 0.52m²Hz⁻¹, 0.54m²Hz⁻¹, 0.64m²Hz⁻¹, 0.54m²Hz⁻¹, 0.37m²Hz⁻¹, 0.26m²Hz⁻¹, 0.27m²Hz⁻¹, 0.29m²Hz⁻¹, 0.21m²Hz⁻¹, 0.22m²Hz⁻¹, 0.21m²Hz⁻¹, 0.14m²Hz⁻¹, 0.14m²Hz⁻¹, 0.05m²Hz⁻¹, 0.11m²Hz⁻¹, 0.05m²Hz⁻¹, 0.06m²Hz⁻¹, 0.07m²Hz⁻¹, 0.06m²Hz⁻¹, 0.02m²Hz⁻¹, 0.02m²Hz⁻¹, 0.02m²Hz⁻¹, 0.01m²Hz⁻¹]
# f1 = [0.02f0Hz, 0.0325f0Hz, 0.0375f0Hz, 0.0425f0Hz, 0.0475f0Hz, 0.0525f0Hz, 0.0575f0Hz, 0.0625f0Hz, 0.0675f0Hz, 0.0725f0Hz, 0.0775f0Hz, 0.0825f0Hz, 0.0875f0Hz, 0.0925f0Hz, 0.1f0Hz, 0.11f0Hz, 0.12f0Hz, 0.13f0Hz, 0.14f0Hz, 0.15f0Hz, 0.16f0Hz, 0.17f0Hz, 0.18f0Hz, 0.19f0Hz, 0.2f0Hz, 0.21f0Hz, 0.22f0Hz, 0.23f0Hz, 0.24f0Hz, 0.25f0Hz, 0.26f0Hz, 0.27f0Hz, 0.28f0Hz, 0.29f0Hz, 0.3f0Hz, 0.31f0Hz, 0.32f0Hz, 0.33f0Hz, 0.34f0Hz, 0.35f0Hz, 0.365f0Hz, 0.385f0Hz, 0.405f0Hz, 0.425f0Hz, 0.445f0Hz, 0.465f0Hz, 0.485f0Hz]