
# split spectrum into omnidirectional & spread function
function omnidirectional_spectrum(x::Spectrum)
    ispolar(x) && return integrate(x, :direction)
    throw(ArgumentError("Spectrum must be in polar coordinates."))
end

function spread_function(x::Spectrum)
    return Spectrum(x.data ./ omnidirectional_spectrum(x), AxisArrays.axes(x)...)
end

function split_spectrum(x::Spectrum)
    omnidirectional = omnidirectional_spectrum(x)
    spread = Spectrum(x ./ omnidirectional, AxisArrays.axes(x)...)
    return (omnidirectional, spread)
end

# parametric spectra

# spectral moments

# dispersion relations
deepwater = DispersionGradient(;
    dispersion = (k -> √(g * k * θ₀)),
    dispersion_inverse = (ω -> ω^2 / (g * θ₀)),
    gradient=(k -> √(g * θ₀ / k) / 2),
    gradient_inverse = (ω -> 2ω / (g * θ₀)),
)
