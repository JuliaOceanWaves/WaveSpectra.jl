# split spectrum into omnidirectional & spread function

function OmnidirectionalSpectrum(x::AbstractSpectrum)
    ispolar(x) && return integrate(x, :direction)
    throw(ArgumentError("Spectrum must be in polar coordinates."))
end

"""
    spread_function(x::AbstractSpectrum)

Directional spread function corresponding to the polar spectrum `x`.
"""
function spread_function(x::AbstractSpectrum)
    data = x ./ OmnidirectionalSpectrum(x)
    return _rebuild_spectrum(x, data, AxisArrays.axes(x)...)
end

"""
    split_spectrum(x::AbstractSpectrum)

Split a polar spectrum into an omnidirectional spectrum and a directional spread function.
"""
split_spectrum(x::AbstractSpectrum) = (OmnidirectionalSpectrum(x), spread_function(x))

function Spectrum(omni::AbstractOmnidirectionalSpectrum, spread::AbstractSpectrum)
    ispolar(spread) || throw(ArgumentError("Spread function must be in polar coordinates."))
    (omni.axis ≉ spread.axis1) && throw(ArgumentError(
        "Spectral-variable axis does not match between omnidirectional spectrum and " *
        "spread function."
    ))
    return _rebuild_spectrum(spread, spread.data .* omni.data, AxisArrays.axes(spread)...)
end

"""
    isspread(x::AbstractSpectrum)

Check whether a `Spectrum` is a spread function.

It must be in polar form and integrate to 1 (dimensionless) over the direction axis at every
spectral-variable sample.
"""
function isspread(x::AbstractSpectrum)
    ispolar(x) || return false
    return all(integrate(x, :direction).data .≈ 1)
end
