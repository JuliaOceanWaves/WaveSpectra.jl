# split spectrum into omnidirectional & spread function

function OmnidirectionalSpectrum(x::AbstractSpectrum)
    ispolar(x) && return integrate(x, :direction)
    throw(ArgumentError("Spectrum must be in polar coordinates."))
end

"""
    spread_function(x::AbstractSpectrum)

Return the directional spreading function corresponding to the polar spectrum `x`.
"""
function spread_function(x::AbstractSpectrum)
    return Spectrum(x ./ OmnidirectionalSpectrum(x), AxisArrays.axes(x)...)
end

"""
    split_spectrum(x::AbstractSpectrum)

Split a polar spectrum into an omnidirectional spectrum and a spread function.
"""
split_spectrum(x::AbstractSpectrum) = (OmnidirectionalSpectrum(x), spread_function(x))

function Spectrum(omni::AbstractOmnidirectionalSpectrum, spread::AbstractSpectrum)
    ispolar(spread) || throw(ArgumentError("Spread function must be in polar coordinates."))
    (omni.axis ≉ spread.axis1) && throw(ArgumentError(
        "Frequency axis do not match between omnidirectional spectrum and spread function."
    ))
    return Spectrum(spread.data .* omni.data, AxisArrays.axes(spread)...)
end

"""
    isspread(x::AbstractSpectrum)

Check whether a `Spectrum` is a spread function.
It must be in polar form and integrate to 1 (dimensionless) over the direction axis at every
frequency.
"""
function isspread(x::AbstractSpectrum)
    !ispolar(x) && return false
    return all(integrate(x, :direction) .≈ 1)
end
