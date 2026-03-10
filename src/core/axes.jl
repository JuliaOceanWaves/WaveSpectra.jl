# Axes utilities

# `AxisArrays.axes` function
AxisArrays.axes(x::AbstractSpectrum) = (x.axis1, x.axis2)
AxisArrays.axes(x::AbstractOmnidirectionalSpectrum) = (x.axis,)


# `axesinfo` function
"""
    axesinfo()
    axesinfo(s::Symbol)
    axesinfo(x::AbstractSpectrum)
    axesinfo(x::AbstractOmnidirectionalSpectrum)

Return axis information including domain (spatial or temporal, or direction), type
(e.g. angular frequency, period, etc.), and its physical dimensions.

- `axesinfo()` returns a dictionary with information for all 9 possible axis types.
- `axesinfo(s::Symbol)` returns information for a specific axis type, e.g.
  `axesinfo(:wavenumber)`
- `axesinfo(x)` return the axis information for the axes of x, where the input x is either
  a `Spectrum` or an `OmnidirectionalSpectrum`.
"""
axesinfo() = Dict(
    :direction => ((:direction, :angle), 𝐀),
    :frequency => ((:temporal, :frequency), 𝐓^-1),
    :angular_frequency => ((:temporal, :angular_frequency), 𝐀 * 𝐓^-1),
    :period => ((:temporal, :period), 𝐓),
    :angular_period => ((:temporal, :angular_period), 𝐓 * 𝐀^-1),
    :wavenumber => ((:spatial, :frequency), 𝐋^-1),
    :angular_wavenumber => ((:spatial, :angular_frequency), 𝐀 * 𝐋^-1),
    :wavelength => ((:spatial, :period), 𝐋),
    :angular_wavelength => ((:spatial, :angular_period), 𝐋 * 𝐀^-1)
)
axesinfo(s::Symbol) = axesinfo()[s]
axesinfo(x) = axesinfo.(axestypes(x))


# `axestypes` function
axestypes(dim::Dimensions) = Dict(v[2] => k for (k, v) in axesinfo())[dim]

function axestypes(domain::Symbol, quantity::Symbol)
    return Dict(v[1] => k for (k, v) in axesinfo())[(domain, quantity)]
end

axestypes(x::Union{Quantity,Units}) = axestypes(dimension(x))
axestypes(x::AbstractArray{<:Quantity}) = axestypes(dimension(eltype(x)))
axestypes(x::AbstractSpectrum) = x.axestypes
axestypes(x::AbstractOmnidirectionalSpectrum) = x.axistype


# utility functions to check type of axes
istemporal(x) = (axesinfo()[axestypes(x)][1][1] == :temporal)
isspatial(x) = (axesinfo()[axestypes(x)][1][1] == :spatial)
isdirection(x) = (axesinfo()[axestypes(x)][1][1] == :direction)
ispolar(x::AbstractSpectrum) = (x.coordinates == :polar)
iscartesian(x::AbstractSpectrum) = (x.coordinates == :cartesian)


# alternate ways to obtain axes properties
coordinates(x::AbstractSpectrum) = x.coordinates
axesnames(x::AbstractSpectrum) = x.axesnames
axesnames(x::AbstractOmnidirectionalSpectrum) = x.axisname


# check if axes are evenly spaced & obtain spacing information
"""
    isevenlyspaced(x::AbstractVector)
    isevenlyspaced(x::AbstractOmnidirectionalSpectrum)
    isevenlyspaced(x::AbstractSpectrum)

Return `true` when axis values are evenly spaced.
For `Spectrum` both axes are checked.
"""
function isevenlyspaced(x::AbstractVector)
    x_range = _convert_to_range(x)
    return ((length(x) == length(x_range)) && isapprox(x, x_range))
end

isevenlyspaced(x::AbstractOmnidirectionalSpectrum) = isevenlyspaced(x.axis)
isevenlyspaced(x::AbstractSpectrum) = (isevenlyspaced(x.axis1) && isevenlyspaced(x.axis2))

"""
    evenspacing(x::AbstractSpectrum)
    evenspacing(x::AbstractOmnidirectionalSpectrum)

Return evenly spaced axis parameters as tuples `(start, step, length)`.

For a `Spectrum` two such tuples are returned, one for each axis.
Throws `ArgumentError` if the axis or axes are not evenly spaced.
"""
function evenspacing(x::AbstractSpectrum)
    !isevenlyspaced(x) && throw(ArgumentError("Axes must be evenly spaced."))
    r1 = _convert_to_range(x.axis1)
    r2 = _convert_to_range(x.axis2)
    return (r1.ref.hi, r1.step.hi, r1.len), (r2.ref.hi, r2.step.hi, r2.len)
end

function evenspacing(x::AbstractOmnidirectionalSpectrum)
    !isevenlyspaced(x) && throw(ArgumentError("Axes must be evenly spaced."))
    r = _convert_to_range(x.axis)
    return (r.ref.hi, r.step.hi, r.len)
end
