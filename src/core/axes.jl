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

Query the axis type and its physical dimensions.

There are 8 supported spectral-variable types, formed from the combinations of
spatial/temporal domain, linear/angular geometry, and frequency/period quantity.
There's one direction type `:direction`, for a total of 9 possible axis types.

- `axesinfo()` returns a dictionary with information for all 9 possible axis types.
- `axesinfo(s::Symbol)` returns information for a specific axis type, e.g.
  `axesinfo(:wavenumber)`
- `axesinfo(x)` return the axis information for the axes of a spectrum x.
"""
function axesinfo end

const _AXESINFO = Dict(
    :direction => ((:direction,), 𝐀),
    :frequency => ((:temporal, :linear, :frequency), 𝐓^-1),
    :angular_frequency => ((:temporal, :angular, :frequency), 𝐀 * 𝐓^-1),
    :period => ((:temporal, :linear, :period), 𝐓),
    :angular_period => ((:temporal, :angular, :period), 𝐓 * 𝐀^-1),
    :wavenumber => ((:spatial, :linear, :frequency), 𝐋^-1),
    :angular_wavenumber => ((:spatial, :angular, :frequency), 𝐀 * 𝐋^-1),
    :wavelength => ((:spatial, :linear, :period), 𝐋),
    :angular_wavelength => ((:spatial, :angular, :period), 𝐋 * 𝐀^-1)
)

const _AXESTYPES_BY_DIM = Dict(v[2] => k for (k, v) in _AXESINFO)
const _AXESTYPES_BY_INFO = Dict(v[1] => k for (k, v) in _AXESINFO)

axesinfo() = _AXESINFO
axesinfo(s::Symbol) = axesinfo()[s]
axesinfo(x) = axesinfo.(axestypes(x))

# `axestypes` function
axestypes(dim::Dimensions) = _AXESTYPES_BY_DIM[dim]

function axestypes(domain::Symbol, geometry::Symbol, quantity::Symbol)
    _AXESTYPES_BY_INFO[(domain, geometry, quantity)]
end

axestypes(x::Union{Quantity, Units}) = axestypes(dimension(x))
axestypes(x::AbstractVector{<:Quantity}) = axestypes(dimension(eltype(x)))
axestypes(x::AbstractSpectrum) = x.axestypes
axestypes(x::AbstractOmnidirectionalSpectrum) = x.axistype

# utility functions to check type of axes
istemporal(x) = (axesinfo()[axestypes(x)][1][1] == :temporal)
isspatial(x) = (axesinfo()[axestypes(x)][1][1] == :spatial)
islinear(x) = (axesinfo()[axestypes(x)][1][2] == :linear)
isangular(x) = (axesinfo()[axestypes(x)][1][2] == :angular)
isfrequency(x) = (axesinfo()[axestypes(x)][1][3] == :frequency)
isperiod(x) = (axesinfo()[axestypes(x)][1][3] == :period)
isdirection(x) = (axesinfo()[axestypes(x)][1][1] == :direction)
isspectralvariable(x) = !isdirection(x)
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
    evenspacing(x::AbstractVector)
    evenspacing(x::AbstractSpectrum)
    evenspacing(x::AbstractOmnidirectionalSpectrum)

Return evenly spaced axis parameters as tuples `(start, step, length)`.

For a `Spectrum` two such tuples are returned, one for each axis.
Throws `ArgumentError` if the axis or axes are not evenly spaced.
"""
function evenspacing(x::AbstractVector)
    isevenlyspaced(x) || throw(ArgumentError("Vector `x` must be evenly spaced."))
    r = _convert_to_range(x)
    return (r.ref.hi, r.step.hi, r.len)
end

evenspacing(x::AbstractSpectrum) = (evenspacing(x.axis1), evenspacing(x.axis2))
evenspacing(x::AbstractOmnidirectionalSpectrum) = evenspacing(x.axis)
