
# units
function unit(x::Spectrum, quantity::Symbol)::Units
    ux, u1, u2 = unit(eltype(x)), unit(eltype(x.axis1)), unit(eltype(x.axis2))
    (quantity == :axis1) && return u1
    (quantity == :axis2) && return u2
    (quantity == x.axesnames[1]) && return u1
    (quantity == x.axesnames[2]) && return u2
    (quantity == :integral) && return ux * u1 * u2
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::Spectrum) = unit(x, :spectrum)

function unit(x::OmnidirectionalSpectrum, quantity::Symbol)::Units
    ux, ua = unit(eltype(x)), unit(eltype(x.axis))
    (quantity == :axis) && return ua
    (quantity == axestypes(x.axis)) && return ua
    (quantity == :integral) && return ux * ua
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::OmnidirectionalSpectrum) = unit(x, :spectrum)


# Axes
AxisArrays.axes(x::Spectrum) = (x.axis1, x.axis2)
AxisArrays.axes(x::OmnidirectionalSpectrum) = (x.axis,)

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

axestypes(dim::Dimensions) = Dict(v[2] => k for (k, v) in axesinfo())[dim]
function axestypes(domain::Symbol, quantity::Symbol)
    return Dict(v[1] => k for (k, v) in axesinfo())[(domain, quantity)]
end
axestypes(x::Union{Quantity,Units}) = axestypes(dimension(x))
axestypes(x::AbstractArray{<:Quantity}) = axestypes(dimension(eltype(x)))
axestypes(x::Spectrum) = x.axestypes
axestypes(x::OmnidirectionalSpectrum) = x.axistype

istemporal(x) = (axesinfo()[axestypes(x)][1][1] == :temporal)
isspatial(x) = (axesinfo()[axestypes(x)][1][1] == :spatial)
isdirection(x) = (axesinfo()[axestypes(x)][1][1] == :direction)
ispolar(x::Spectrum) = (x.coordinates == :polar)
iscartesian(x::Spectrum) = (x.coordinates == :cartesian)

coordinates(x::Spectrum) = x.coordinates
axesnames(x::Spectrum) = x.axesnames
axesnames(x::OmnidirectionalSpectrum) = x.axisname

function isevenlyspaced(x::AbstractVector)
    x_range = _convert_to_range(x)
    return ((length(x) == length(x_range)) && isapprox(x, x_range))
end

isevenlyspaced(x::OmnidirectionalSpectrum) = isevenlyspaced(x.axis)
isevenlyspaced(x::Spectrum) = (isevenlyspaced(x.axis1) && isevenlyspaced(x.axis2))


# convert to/from AxisArrays
function AxisArray(x::Spectrum)
    axis1 = Axis{x.axesnames[1]}(x.axis1)
    axis2 = Axis{x.axesnames[2]}(x.axis2)
    return AxisArray(x, axis1, axis2)
end

function Spectrum(x::AxisArray)
    axes = axisvalues(x)
    if length(axes) ≠ 2
        throw(DimensionMismatch("Exactly two axes are required for 'Spectrum'."))
    end
    axis1, axis2 = axes
    return Spectrum(x.data, axis1, axis2)
end

function AxisArray(x::OmnidirectionalSpectrum)
    name = x.axisname
    axis = Axis{name}(x.axis)
    return AxisArray(x, axis)
end

function OmnidirectionalSpectrum(x::AxisArray)
    axis = axisvalues(x)
    if length(axis) ≠ 1
        throw(DimensionMismatch(
            "Exactly one axis is required for 'OmnidirectionalSpectrum'.")
        )
    end
    return OmnidirectionalSpectrum(x.data, axis[1])
end

function DataFrame(x::Spectrum)
    n1 = length(x.axis1)
    n2 = length(x.axis2)
    axis1 = repeat(x.axis1, outer=n2)
    axis2 = repeat(x.axis2, inner=n1)
    return DataFrame(
        x.axesnames[1] => axis1,
        x.axesnames[2] => axis2,
        :spectrum => vec(x.data),
    )
end

function Spectrum(df::DataFrame)
    cols = names(df)
    spec_col = :spectrum
    axis_cols = (spec_col in cols) ? filter(!=(spec_col), cols) : cols[1:2]
    data_col = (spec_col in cols) ? spec_col : cols[3]
    axis1 = unique(df[!, axis_cols[1]])
    axis2 = unique(df[!, axis_cols[2]])
    data = reshape(df[!, data_col], length(axis1), length(axis2))
    return Spectrum(data, axis1, axis2)
end

function DataFrame(x::OmnidirectionalSpectrum)
    return DataFrame(x.axisname => x.axis, :spectrum => x.data)
end

function OmnidirectionalSpectrum(df::DataFrame)
    cols = names(df)
    spec_col = :spectrum
    axis_col = (spec_col in cols) ? first(filter(!=(spec_col), cols)) : cols[1]
    data_col = (spec_col in cols) ? spec_col : cols[2]
    return OmnidirectionalSpectrum(df[!, data_col], df[!, axis_col])
end
