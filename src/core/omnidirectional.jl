# Define `OmnidirectionalSpectrum` struct and its default behavior

@doc """
    OmnidirectionalSpectrum(x::AxisArray)

Convert from an `AxisArray` to an `OmnidirectionalSpectrum`.

The input must have exactly one axis.
""" OmnidirectionalSpectrum(x::AxisArray)

@doc """
    OmnidirectionalSpectrum(x::AbstractSpectrum)

Calculate the omnidirectional spectrum of a given `Spectrum` by integrating over its
directional axis.
""" OmnidirectionalSpectrum(x::AbstractSpectrum)

"""
    OmnidirectionalSpectrum(
        data::AbstractVector,
        axis::AbstractVector{<:Quantity}
    )

One-dimensional omnidirectional wave spectrum with a physical-unit axis.

The inputs `data` and `axis` must have equal length.
"""
struct OmnidirectionalSpectrum{TDAT, TAX} <: AbstractOmnidirectionalSpectrum{TDAT}
    data::Vector{TDAT}
    axis::Vector{TAX}
    axistype::Symbol
    axisname::Symbol

    function OmnidirectionalSpectrum(
            data::AbstractVector{},
            axis::AbstractVector{<:Quantity}
    )
        # checks
        @assert length(data)==length(axis) "Data and axis lengths do not match!"
        _check_typeconsistency(data)
        _check_typeconsistency(axis)
        data, axis = _ensure_increasing_axis(data, axis)
        if !(istemporal(axis) || isspatial(axis))
            throw(ArgumentError("Invalid axis type for an 'OmnidirectionalSpectrum'."))
        end

        # assign axis type and name
        axistype = axisname = axestypes(axis)

        return new{eltype(data), eltype(axis)}(data, axis, axistype, axisname)
    end
end

# array interface: behave like a matrix
Base.size(x::AbstractOmnidirectionalSpectrum) = size(x.data)
Base.eltype(x::AbstractOmnidirectionalSpectrum) = eltype(x.data)

function Base.copy(x::AbstractOmnidirectionalSpectrum)
    return OmnidirectionalSpectrum(copy(x.data), copy(x.axis))
end

Base.getindex(x::AbstractOmnidirectionalSpectrum, i::Int) = getindex(x.data, i)
Base.setindex!(x::AbstractOmnidirectionalSpectrum, v::Any, i::Int) = setindex!(x.data, v, i)

function Base.BroadcastStyle(::Type{<:AbstractOmnidirectionalSpectrum})
    return Broadcast.ArrayStyle{AbstractOmnidirectionalSpectrum}()
end

function Base.similar(x::AbstractOmnidirectionalSpectrum, ::Type{S}, dims::Dims) where {S}
    (dims ≠ size(x)) && return similar(x.data, S, dims)
    return OmnidirectionalSpectrum(similar(x.data, S, dims), x.axis)
end

function Base.similar(
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AbstractOmnidirectionalSpectrum}},
        ::Type{S}
) where {S}
    sp = _find_first_in_broadcast(bc.args, AbstractOmnidirectionalSpectrum)
    sp === nothing && return similar(Array{S}, axes(bc))
    _check_axes_in_broadcast(bc.args, sp)
    shape = Base.to_shape(axes(bc))
    return similar(sp, S, shape)
end

# fancy indexing using AxisArray
function Base.getindex(x::AbstractOmnidirectionalSpectrum; kwargs...)
    return OmnidirectionalSpectrum(getindex(AxisArray(x); _update_kwargs(kwargs)...))
end

function Base.setindex!(x::AbstractOmnidirectionalSpectrum, v::Any; kwargs...)
    y = AxisArray(x)
    setindex!(y, v; _update_kwargs(kwargs)...)
    x = OmnidirectionalSpectrum(y)
    return nothing
end

function Base.getindex(
        x::AbstractOmnidirectionalSpectrum,
        i::Union{Quantity, ClosedInterval{<:Quantity}}
)
    kwargs = Dict(x.axisname => i,)
    return getindex(x; kwargs...)
end

function Base.setindex!(
        x::AbstractOmnidirectionalSpectrum,
        v::Any,
        i::Union{Quantity, ClosedInterval{<:Quantity}}
)
    kwargs = Dict(x.axisname => i,)
    y = AxisArray(x)
    setindex!(y, v; kwargs...)
    x = OmnidirectionalSpectrum(y)
    return nothing
end

# units: extend `Unitful.unit` function
"""
    unit(x::AbstractOmnidirectionalSpectrum, quantity::Symbol)
    unit(x::AbstractOmnidirectionalSpectrum)

Extend `Unitful.unit` for omnidirectional spectra.
The `quantity` can be `:axis` (or the axis name), `:integral`, or `:spectrum`.
These return the units of the frequency axis, the integral quantity, and the spectral
density (unit of integral quantity / unit of axis), respectively.
The default is `quantity=:spectrum`.
"""
function unit(x::AbstractOmnidirectionalSpectrum, quantity::Symbol)::Units
    ux, ua = unit(eltype(x)), unit(eltype(x.axis))
    (quantity == :axis) && return ua
    (quantity == axestypes(x.axis)) && return ua
    (quantity == :integral) && return ux * ua
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::AbstractOmnidirectionalSpectrum) = unit(x, :spectrum)

# convert to/from AxisArray
function AxisArray(x::AbstractOmnidirectionalSpectrum)
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
