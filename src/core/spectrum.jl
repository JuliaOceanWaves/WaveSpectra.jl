# Define `Spectrum` struct and its default behavior

@doc """
    Spectrum(x::AxisArray)

Construct a `Spectrum` from an `AxisArray`.

The input must have exactly two axes.
""" Spectrum(x::AxisArray)

@doc """
    Spectrum(omni::AbstractOmnidirectionalSpectrum, spread::AbstractSpectrum)

Reconstruct a directional `Spectrum` by combining an omnidirectional spectrum
with a directional spread function.

`spread` must use polar coordinates and its first axis must match `omni.axis`.
""" Spectrum(omni::AbstractOmnidirectionalSpectrum, spread::AbstractSpectrum)

"""
    Spectrum(
        data::AbstractMatrix,
        axis1::AbstractVector{<:Quantity},
        axis2::AbstractVector{<:Quantity}
    )

Directional wave spectrum with physical-unit axes and data.

The `data` matrix must have size `(length(axis1), length(axis2))`.
The axes vectors must have units corresponding to one of 8 acceptable frequency types:
temporal/spatial frequency, period, angular frequency, or angular period.
The coordinate system is inferred from axis types and is either `:cartesian` or `:polar`.
"""
struct Spectrum{
    TDAT,
    TAX1 <: AbstractVector{<:Quantity},
    TAX2 <: AbstractVector{<:Quantity}
} <: AbstractSpectrum{TDAT}
    data::Matrix{TDAT}
    axis1::TAX1 # AbstractVector{TAX1}
    axis2::TAX2 # AbstractVector{TAX2}
    coordinates::Symbol
    axestypes::Tuple{Symbol, Symbol}
    axesnames::Tuple{Symbol, Symbol}

    function Spectrum(
            data::AbstractMatrix{},
            axis1::AbstractVector{<:Quantity},
            axis2::AbstractVector{<:Quantity}
    )
        # perform checks
        @assert(size(data)==(length(axis1), length(axis2)),
            "Data and axes sizes do not match!")
        _check_typeconsistency(data)
        _check_typeconsistency(axis1)
        _check_typeconsistency(axis2)
        data, axis1, axis2 = _ensure_increasing_axes(data, axis1, axis2)

        # determine coordinates type
        cartesian_1 = (istemporal(axis1) || isspatial(axis1))
        cartesian_2 = (istemporal(axis2) || isspatial(axis2))
        if cartesian_1 && cartesian_2
            coordinates = :cartesian
        elseif cartesian_1 && isdirection(axis2)
            coordinates = :polar
        elseif cartesian_2 && isdirection(axis1)
            coordinates = :polar
            @warn "Swapping order of axes to have direction as second axis."
            axis1, axis2 = axis2, axis1
        else
            throw(ArgumentError("Axes must define a cartesian or polar coordinate."))
        end

        # assign axes types and names
        axes_types = (axestypes(axis1), axestypes(axis2))
        if axes_types[1] == axes_types[2]
            name1 = Symbol(string(axes_types[1]) * "_1")
            name2 = Symbol(string(axes_types[1]) * "_2")
            axesnames = (name1, name2)
        else
            axesnames = axes_types
        end

        return new{eltype(data), typeof(axis1), typeof(axis2)}(
            data, axis1, axis2, coordinates, axes_types, axesnames
        )
    end
end

# array interface: behave like a matrix
Base.size(x::AbstractSpectrum) = size(x.data)
Base.eltype(x::AbstractSpectrum) = eltype(x.data)
Base.copy(x::AbstractSpectrum) = Spectrum(copy(x.data), copy(x.axis1), copy(x.axis2))

Base.getindex(x::AbstractSpectrum, i::Int) = getindex(x.data, i)
Base.getindex(x::AbstractSpectrum, I::Vararg{Int, 2}) = getindex(x.data, I...)
Base.setindex!(x::AbstractSpectrum, v, i::Int) = setindex!(x.data, v, i)
Base.setindex!(x::AbstractSpectrum, v, I::Vararg{Int, 2}) = (x.data[I...] = v)

Base.BroadcastStyle(::Type{<:AbstractSpectrum}) = Broadcast.ArrayStyle{AbstractSpectrum}()

function Base.similar(x::AbstractSpectrum, ::Type{S}, dims::Dims) where {S}
    (dims ≠ size(x)) && return similar(x.data, S, dims)
    return Spectrum(similar(x.data, S, dims), x.axis1, x.axis2)
end

function Base.similar(
        bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{AbstractSpectrum}},
        ::Type{S}
) where {S}
    sp = _find_first_in_broadcast(bc.args, AbstractSpectrum)
    sp === nothing && return similar(Array{S}, axes(bc))
    _check_axes_in_broadcast(bc.args, sp)
    shape = Base.to_shape(axes(bc))
    return similar(sp, S, shape)
end

# fancy indexing using AxisArray
function Base.getindex(x::AbstractSpectrum; kwargs...)
    return Spectrum(getindex(AxisArray(x); _update_kwargs(kwargs)...))
end

function Base.setindex!(x::AbstractSpectrum, v; kwargs...)
    y = AxisArray(x)
    setindex!(y, v; _update_kwargs(kwargs)...)
    x = Spectrum(y)
    return nothing
end

function Base.getindex(
        x::AbstractSpectrum,
        i::Vararg{Union{Quantity, ClosedInterval{<:Quantity}}, 2}
)
    kwargs = Dict()
    for (k, v) in zip(x.axesnames, i)
        kwargs[k] = v
    end
    return getindex(x; kwargs...)
end

function Base.setindex!(
        x::AbstractSpectrum,
        v::Any,
        i::Vararg{Union{Quantity, ClosedInterval{<:Quantity}}, 2}
)
    kwargs = Dict()
    for (key, value) in zip(x.axesnames, i)
        kwargs[key] = value
    end
    setindex!(x, v; kwargs...)
    return nothing
end

# units: extend `Unitful.unit` function
"""
    unit(x::AbstractSpectrum, quantity::Symbol)
    unit(x::AbstractSpectrum)

Extend `Unitful.unit` for spectra.
The `quantity` can be `:axis1` (or the axis name), `axis2` (or the axis name), `:integral`,
or `:spectrum`.
These return the units of the axes, the integral quantity, and the spectral
density (unit of integral quantity / units of axes), respectively.
The default is `quantity=:spectrum`.
"""
function unit(x::AbstractSpectrum, quantity::Symbol)::Units
    ux, u1, u2 = unit(eltype(x)), unit(eltype(x.axis1)), unit(eltype(x.axis2))
    (quantity == :axis1) && return u1
    (quantity == :axis2) && return u2
    (quantity == x.axesnames[1]) && return u1
    (quantity == x.axesnames[2]) && return u2
    (quantity == :integral) && return ux * u1 * u2
    (quantity == :spectrum) && return ux
    throw(ArgumentError("Unknown `quantity`."))
end

unit(x::AbstractSpectrum) = unit(x, :spectrum)

# convert to/from AxisArray
function AxisArray(x::AbstractSpectrum)
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
