# 2D directional spectrum
struct Spectrum{
    TDAT,
    TAX1<:AbstractVector{<:Quantity},
    TAX2<:AbstractVector{<:Quantity}
} <: AbstractMatrix{TDAT}

    data :: Matrix{TDAT}
    axis1:: TAX1 # AbstractVector{TAX1}
    axis2:: TAX2 # AbstractVector{TAX2}
    coordinates :: Symbol
    axestypes :: Tuple{Symbol, Symbol}
    axesnames :: Tuple{Symbol, Symbol}

    function Spectrum(
        data::AbstractMatrix{},
        axis1::AbstractVector{<:Quantity},
        axis2::AbstractVector{<:Quantity}
    )
        # perform checks
        @assert(
            size(data) == (length(axis1), length(axis2)),
            "Data and axes sizes do not match!"
        )
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

Base.size(x::Spectrum) = size(x.data)
Base.eltype(x::Spectrum) = eltype(x.data)
Base.copy(x::Spectrum) = Spectrum(copy(x.data), copy(x.axis1), copy(x.axis2))

Base.getindex(x::Spectrum, i::Int) = getindex(x.data, i)
Base.getindex(x::Spectrum, I::Vararg{Int, 2}) = getindex(x.data, I...)
Base.setindex!(x::Spectrum, v, i::Int) = setindex!(x.data, v, i)
Base.setindex!(x::Spectrum, v, I::Vararg{Int, 2}) = (x.data[I...] = v)

function Base.getindex(x::Spectrum; kwargs...)
    return Spectrum(getindex(AxisArray(x); _update_kwargs(kwargs)...))
end

function Base.setindex!(x::Spectrum, v; kwargs...)
    y = AxisArray(x)
    setindex!(y, v; _update_kwargs(kwargs)...)
    x = Spectrum(y)
    return nothing
end

function Base.getindex(
    x::Spectrum,
    i::Vararg{Union{Quantity, ClosedInterval{<:Quantity}}, 2}
)
    kwargs = Dict()
    for (k,v) in zip(x.axesnames, i)
        kwargs[k] = v
    end
    return getindex(x; kwargs...)
end

function Base.setindex!(
    x::Spectrum,
    v::Any,
    i::Vararg{Union{Quantity, ClosedInterval{<:Quantity}}, 2},
)
    kwargs = Dict()
    for (key, value) in zip(x.axesnames, i)
        kwargs[key] = value
    end
    setindex!(x, v; kwargs...)
    return nothing
end

Base.BroadcastStyle(::Type{<:Spectrum}) = Broadcast.ArrayStyle{Spectrum}()

function Base.similar(x::Spectrum, ::Type{S}, dims::Dims) where S
    (dims ≠ size(x)) && return similar(x.data, S, dims)
    return Spectrum(similar(x.data, S, dims), x.axis1, x.axis2)
end

function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{Spectrum}},
    ::Type{S}
) where S
    sp = _find_first_in_broadcast(bc.args, Spectrum)
    sp === nothing && return similar(Array{S}, axes(bc))
    check_axes_in_broadcast(bc.args, sp)
    shape = Base.to_shape(axes(bc))
    return similar(sp, S, shape)
end


# 1D omnidirectional spectrum
struct OmnidirectionalSpectrum{TDAT, TAX}  <: AbstractVector{TDAT}
    data :: Vector{TDAT}
    axis :: Vector{TAX}
    axistype :: Symbol
    axisname :: Symbol

    function OmnidirectionalSpectrum(
        data::AbstractVector{},
        axis::AbstractVector{<:Quantity}
    )
        # checks
        @assert length(data) == length(axis) "Data and axis lengths do not match!"
        _check_typeconsistency(data)
        _check_typeconsistency(axis)
        issorted(axis) || throw(ArgumentError("Axis must be increasing."))
        if !(istemporal(axis) || isspatial(axis))
            throw(ArgumentError("Invalid axis type for an 'OmnidirectionalSpectrum'."))
        end

        # assign axis type and name
        axistype = axisname = axestypes(axis)

        return new{eltype(data), eltype(axis)}(data, axis, axistype, axisname)
    end
end

Base.size(x::OmnidirectionalSpectrum) = size(x.data)
Base.eltype(x::OmnidirectionalSpectrum) = eltype(x.data)
Base.copy(x::OmnidirectionalSpectrum) = OmnidirectionalSpectrum(copy(x.data), copy(x.axis))

Base.getindex(x::OmnidirectionalSpectrum, i::Int) = getindex(x.data, i)
Base.setindex!(x::OmnidirectionalSpectrum, v::Any, i::Int) = setindex!(x.data, v, i)

function Base.getindex(x::OmnidirectionalSpectrum; kwargs...)
    return OmnidirectionalSpectrum(getindex(AxisArray(x); _update_kwargs(kwargs)...))
end

function Base.setindex!(x::OmnidirectionalSpectrum, v::Any; kwargs...)
    y = AxisArray(x)
    setindex!(y, v; _update_kwargs(kwargs)...)
    x = OmnidirectionalSpectrum(y)
    return nothing
end

function Base.getindex(
    x::OmnidirectionalSpectrum,
    i::Union{Quantity, ClosedInterval{<:Quantity}}
)
    kwargs = Dict(x.axisname => i,)
    return getindex(x; kwargs...)
end

function Base.setindex!(
    x::OmnidirectionalSpectrum,
    v::Any,
    i::Union{Quantity, ClosedInterval{<:Quantity}}
)
    kwargs = Dict(x.axisname => i,)
    y = AxisArray(x)
    setindex!(y, v; kwargs...)
    x = OmnidirectionalSpectrum(y)
    return nothing
end

function Base.BroadcastStyle(::Type{<:OmnidirectionalSpectrum})
    return Broadcast.ArrayStyle{OmnidirectionalSpectrum}()
end

function Base.similar(x::OmnidirectionalSpectrum, ::Type{S}, dims::Dims) where S
    (dims ≠ size(x)) && return similar(x.data, S, dims)
    return OmnidirectionalSpectrum(similar(x.data, S, dims), x.axis)
end

function Base.similar(
    bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{OmnidirectionalSpectrum}},
    ::Type{S}
) where S
    sp = _find_first_in_broadcast(bc.args, OmnidirectionalSpectrum)
    sp === nothing && return similar(Array{S}, axes(bc))
    check_axes_in_broadcast(bc.args, sp)
    shape = Base.to_shape(axes(bc))
    return similar(sp, S, shape)
end


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


# Axes information
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
axestypes(x::Union{Quantity, Units}) = axestypes(dimension(x))
axestypes(x::AbstractArray{<:Quantity}) = axestypes(dimension(eltype(x)))
axestypes(x::Spectrum) = x.axestypes
axestypes(x::OmnidirectionalSpectrum) = x.axistype

istemporal(x) = (axesinfo()[axestypes(x)][1][1]==:temporal)
isspatial(x) = (axesinfo()[axestypes(x)][1][1] == :spatial)
isdirection(x) = (axesinfo()[axestypes(x)][1][1] == :direction)

# Support functions
@inline function _check_typeconsistency(x::AbstractArray)::Nothing
    consistent = all(y->typeof(y)==eltype(x), x)
    !consistent && throw(ArgumentError("All elements of array must be of same type."))
    return nothing
end

@inline function _update_kwargs(d::AbstractDict)
    kwargs = Dict()
    for (k, v) in d
        if (typeof(v) <: Integer)
            kwargs[k] = (v:v)
        elseif (typeof(v) <: Quantity)
            kwargs[k] = (v .. v)
        else
            kwargs[k] = v
        end
    end
    return kwargs
end

@inline function _find_first_in_broadcast(args, ::Type{T}) where {T}
    for arg in args
        if arg isa T
            return arg
        elseif arg isa Broadcast.Broadcasted
            sp = _find_first_in_broadcast(arg.args, T)
            sp === nothing || return sp
        end
    end
end

@inline function _ensure_increasing_axes(data, axis1, axis2)
    if issorted(axis1)
        data1 = data
    elseif issorted(axis1; rev=true)
        axis1 = reverse(axis1)
        data1 = reverse(data, dims=1)
    else
        throw(ArgumentError("Axis 1 must be monotonic."))
    end

    if issorted(axis2)
        return data1, axis1, axis2
    elseif issorted(axis2; rev=true)
        axis2 = reverse(axis2)
        return reverse(data1, dims=2), axis1, axis2
    else
        throw(ArgumentError("Axis 2 must be monotonic."))
    end
end

@inline function _axes_match(a::Spectrum, b::Spectrum)
    return ((a.axis1 ≈ b.axis1) && (a.axis2 ≈ b.axis2))
end

@inline function _axes_match(a::OmnidirectionalSpectrum, b::OmnidirectionalSpectrum)
    return (a.axis ≈ b.axis)
end

@inline function check_axes_in_broadcast(args, sp::Union{Spectrum, OmnidirectionalSpectrum})
    type = (typeof(sp) <: Spectrum) ? Spectrum : OmnidirectionalSpectrum
    typename = (typeof(sp) <: Spectrum) ? "Spectrum axes" : "OmnidirectionalSpectrum axis"
    for arg in args
        if arg isa type
            arg === sp && continue
            _axes_match(sp, arg) || throw(DimensionMismatch(
                "$(typename) must match for broadcasting."
            ))
        elseif arg isa Broadcast.Broadcasted
            check_axes_in_broadcast(arg.args, sp)
        end
    end
    return nothing
end
