# Support functions

# Ensure an array contains the same types and units
@inline function _check_typeconsistency(x::AbstractArray)::Nothing
    consistent = all(y -> typeof(y) == eltype(x), x)
    !consistent && throw(ArgumentError("All elements of array must be of same type."))
    return nothing
end

# Different checks on axes
@inline function _ensure_increasing_axes(data, axis1, axis2)
    if issorted(axis1)
        data1 = data
    elseif issorted(axis1; rev = true)
        axis1 = reverse(axis1)
        data1 = reverse(data, dims = 1)
    else
        throw(ArgumentError("Axis 1 must be monotonic."))
    end

    if issorted(axis2)
        return data1, axis1, axis2
    elseif issorted(axis2; rev = true)
        axis2 = reverse(axis2)
        return reverse(data1, dims = 2), axis1, axis2
    else
        throw(ArgumentError("Axis 2 must be monotonic."))
    end
end

@inline function _ensure_increasing_axis(data, axis)
    if issorted(axis)
        return data, axis
    elseif issorted(axis; rev = true)
        axis = reverse(axis)
        return reverse(data), axis
    else
        throw(ArgumentError("Axis must be monotonic."))
    end
end

@inline function _check_strictly_positive_finite_spectral_axis(axis::AbstractVector)
    if first(axis) <= zero(first(axis))
        throw(ArgumentError("Spectral-variable axis values must be positive."))
    elseif !all(isfinite, axis)
        throw(ArgumentError("Spectral-variable axis values must be finite."))
    end
    return nothing
end

@inline function _check_spectrum_axes(axis1::AbstractVector, axis2::AbstractVector)
    axis1_is_spectral = isspectralvariable(axis1)
    axis2_is_spectral = isspectralvariable(axis2)
    if axis1_is_spectral && axis2_is_spectral
        return :cartesian, axis1, axis2
    elseif axis1_is_spectral && isdirection(axis2)
        _check_strictly_positive_finite_spectral_axis(axis1)
        return :polar, axis1, axis2
    elseif axis2_is_spectral && isdirection(axis1)
        @warn "Swapping order of axes to have direction as second axis."
        _check_strictly_positive_finite_spectral_axis(axis2)
        return :polar, axis2, axis1
    end
    throw(ArgumentError("Axes must define a cartesian or polar coordinate."))
end

# Add the zero-spectral-value (`S(0) = 0`)
@inline function _prepend_zero_axis(axis::AbstractVector)
    out = Vector{eltype(axis)}(undef, length(axis) + 1)
    out[1] = zero(first(axis))
    copyto!(out, 2, axis, 1, length(axis))
    return out
end

@inline function _prepend_zero_data(data::AbstractVector)
    out = Vector{eltype(data)}(undef, length(data) + 1)
    out[1] = zero(eltype(data))
    copyto!(out, 2, data, 1, length(data))
    return out
end

@inline function _prepend_zero_data(data::AbstractMatrix; dim::Int = 1)
    if dim == 1
        out = Matrix{eltype(data)}(undef, size(data, 1) + 1, size(data, 2))
        out[1, :] .= zero(eltype(data))
        out[2:end, :] = data
        return out
    elseif dim == 2
        out = Matrix{eltype(data)}(undef, size(data, 1), size(data, 2) + 1)
        out[:, 1] .= zero(eltype(data))
        out[:, 2:end] = data
        return out
    end
    throw(ArgumentError("`dim` must be 1 or 2."))
end

# Create the appropriate indexes for AxisArrays when indexing a spectrum using kwargs.
@inline function _update_kwargs(d::AbstractDict)
    kwargs = Dict()
    for (k, v) in d
        kwargs[k] = _axis_selector(v)
    end
    return kwargs
end

@inline _axis_selector(v::Integer) = (v:v)
@inline _axis_selector(v::Quantity) = (v .. v)
@inline _axis_selector(v) = v

@inline _update_indices(I...) = map(_axis_selector, I)

# Check whether two spectra share the same axes
@inline function _axes_match(a::AbstractSpectrum, b::AbstractSpectrum)
    return ((a.axis1 ≈ b.axis1) && (a.axis2 ≈ b.axis2))
end

@inline function _axes_match(
        a::AbstractOmnidirectionalSpectrum,
        b::AbstractOmnidirectionalSpectrum
)
    return (a.axis ≈ b.axis)
end

# Create a new spectrum of the same type but with different data.
@inline _typewrapper(::Type{T}) where {T} = Base.typename(T).wrapper

@inline function _rebuild_spectrum(
        x::AbstractSpectrum,
        data::AbstractMatrix,
        axis1::AbstractVector,
        axis2::AbstractVector
)
    return _typewrapper(typeof(x))(data, axis1, axis2)
end

@inline function _rebuild_spectrum(
        x::AbstractOmnidirectionalSpectrum,
        data::AbstractVector,
        axis::AbstractVector
)
    return _typewrapper(typeof(x))(data, axis)
end

# Make broadcasting work.
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

@inline function _check_axes_in_broadcast(
        args,
        sp::Union{AbstractSpectrum, AbstractOmnidirectionalSpectrum}
)
    type = (typeof(sp) <: AbstractSpectrum) ?
           AbstractSpectrum :
           AbstractOmnidirectionalSpectrum
    typename = (typeof(sp) <: Spectrum) ? "Spectrum axes" : "OmnidirectionalSpectrum axis"
    for arg in args
        if arg isa type
            arg === sp && continue
            _axes_match(sp, arg) || throw(DimensionMismatch(
                "$(typename) must match for broadcasting."
            ))
        elseif arg isa Broadcast.Broadcasted
            _check_axes_in_broadcast(arg.args, sp)
        end
    end
    return nothing
end

# Convert vector to range by assuming it is evenly spaced. Will not fail if it isn't.
function _convert_to_range(x::AbstractVector)
    (length(x) == 0) && return nothing
    (length(x) == 1) && return (x:x)
    start_val = x[begin]
    end_val = x[end]
    step_val = x[2] - x[1]
    return range(start_val, stop = end_val, step = step_val)
end
