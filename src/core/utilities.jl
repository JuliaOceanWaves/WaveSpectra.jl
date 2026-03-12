# Support functions

@inline function _check_typeconsistency(x::AbstractArray)::Nothing
    consistent = all(y -> typeof(y) == eltype(x), x)
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

@inline function _axes_match(a::AbstractSpectrum, b::AbstractSpectrum)
    return ((a.axis1 ≈ b.axis1) && (a.axis2 ≈ b.axis2))
end

@inline function _axes_match(
        a::AbstractOmnidirectionalSpectrum,
        b::AbstractOmnidirectionalSpectrum
)
    return (a.axis ≈ b.axis)
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
