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
