# integrate spectra

abstract type AbstractSampledIntegralAlgorithm end
struct TrapezoidalRule <: AbstractSampledIntegralAlgorithm end

# centered rectangular integration for evenly spaced spectra
"""
    RectangularRule

Rectangular integration method for evenly spaced samples.

Every sample has a weight equal to the sample spacing.
"""
struct RectangularRule <: AbstractSampledIntegralAlgorithm end

@inline _axis_index_range(axis::AbstractVector, start::Int, stop::Int) = start:stop
@inline function _axis_index_range(axis::AbstractVector, axrange::Tuple)
    return _axis_index_range(axis, axrange...)
end

@inline function _axis_index_range(
        axis::AbstractVector{<:Quantity},
        start::Quantity,
        stop::Quantity
)
    i = searchsortedfirst(axis, start)
    j = searchsortedlast(axis, stop)
    (i > j) && throw(ArgumentError("Integration interval does not overlap the axis."))
    return i:j
end

@inline _print_integration_bounds(axis::AbstractVector, idx) = nothing
@inline function _print_integration_bounds(axis::AbstractVector{<:Quantity}, idx)
    println("Integrating over ", axis[first(idx)], " to ", axis[last(idx)], ".")
    return nothing
end

function _convert_to_range(x::AbstractVector)
    (length(x) == 0) && return nothing
    (length(x) == 1) && return (x:x)
    start_val = x[begin]
    end_val = x[end]
    step_val = x[2] - x[1]
    return range(start_val, stop = end_val, step = step_val)
end

function _trapezoidal_weights(axis::AbstractVector)
    n = length(axis)
    (n < 2) && throw(ArgumentError("Integration axis must contain at least two samples."))
    dx = diff(axis)
    w = Vector{eltype(dx)}(undef, n)
    w[begin] = dx[begin] / 2
    for i in 2:(n - 1)
        w[i] = (dx[i - 1] + dx[i]) / 2
    end
    w[end] = dx[end] / 2
    return w
end

function _rectangular_weights(axis::AbstractVector)
    n = length(axis)
    (n < 2) && throw(ArgumentError("Integration axis must contain at least two samples."))
    if axis isa AbstractRange
        h = step(axis)
        return fill(h, n)
    end
    axis_range = _convert_to_range(axis)
    if (length(axis) == length(axis_range)) && isapprox(axis, axis_range)
        h = step(axis_range)
        return fill(h, n)
    end
    throw(ArgumentError("Integrand must be evenly spaced."))
end

_integration_weights(axis::AbstractVector, ::TrapezoidalRule) = _trapezoidal_weights(axis)
_integration_weights(axis::AbstractVector, ::RectangularRule) = _rectangular_weights(axis)

function _sampled_integral(data::AbstractVector, axis::AbstractVector, method)
    w = _integration_weights(axis, method)
    return sum(data .* w)
end

function _sampled_integral(data::AbstractMatrix, axis::AbstractVector, dim::Int, method)
    w = _integration_weights(axis, method)
    if dim == 1
        length(w) == size(data, 1) || throw(DimensionMismatch("axis length must match data rows"))
        return [sum(@view(data[:, j]) .* w) for j in axes(data, 2)]
    elseif dim == 2
        length(w) == size(data, 2) || throw(DimensionMismatch("axis length must match data columns"))
        return [sum(@view(data[i, :]) .* w) for i in axes(data, 1)]
    end
    throw(ArgumentError("Integration dimension must be 1 or 2."))
end

"""
    integrate(x::AbstractSpectrum; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    integrate(x::AbstractSpectrum, ax::Symbol; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    integrate(x::AbstractSpectrum, ax::Symbol, axrange::Tuple; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    integrate(x::AbstractSpectrum, axrange1::Tuple, axrange2::Tuple; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())

Integrate a 2D [`Spectrum`] along one axis.
"""
function integrate(x::AbstractSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    ax = (x.axestypes[1] == :direction) ? :axis1 : :axis2
    return integrate(integrate(x, ax; method); method)
end

function integrate(x::AbstractSpectrum, ax::Symbol;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    if (ax == :axis1) || (ax == x.axesnames[1])
        values = _sampled_integral(x.data, x.axis1, 1, method)
        axis = x.axis2
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        values = _sampled_integral(x.data, x.axis2, 2, method)
        axis = x.axis1
    else
        throw(ArgumentError("Unknown axis."))
    end
    if axestypes(axis) == :direction
        @warn ("Integration: The resulting vector is a function of direction, and is not " *
               "returned as an 'OmnidirectionalSpectrum'.")
        return values
    end
    return OmnidirectionalSpectrum(values, axis)
end

function integrate(
        x::AbstractSpectrum,
        ax::Symbol,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    if (ax == :axis1) || (ax == x.axesnames[1])
        idx = _axis_index_range(x.axis1, axrange)
        _print_integration_bounds(x.axis1, idx)
        values = _sampled_integral(@view(x.data[idx, :]), @view(x.axis1[idx]), 1, method)
        axis = x.axis2
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        idx = _axis_index_range(x.axis2, axrange)
        _print_integration_bounds(x.axis2, idx)
        values = _sampled_integral(@view(x.data[:, idx]), @view(x.axis2[idx]), 2, method)
        axis = x.axis1
    else
        throw(ArgumentError("Unknown axis."))
    end
    if axestypes(axis) == :direction
        @warn ("Integration: The resulting vector is a function of direction, and is not " *
               "returned as an 'OmnidirectionalSpectrum'.")
        return values
    end
    return OmnidirectionalSpectrum(values, axis)
end

function integrate(
        x::AbstractSpectrum,
        axrange1::Tuple,
        axrange2::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    ax = (x.axestypes[1] == :direction) ? :axis1 : :axis2
    if ax == :axis1
        return integrate(integrate(x, :axis1, axrange1; method), axrange2; method)
    end
    return integrate(integrate(x, :axis2, axrange2; method), axrange1; method)
end

"""
    integrate(x::AbstractOmnidirectionalSpectrum; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    integrate(x::AbstractOmnidirectionalSpectrum, axrange::Tuple; method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())

Integrate an `OmnidirectionalSpectrum` and return a scalar quantity.
"""
function integrate(x::AbstractOmnidirectionalSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    return _sampled_integral(x.data, x.axis, method)
end

function integrate(
        x::AbstractOmnidirectionalSpectrum,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    idx = _axis_index_range(x.axis, axrange)
    _print_integration_bounds(x.axis, idx)
    return _sampled_integral(@view(x.data[idx]), @view(x.axis[idx]), method)
end
