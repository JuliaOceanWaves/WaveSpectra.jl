# integrate spectra

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

"""
    integrate(
        x::AbstractSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )
    integrate(
        x::AbstractSpectrum,
        ax::Symbol;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )
    integrate(
        x::AbstractSpectrum,
        ax::Symbol,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
    )
    integrate(
        x::AbstractSpectrum,
        axrange1::Tuple,
        axrange2::Tuple;
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
    )

Integrate a 2D [`Spectrum`] along one axis.

The axis of integration `ax` can be `:axis1`, `:axis2`, or the corresponding entry in
`x.axesnames`.
If the remaining axis describes a spatial or temporal frequency, the result is returned as
an `OmnidirectionalSpectrum`, else it is returned as a vector and a warning is raised.
If no axis is specified a double integration, over both axis is performed returning a
scalar quantity.

When provided, `axrange`, `axrange1`, and `axrange2` are `(start, stop)` tuples selecting
the closed interval of the corresponding axis.
The bounds can be integer indices or `Quantity`.
If the bounds are quantities, the next value equal to or larger than `start`, and the next
value equal to or smaller than `end` are used as the bounds of integration.
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
        problem = SampledIntegralProblem(x.data, x.axis1; dim = 1)
        axis = x.axis2
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        problem = SampledIntegralProblem(x.data, x.axis2; dim = 2)
        axis = x.axis1
    else
        throw(ArgumentError("Unknown axis."))
    end
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    if axestypes(axis) == :direction
        @warn ("Integration: The resulting vector is a function of direction, and is not " *
               "returned as an 'OmnidirectionalSpectrum'."
        )
        return sol.u
    end
    return OmnidirectionalSpectrum(sol.u, axis)
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
        problem = SampledIntegralProblem(
            @view(x.data[idx, :]), @view(x.axis1[idx]); dim = 1)
        axis = x.axis2
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        idx = _axis_index_range(x.axis2, axrange)
        _print_integration_bounds(x.axis2, idx)
        problem = SampledIntegralProblem(
            @view(x.data[:, idx]), @view(x.axis2[idx]); dim = 2)
        axis = x.axis1
    else
        throw(ArgumentError("Unknown axis."))
    end
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    if axestypes(axis) == :direction
        @warn ("Integration: The resulting vector is a function of direction, and is not " *
               "returned as an 'OmnidirectionalSpectrum'."
        )
        return sol.u
    end
    return OmnidirectionalSpectrum(sol.u, axis)
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
    integrate(
        x::AbstractOmnidirectionalSpectrum;
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
    )

    integrate(
        x::AbstractOmnidirectionalSpectrum,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
    )

Integrate an `OmnidirectionalSpectrum` and return a scalar quantity.

When provided, `axrange` is a `(start, stop)` tuple selecting the closed interval of the
axis.

Integer bounds are treated as indices; `Quantity` bounds are matched against the
axis values.

When provided, `axrange`is a `(start, stop)` tuplee selecting the closed interval of the
axis.
The bounds can be integer indices or `Quantity`.
If the bounds are quantities, the next value equal to or larger than `start`, and the next
value equal to or smaller than `end` are used as the bounds of integration.
"""
function integrate(x::AbstractOmnidirectionalSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    problem = SampledIntegralProblem(x.data, x.axis)
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    return sol.u
end

function integrate(
        x::AbstractOmnidirectionalSpectrum,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    idx = _axis_index_range(x.axis, axrange)
    _print_integration_bounds(x.axis, idx)
    problem = SampledIntegralProblem(@view(x.data[idx]), @view(x.axis[idx]))
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    return sol.u
end

# centered rectangular integration for evenly spaced Spectra
"""
    RectangularRule

Rectangular integration method for evenly spaced samples.

Every sample has a weight equal to the sample spacing.
Meant to work with `SampledIntegralProblem` from the `Integrals.jl` package.
"""
struct RectangularRule <: AbstractSampledIntegralAlgorithm end

struct RectangularUniformWeights{T} <: UniformWeights
    n::Int
    h::T
end

@inline Base.getindex(w::RectangularUniformWeights, i) = w.h

function find_weights(x::AbstractVector, ::RectangularRule)
    (x isa AbstractRange) && return RectangularUniformWeights(length(x), step(x))
    x_range = _convert_to_range(x)
    if (length(x) == length(x_range)) && isapprox(x, x_range)
        return RectangularUniformWeights(length(x_range), step(x_range))
    else
        throw(ArgumentError("Integrand must be evenly spaced."))
    end
    return TrapezoidalNonuniformWeights(x)
end

function _convert_to_range(x::AbstractVector)
    (length(x) == 0) && return nothing
    (length(x) == 1) && return (x:x)
    start_val = x[begin]
    end_val = x[end]
    step_val = x[2] - x[1]
    return range(start_val, stop = end_val, step = step_val)
end
