# integrate spectra

"""
    integrate(
        x::AbstractSpectrum;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )
    integrate(
        x::AbstractSpectrum,
        ax::Symbol;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )
    integrate(
        x::AbstractSpectrum,
        ax::Symbol,
        axrange::Tuple;
        include_zero::Bool = false,
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
For polar spectra, integrating over direction returns an `OmnidirectionalSpectrum`, while
integrating over the spectral-variable axis returns an `AxisArray` indexed by direction.
For cartesian spectra, integrating over a single axis returns an `AxisArray` describing a
marginal spectrum.
If no axis is specified a double integration, over both axis is performed returning a
scalar quantity.

When provided, `axrange`, `axrange1`, and `axrange2` are `(start, stop)` tuples selecting
the closed interval of the corresponding axis.
The bounds can be integer indices or `Quantity`.
If the bounds are quantities, the next value equal to or larger than `start`, and the next
value equal to or smaller than `end` are used as the bounds of integration.

The option `include_zero` will add the point `(S(0,0) = 0)` before integrating, if the
spectrum is in polar coordinates.
"""
function integrate(x::AbstractSpectrum;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    ax = ispolar(x) ? :axis2 : :axis1
    y = _integrate_spectrum_axis(x, ax; method, warn_non_spectrum = false, include_zero)
    ispolar(x) && return integrate(y; method, include_zero)
    return integrate(y; method)
end

function integrate(x::AbstractSpectrum, ax::Symbol;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    return _integrate_spectrum_axis(x, ax; method, warn_non_spectrum = true, include_zero)
end

function integrate(
        x::AbstractSpectrum,
        ax::Symbol,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    return _integrate_spectrum_axis(x, ax, axrange; method, warn_non_spectrum = true)
end

function integrate(
        x::AbstractSpectrum,
        axrange1::Tuple,
        axrange2::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    ax = (x.axestypes[1] == :direction) ? :axis1 : :axis2
    if ax == :axis1
        y = _integrate_spectrum_axis(x, :axis1, axrange1; method, warn_non_spectrum = false)
        return integrate(y, axrange2; method)
    end
    y = _integrate_spectrum_axis(x, :axis2, axrange2; method, warn_non_spectrum = false)
    return integrate(y, axrange1; method)
end

"""
    integrate(
        x::AxisArray;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )
    integrate(
        x::AxisArray,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Integrate a 1D `AxisArray` and return a scalar quantity.

When provided, `axrange` is a `(start, stop)` tuple selecting the closed interval of the
axis.
The bounds can be integer indices or `Quantity`.
If the bounds are quantities, the next value equal to or larger than `start`, and the next
value equal to or smaller than `stop` are used as the bounds of integration.
"""
function integrate(x::AxisArray;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    axis = _only_axis(x)
    return _solve_sampled_integral(x.data, axis; method)
end

function integrate(x::AxisArray, axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    axis = _only_axis(x)
    idx = _axis_index_range(axis, axrange...)
    @info string("Integrating over $(axis[first(idx)]) to $(axis[last(idx)]).")
    return _solve_sampled_integral(@view(x.data[idx]), @view(axis[idx]); method)
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
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    data = include_zero ? _prepend_zero_data(x.data) : x.data
    axis = include_zero ? _prepend_zero_axis(x.axis) : x.axis
    return _solve_sampled_integral(data, axis; method)
end

function integrate(
        x::AbstractOmnidirectionalSpectrum,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
)
    idx = _axis_index_range(x.axis, axrange...)
    @info string("Integrating over $(x.axis[first(idx)]) to $(x.axis[last(idx)]).")
    return _solve_sampled_integral(@view(x.data[idx]), @view(x.axis[idx]); method)
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

# utilities
function _solve_sampled_integral(data, axis;
        method::AbstractSampledIntegralAlgorithm,
        dim::Union{Nothing, Int} = nothing
)
    problem = isnothing(dim) ? SampledIntegralProblem(data, axis) :
              SampledIntegralProblem(data, axis; dim)
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    return sol.u
end

function _only_axis(x::AxisArray)
    axis = axisvalues(x)
    if length(axis) ≠ 1
        throw(DimensionMismatch(
            "Exactly one axis is required for integration of an 'AxisArray'.")
        )
    end
    return axis[1]
end

@inline function _axis_index_range(axis::AbstractVector, start::Int, stop::Int)
    (start > stop) && throw(ArgumentError("Invalid index range."))
    return start:stop
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

function _spectrum_axis_info(x::AbstractSpectrum, ax::Symbol)
    if (ax == :axis1) || (ax == x.axesnames[1])
        return (
            dim = 1,
            integrated_axis = x.axis1,
            remaining_axis = x.axis2,
            remaining_axis_name = x.axesnames[2]
        )
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        return (
            dim = 2,
            integrated_axis = x.axis2,
            remaining_axis = x.axis1,
            remaining_axis_name = x.axesnames[1]
        )
    end
    throw(ArgumentError("Unknown axis."))
end

function _slice_spectrum_data(data, idx, dim::Int)
    (dim == 1) && return @view(data[idx, :])
    (dim == 2) && return @view(data[:, idx])
    throw(ArgumentError("Unknown axis dimension."))
end

function _integrate_spectrum_axis(
        x::AbstractSpectrum,
        ax::Symbol;
        method::AbstractSampledIntegralAlgorithm,
        warn_non_spectrum::Bool,
        include_zero::Bool
)
    axis_info = _spectrum_axis_info(x, ax)
    return _solve_spectrum_axis_integration(
        x, x.data, axis_info; method, warn_non_spectrum, include_zero
    )
end

function _integrate_spectrum_axis(
        x::AbstractSpectrum,
        ax::Symbol,
        axrange::Tuple;
        method::AbstractSampledIntegralAlgorithm,
        warn_non_spectrum::Bool
)
    axis_info = _spectrum_axis_info(x, ax)
    integrated_axis = axis_info.integrated_axis
    idx = _axis_index_range(integrated_axis, axrange...)
    @info string(
        "Integrating over $(integrated_axis[first(idx)]) to $(integrated_axis[last(idx)])."
    )
    return _solve_spectrum_axis_integration(
        x,
        _slice_spectrum_data(x.data, idx, axis_info.dim),
        merge(axis_info, (; integrated_axis = @view(integrated_axis[idx])));
        method,
        warn_non_spectrum
    )
end

function _solve_spectrum_axis_integration(
        x::AbstractSpectrum,
        data,
        axis_info;
        method::AbstractSampledIntegralAlgorithm,
        warn_non_spectrum::Bool,
        include_zero::Bool = false
)
    integrated_axis = axis_info.integrated_axis
    if include_zero && ispolar(x) && (axis_info.dim == 1)
        data = _prepend_zero_data(data; dim = axis_info.dim)
        integrated_axis = _prepend_zero_axis(integrated_axis)
    end
    result = _solve_sampled_integral(data, integrated_axis; method, dim = axis_info.dim)
    remaining_axis = axis_info.remaining_axis
    if ispolar(x) && isdirection(remaining_axis)
        warn_non_spectrum &&
            @warn("Integration: The result is a function of direction, and is returned as " *
                  "an 'AxisArray' rather than an 'OmnidirectionalSpectrum'.")
        return AxisArray(result, Axis{axis_info.remaining_axis_name}(remaining_axis))
    elseif ispolar(x)
        return OmnidirectionalSpectrum(result, remaining_axis)
    end
    warn_non_spectrum &&
        @warn("Integration: The result is a marginal spectrum, and is returned as an " *
              "'AxisArray' rather than an 'OmnidirectionalSpectrum'.")
    return AxisArray(result, Axis{axis_info.remaining_axis_name}(remaining_axis))
end
