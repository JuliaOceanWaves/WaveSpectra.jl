# integrate
function integrate(
    x::Spectrum,
    ax::Symbol,
    method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
)
    if (ax == :axis1) || (ax == x.axesnames[1])
        problem = SampledIntegralProblem(x.data, x.axis1; dim=1)
        axis = x.axis2
    elseif (ax == :axis2) || (ax == x.axesnames[2])
        problem = SampledIntegralProblem(x.data, x.axis2; dim=2)
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

function integrate(x::Spectrum, method::AbstractSampledIntegralAlgorithm=TrapezoidalRule())
    ax = (x.axestypes[1] == :direction) ? :axis1 : :axis2
    return integrate(integrate(x, ax, method), method)
end

function integrate(
    x::OmnidirectionalSpectrum,
    method::AbstractSampledIntegralAlgorithm=TrapezoidalRule()
)
    problem = SampledIntegralProblem(x.data, x.axis)
    sol = solve(problem, method)
    if Int(sol.retcode) ≠ 1
        throw(ProcessFailedException("Integration failed with return code: $(sol.retcode)"))
    end
    return sol.u
end

# Centered rectangular integration for evenly spaced Spectra
struct RectangularRule <: AbstractSampledIntegralAlgorithm end

struct RectangularUniformWeights{T} <: UniformWeights
    n::Int
    h::T
end

@inline Base.getindex(w::RectangularUniformWeights, i) = w.h

function find_weights(x::AbstractVector, ::RectangularRule)
    (x isa AbstractRange) && return RectangularUniformWeights(length(x), step(x))
    x_range = _convert_to_range(x)
    if (length(x)==length(x_range)) && isapprox(x, x_range)
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
    return range(start_val, stop=end_val, step=step_val)
end
