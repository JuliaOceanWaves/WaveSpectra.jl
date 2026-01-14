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


# split into omnidirectional & spread function
function omnidirectional_spectrum(x::Spectrum)
    ispolar(x) && return integrate(x, :direction)
    throw(ArgumentError("Spectrum must be in polar coordinates."))
end

function spread_function(x::Spectrum)
    return Spectrum(x.data ./ omnidirectional_spectrum(x), AxisArrays.axes(x)...)
end

function split_spectrum(x::Spectrum)
    omnidirectional = omnidirectional_spectrum(x)
    spread = Spectrum(x ./ omnidirectional, AxisArrays.axes(x)...)
    return (omnidirectional, spread)
end
