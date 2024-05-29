module Spectra

using Unitful: upreferred, uconvert, m, Hz, unit, Length, Time, Frequency, ustrip
using DimensionfulAngles: rad
using PhysicalConstants.CODATA2018: g_n as g
using Integrals: SampledIntegralProblem, TrapezoidalRule, IntegralProblem, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode
using Distributions: Chisq, quantile
using Interpolations: interpolate, extrapolate, scale, Gridded, Linear
using Plots
# Write your package code here.
struct spectra
    spectra_func::Union{Nothing, Function}
    frequencies::Union{Nothing, AbstractVector{<:Real}}
    values::Union{Nothing, AbstractVector{<:Real}}
end

spectra(freqs::AbstractVector{<:Real}, vals::AbstractVector{<:Real}) = spectra(nothing, freqs, vals)
spectra(func::Function) = spectra(func, nothing, nothing)
spectra(func::Function, freqs::AbstractVector{<:Real}) = spectra(func, freqs, func(freqs))

@recipe function f(S::spectra)
    marker := :circle
    (S.frequencies, S.values)
end
@recipe function f(S::spectra, freq::AbstractVector)
    marker := :circle
    (freq, S(freq))
end

function (S::spectra)(freq::Real)
    if freq in S.frequencies 
        return S.values[findfirst(S.frequencies .== freq)]
    else
        b = findfirst(S.frequencies .> freq)
        if (b == 1 || isnothing(b))
            return 0
        else
            itp = interpolate((S.frequencies[b-1: b],), S.values[b-1:b], Gridded(Linear()))
            return itp(freq)
        end
    end
end

function (S::spectra)(freqs...)
    !(typeof(freqs) <: Tuple{Vararg{Real}}) && error("Must pass arguments of type <: Real!")
    extp = extrapolate(interpolate((S.frequencies,), S.values, Gridded(Linear())), 0)
    return extp(collect(freqs))
end

function (S::spectra)(freqs::AbstractVector{<:Real})
    extp = extrapolate(interpolate((S.frequencies,), S.values, Gridded(Linear())), 0)
    return extp(collect(freqs))
end

# spectral statistics
function integrate(y::AbstractVector, x::AbstractVector, method::AbstractIntegralAlgorithm=TrapezoidalRule())
    sol = solve(SampledIntegralProblem(y, x), method)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractMatrix, x::AbstractVector, method::AbstractIntegralAlgorithm=TrapezoidalRule())
    sol = solve(SampledIntegralProblem(y, x; dim=2), method)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractVector, x::AbstractVector, a::Number, b::Number, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    f = extrapolate(interpolate((x,), y, Gridded(Linear())), 0.0)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = solve(IntegralProblem((x, p)->f(x), a, b), method; abstol)
    sol.retcode ≠ ReturnCode.Success && error("solution unsuccessful with code: $(sol.retcode)")
    return sol.u
end

function integrate(y::AbstractMatrix, x::AbstractVector, a::Number, b::Number, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = ones(eltype(x[1]*y[1,1]), size(y)[1])
    for (i, yᵢ) in enumerate(eachrow(y))
        fᵢ = extrapolate(interpolate((x,), yᵢ, Gridded(Linear())), 0.0)
        solᵢ = solve(IntegralProblem((x, p)->fᵢ(x), a, b), method; abstol)
        solᵢ.retcode ≠ ReturnCode.Success && error("solution $i unsuccessful with code: $(sol.retcode)")
        sol[i] = solᵢ.u
    end
    return sol
end

function integrate(y::AbstractMatrix, x::AbstractVector, a::Union{Number, AbstractVector}, b::Union{Number, AbstractVector}, method::AbstractIntegralAlgorithm=QuadGKJL(); abstol=nothing)
    isnothing(abstol) && (abstol=1.0e-8*unit(eltype(x))*unit(eltype(y)))
    sol = ones(eltype(x[1]*y[1,1]), size(y)[1])
    a = (length(a)≠1) ? a : ones(size(y)[1])*a
    b = (length(b)≠1) ? b : ones(size(y)[1])*b
    for (i, yᵢ) in enumerate(eachrow(y))
        fᵢ = extrapolate(interpolate((x,), yᵢ, Gridded(Linear())), 0.0)
        solᵢ = solve(IntegralProblem((x, p)->fᵢ(x), a[i], b[i]), method; abstol)
        solᵢ.retcode ≠ ReturnCode.Success && error("solution $i unsuccessful with code: $(sol.retcode)")
        sol[i] = solᵢ.u
    end
    return sol
end

∫ = integrate

function spectral_moment(S::AbstractVector, f::AbstractVector, n::Integer; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    return ∫((S.*(f.^n)), f, method)
end

function spectral_moment(S::AbstractMatrix, f::AbstractVector, n::Integer; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    return ∫((S.*(f.^n)'), f, method)
end

function spectral_moment(S::AbstractVector, f::AbstractVector, n::Integer, a::Number, b::Number; method::AbstractIntegralAlgorithm=QuadGKJL())
    return ∫((S.*(f.^n)), f, a, b, method)
end

function spectral_moment(S::AbstractMatrix, f::AbstractVector, n::Integer, a::Union{Number, AbstractVector}, b::Union{Number, AbstractVector}; method::AbstractIntegralAlgorithm=QuadGKJL())
    return ∫((S.*(f.^n)'), f, a, b, method)
end

function energy_period(S::AbstractVecOrMat, f::AbstractVector; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(S, f, -1; method)
    m_0 = spectral_moment(S, f, 0; method)
    return upreferred.(m_n1./m_0)
end

function energy_period(S::spectra; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(S.values, S.frequencies, -1; method)
    m_0 = spectral_moment(S.values, S.frequencies, 0; method)
    return upreferred.(m_n1./m_0)
end

function significant_waveheight(S::AbstractVecOrMat, f::AbstractVector; method::AbstractIntegralAlgorithm=TrapezoidalRule())
	m_0 = spectral_moment(S, f, 0; method)
    return 4(.√m_0)
end

function significant_waveheight(S::spectra; method::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_0 = spectral_moment(S.values, S.frequencies, 0; method)
    return 4(.√m_0)
end

export spectra 
export spectral_moment
export energy_period
export significant_waveheight
end
