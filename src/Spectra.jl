module Spectra

export spectrum
export spectral_moment
export energy_period
export significant_waveheight

using Unitful: upreferred, uconvert, m, Hz, unit, Length, Time, Frequency, ustrip
using DimensionfulAngles: rad
using PhysicalConstants.CODATA2018: g_n as g
using Integrals: SampledIntegralProblem, TrapezoidalRule, IntegralProblem, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode
using Distributions: Chisq, quantile
using Interpolations: interpolate, extrapolate, scale, Gridded, Linear
using Plots

# Write your package code here.
struct spectrum
    func::Union{Nothing, Function}
    frequency::Union{Nothing, AbstractVector{<:Number}}
    spectrum::Union{Nothing, AbstractVector{<:Number}}
end

function spectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number})
    (size(freqs) == size(vals)) || error("Frequencies and Value Arrays are not same size!")
    extp_func = extrapolate(interpolate((freqs,), vals, Gridded(Linear())), 0)
    return spectrum(x->extp_func(x), freqs, vals) # TODO: add interpolation function
end
spectrum(func::Function) = spectrum(func, nothing, nothing)
spectrum(func::Function, freqs::AbstractVector{<:Number}) = spectrum(func, freqs, func(freqs))

# https://docs.juliaplots.org/latest/recipes/#series-recipes
# @recipe f(S::spectrum) = (S.frequencies, S.values)
# @recipe f(S::spectrum, freq::AbstractVector) = (freq, S(freq))

@recipe function f(S::spectrum)
    if !isnothing(S.frequency) && !isnothing(S.spectrum)
        marker := :circle
        (S.frequency, S.spectrum)
    else
        (S.func, 0, 10)
    end
end
@recipe function f(S::spectrum, args...)
    if !isnothing(S.frequency) && !isnothing(S.spectrum)
        marker := :circle
        (S.frequency, S.spectrum, args...)
    else
        (S.func, args...)
    end
end

@recipe function f(S::spectrum, freq::AbstractVector{<:Real})
    marker := :circle
    (freq, S(freq))
end
@recipe function f(S::spectrum, freq::AbstractVector{<:Real}, args...)
    marker := :circle
    (freq, S(freq), args...)
end

(S::spectrum)(freq::Real) =  S.func(freq)
(S::spectrum)(freqs::AbstractVecOrMat{<:Real}) = map(S.func, freqs)

include("statistics.jl")

end
