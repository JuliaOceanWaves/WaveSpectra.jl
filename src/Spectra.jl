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
    frequency::Union{Nothing, AbstractVector{<:Number}} # TODO fix for Units
    spectrum::Union{Nothing, AbstractVector{<:Number}}
end

function spectrum(freqs::AbstractVector{<:Real}, vals::AbstractVector{<:Real})
    (size(freqs) == size(vals)) || error("Frequencies and Value Arrays are not same size!")
    extp_func = extrapolate(interpolate((freqs,), vals, Gridded(Linear())), 0)
    return spectrum(x->extp_func(x), freqs, vals) # TODO: add interpolation function
end
spectrum(func::Function) = spectrum(func, nothing, nothing)
spectrum(func::Function, freqs::AbstractVector{<:Real}) = spectrum(func, freqs, func(freqs))

# https://docs.juliaplots.org/latest/recipes/#series-recipes
# @recipe f(S::spectrum) = (S.frequencies, S.values)
# @recipe f(S::spectrum, freq::AbstractVector) = (freq, S(freq))

@recipe function f(S::spectrum)
    # marker := :circle
    (S.frequency, S.spectrum)
    
    # TODO: 2 options if vectors exist or not
end
@recipe function f(S::spectrum, freq::AbstractVector{<:Real})
    marker := :circle
    (freq, S(freq))
end

(S::spectrum)(freq::Real) =  S.func(freq)
# function (S::spectrum)(freq::Real)
#     if freq in S.frequencies return S.values[findfirst(S.frequencies .== freq)] end

#     b = findfirst(S.frequencies .> freq)
#     if (b == 1 || isnothing(b)) return 0 end

#     itp = interpolate((S.frequencies[b-1: b],), S.values[b-1:b], Gridded(Linear()))
#     return itp(freq)

# end

# TODO: move to struct def
function (S::spectrum)(freqs::AbstractVecOrMat{<:Real})
    #extp = extrapolate(interpolate((S.frequencies,), S.values, Gridded(Linear())), 0)
    return map(S.func, freqs)
end

# TODO: Remove and check dot notation
# function (S::spectrum)(freqs...)
#     !(typeof(freqs) <: Tuple{Vararg{Real}}) && error("Must pass arguments of type <: Real!")
#     extp = extrapolate(interpolate((S.frequencies,), S.values, Gridded(Linear())), 0)
#     return map(extp, freqs)
# end

include("statistics.jl")

end
