module WaveSpectra

export WaveSpectrum, isdiscrete

using Interpolations: interpolate, extrapolate, Gridded, Linear, InterpolationType #, scale
using Plots: @recipe

# Struct definition
struct WaveSpectrum
    func::Union{Nothing, Function}
    frequency::Union{Nothing, AbstractVector{<:Number}}
    density::Union{Nothing, AbstractVector{<:Number}}

    function WaveSpectrum(func, frequency, density)
        if !(isnothing(frequency) && isnothing(density)) && (size(frequency) !== size(density))
            error("Frequencies and Value Arrays are not same size!")
        end
        return new(func, frequency, density)
    end
end

# Constructors
function WaveSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; interpmode::InterpolationType=Gridded(Linear()))
    extrap_func = extrapolate(interpolate((freqs,), vals, interpmode), 0)
    return WaveSpectrum(x->extrap_func(x), freqs, vals)
end

WaveSpectrum(func::Function) = WaveSpectrum(func, nothing, nothing)
WaveSpectrum(func::Function, freqs::AbstractVector{<:Number}) = WaveSpectrum(func, freqs, func(freqs))

# Interface
# TODO: Add costum error for discrete spectra interfaces if Spectrum has no vectors
(S::WaveSpectrum)(freq::Real) =  S.func(freq)
(S::WaveSpectrum)(freqs::AbstractVecOrMat{<:Real}) = map(S.func, freqs)

Base.getindex(S::WaveSpectrum, i::Int) = S.density[i]
Base.firstindex(S::WaveSpectrum) = 1
Base.lastindex(S::WaveSpectrum) = length(S.density)
Base.length(S::WaveSpectrum) = length(S.density)
Base.size(S::WaveSpectrum) = size(S.density)

Base.iterate(S::WaveSpectrum) = (S.density[firstindex(S)], 2)
function Base.iterate(S::WaveSpectrum, state::Int)
    return (state<=length(S)) ? (S.density[state], state+1) : nothing
end

isdiscrete(S::WaveSpectrum) = !isnothing(S.density) && !isnothing(S.frequency)

# Plots recipes
# TODO: add xlabel and ylabel
@recipe function f(S::WaveSpectrum, args...)
    if isdiscrete(S)
        marker := :auto
        (S.frequency, S.density, args...)
    else
        isempty(args) && (args=(0, 10))
        (S.func, args...)
    end
end

@recipe function f(S::WaveSpectrum, freq::AbstractVector{<:Real}, args...)
    marker := :auto
    (freq, S(freq), args...)
end

end
