module WaveSpectra

export WaveSpectrum, isdiscrete

using Interpolations: interpolate, extrapolate, Gridded, Linear, InterpolationType #, scale
using Plots: @recipe
using Unitful: Quantity, Frequency, unit, Hz

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
isdiscrete(S::WaveSpectrum) = !isnothing(S.density) && !isnothing(S.frequency)
isunitful(S::WaveSpectrum) = hasmethod(S.func, (Frequency,))

(S::WaveSpectrum)(freq::Number) =  S.func(freq)
(S::WaveSpectrum)(freqs::AbstractVecOrMat{<:Number}) = map(S.func, freqs)

# TODO: Add costum error for discrete spectra interfaces if Spectrum has no vectors
_firststate() = 1
Base.getindex(S::WaveSpectrum, i::Int) = S.density[i]
Base.firstindex(S::WaveSpectrum) = _firststate()
Base.lastindex(S::WaveSpectrum) = length(S.density)
Base.length(S::WaveSpectrum) = length(S.density)
Base.size(S::WaveSpectrum) = size(S.density)
Base.iterate(S::WaveSpectrum) = (S.density[firstindex(S)], _firststate()+1)
function Base.iterate(S::WaveSpectrum, state::Int)
    return (state<=length(S)) ? (S.density[state], state+1) : nothing
end
Base.eltype(S::WaveSpectrum) = eltype(S.density)

# Plots recipes
function _labels(S::WaveSpectrum)
    xlabel = "Frequency"
    ylabel = "Spectral density"
    return (xlabel, ylabel)
end

@recipe function f(S::WaveSpectrum, args...)
    _xlabel, _ylabel = _labels(S)
    xlabel --> _xlabel
    ylabel --> _ylabel
    if isdiscrete(S)
        marker := :auto
        (S.frequency, S.density, args...)
    else
        if isempty(args)
            isunitful(S) ? (args=(0Hz, 10Hz)) : (args=(0, 10))
        end
        (S.func, args...)
    end
end

@recipe function f(S::WaveSpectrum, freq::AbstractVector{<:Number}, args...)
    _xlabel, _ylabel = _labels(S)
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (freq, S(freq), args...)
end

end
