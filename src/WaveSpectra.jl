module WaveSpectra

export WaveSpectrum, isdiscrete

using Interpolations: interpolate, extrapolate, AbstractExtrapolation, linear_interpolation, cubic_spline_interpolation
using Plots: @recipe
using Unitful: Quantity, Frequency, unit, Hz
import Base: Exception, showerror

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
# function WaveSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; interpmode::InterpolationType=Gridded(Linear()))
#     if (typeof(interpmode) <: Gridded) 
#         extrap_func = extrapolate(interpolate((freqs,), vals, interpmode), 0)
#     elseif (typeof(interpmode) <: BSpline)
#         (typeof(freqs) <: AbstractRange) || error("If using Bspline, frequencies must be <:AbstractRange")
#         extrap_func = extrapolate(scale(interpolate(vals, interpmode), freqs), 0)
#     else
#         error("Unrecognized interpolation mode! Please use Gridded or Bspline.")
#     end
#     return WaveSpectrum(x->extrap_func(x), freqs, vals)
# end
function WaveSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; interpmode::Function=linear_interpolation)
    if !(hasmethod(interpmode, (AbstractVector, AbstractVector)) || hasmethod(interpmode, (AbstractRange, AbstractVector)))
        error("Must be a proper interpolation function whose positional arguments are (AbstractVector, AbstractVector) or (AbstractRange, AbstractVector)")
    end
    if (interpmode == cubic_spline_interpolation) 
        (typeof(freqs) <: AbstractRange) || error("If using cubic spline, frequencies must be <:AbstractRange")  
    end
    extrap_func = interpmode(freqs, vals; extrapolation_bc = 0)
    return WaveSpectrum(x->extrap_func(x), freqs, vals)
end

WaveSpectrum(func::Function) = WaveSpectrum(func, nothing, nothing)
WaveSpectrum(func::Function, freqs::AbstractVector{<:Number}) = WaveSpectrum(func, freqs, func(freqs))

# Interface
isdiscrete(S::WaveSpectrum) = !isnothing(S.density) && !isnothing(S.frequency)
isunitful(S::WaveSpectrum) = hasmethod(S.func, (Frequency,)) && !(hasmethod(S.func, (Any,)))

(S::WaveSpectrum)(freq::Number) =  S.func(freq)
(S::WaveSpectrum)(freqs::AbstractVecOrMat{<:Number}) = map(S.func, freqs)

mutable struct DiscreteError <: Exception
    message::String
    DiscreteError() = new("WaveSpectrum instance was treated as discrete when it is continuous!")
end

_firststate() = 1
Base.showerror(io::IO, e::DiscreteError) = print(io, "DiscreteError: $(e.message)")
Base.getindex(S::WaveSpectrum, i::Int) = isdiscrete(S) ? S.density[i] : throw(DiscreteError()) 
Base.firstindex(S::WaveSpectrum) = isdiscrete(S) ? _firststate() : throw(DiscreteError()) 
Base.lastindex(S::WaveSpectrum) = isdiscrete(S) ? length(S.density) : throw(DiscreteError()) 
Base.length(S::WaveSpectrum) = isdiscrete(S) ? length(S.density) : throw(DiscreteError()) 
Base.size(S::WaveSpectrum) = isdiscrete(S) ? size(S.density) : throw(DiscreteError()) 
Base.iterate(S::WaveSpectrum) = isdiscrete(S) ? (S.density[firstindex(S)], _firststate()+1) : throw(DiscreteError()) 
Base.eltype(S::WaveSpectrum) = isdiscrete(S) ? eltype(S.density) : throw(DiscreteError())
function Base.iterate(S::WaveSpectrum, state::Int)
    isdiscrete(S) || throw(DiscreteError())
    return (state<=length(S)) ? (S.density[state], state+1) : nothing
end 

# Plots recipes
_labels() = ("Frequency", "Spectral density")
@recipe function f(S::WaveSpectrum, args...)
    _xlabel, _ylabel = _labels()
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
    _xlabel, _ylabel = _labels()
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (freq, S(freq), args...)
end

end #Module end
