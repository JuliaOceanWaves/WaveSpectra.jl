module WaveSpectra

export WaveSpectrum, isdiscrete, isunitful, quantity
# export integrate, ∫, spectral_moment, energy_period, significant_waveheight

using Interpolations: linear_interpolation
using Plots: @recipe
using Unitful: Frequency, Hz, upreferred
import Unitful: unit, dimension, Quantity
using Integrals: IntegralProblem, solve, QuadGKJL
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode
import Base: Exception, showerror


# Struct definition
struct WaveSpectrum{D, U}
    func::Union{Nothing, Function}
    frequency::Union{Nothing, AbstractVector{<:Number}}
    value::Union{Nothing, AbstractVector{<:Number}}

    function WaveSpectrum{D,U}(func, frequency, value) where {D,U}
        # check parameters
        if !(isa(D, Bool) && isa(U, Bool))
            error("parameters must be booleans")
        end
        # check arguments
        option_1 = !isnothing(func) && isnothing(frequency) && isnothing(value) # T F F
        option_2 = !isnothing(func) && !isnothing(frequency) && isnothing(value) # T T F
        option_3 = isnothing(func) && !isnothing(frequency) && !isnothing(value) # F T T
        function check_frequency(frequency)
            (eltype(frequency)<:Frequency) || error("unknown units for `frequency` vector")
        end
        if !(option_1 || option_2 || option_3)
            error("bad combination of input arguments")
        end
        if option_3 && (size(frequency) !== size(value))
            error("`frequency` and `value` arrays are not same size!")
        end
        (option_2||option_3) && U && check_frequency(frequency)
        # create parametric object
        option_2 && (value = func.(frequency))
        return new{D,U}(func, frequency, value)
    end
end

# Constructors
function WaveSpectrum(func::Function; density=true, units=true)
    # option 1
    return WaveSpectrum{density}{units}(func, nothing, nothing)
end

function WaveSpectrum(func::Function, freqs::AbstractVector{<:Number}; density=true)
    # option 2
    units = (eltype(freqs) <: Quantity) ? true : false
    return WaveSpectrum{density}{units}(func, freqs)
end

function WaveSpectrum(S::WaveSpectrum, freqs::AbstractVector{<:Number})
    # option 2b
    density, units = typeof(S).parameters
    return WaveSpectrum(S.func, freqs; density, units)
end

function WaveSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; density=true, interpolation::Function=linear_interpolation)
    # option 3
    units = (eltype(vals) <: Quantity) ? true : false
    func = x->interpolation(freqs, vals; extrapolation_bc = 0)(x)
    S = WaveSpectrum{density, units}(func, freqs, nothing)
    @assert (S.value  ≈ vals)
    return S
end

# Interface
isdiscrete(S::WaveSpectrum) = !isnothing(S.value) && !isnothing(S.frequency)
isdensity(S::WaveSpectrum) = typeof(S) <: WaveSpectrum{true, U} where U
isunitful(S::WaveSpectrum) = typeof(S) <: WaveSpectrum{D, true} where D
unit(S::WaveSpectrum) = unit(eltype(S.value))
dimension(S::WaveSpectrum) = dimension(eltype(S.value))
function quantity(S::WaveSpectrum{true,true})
    dims = dimension(S) * dimension(eltype(S.frequency))
    units = upreferred(unit(S) * unit(eltype(S.frequency)))
    return dims, units
end

(S::WaveSpectrum)(freq::Number) =  S.func(freq)

mutable struct DiscreteError <: Exception
    message::String
    DiscreteError() = new("continuous `WaveSpectrum` instance was treated as discrete")
end
Base.showerror(io::IO, e::DiscreteError) = print(io, "DiscreteError: $(e.message)")

_firststate() = 1
Base.getindex(S::WaveSpectrum, i::Int) = isdiscrete(S) ? S.value[i] : throw(DiscreteError())
Base.firstindex(S::WaveSpectrum) = isdiscrete(S) ? _firststate() : throw(DiscreteError())
Base.lastindex(S::WaveSpectrum) = isdiscrete(S) ? length(S.value) : throw(DiscreteError())
Base.length(S::WaveSpectrum) = isdiscrete(S) ? length(S.value) : throw(DiscreteError())
Base.size(S::WaveSpectrum) = isdiscrete(S) ? size(S.value) : throw(DiscreteError())
Base.eltype(S::WaveSpectrum) = isdiscrete(S) ? eltype(S.value) : throw(DiscreteError())
Base.iterate(S::WaveSpectrum) = isdiscrete(S) ? (S.value[firstindex(S)], _firststate()+1) : throw(DiscreteError())
function Base.iterate(S::WaveSpectrum, state::Int)
    isdiscrete(S) || throw(DiscreteError())
    return (state<=length(S)) ? (S.value[state], state+1) : nothing
end

# Plots recipes
function _labels(S::WaveSpectrum)
    x_label = "frequency"
    y_label = isdensity(S) ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(S::WaveSpectrum, args...)
    _xlabel, _ylabel = _labels(S)
    xlabel --> _xlabel
    ylabel --> _ylabel
    if isdiscrete(S)
        marker := :auto
        (S.frequency, S.value, args...)
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
    (freq, S.(freq), args...)
end

# include("statistics.jl")

end #Module end
