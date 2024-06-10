
function _check_omnidirectional_dimensions(TS, TF)
    if !(isa(1*upreferred(TF()), Frequency) || (TF == typeof(NoDims)))
        error("parameter `TF` must be a `Unitful.Frequency` or `Unitful.NoDims`")
    end
<<<<<<< HEAD
    if (TF == typeof(NoDims)) && !(TS == typeof(NoDims))
        error("if `TF` is dimensionless `TS` must be too")
=======
end

# Constructors
function OmnidirectionalSpectrum(func::Function; dims_frequency::Dimensions=𝐓^-1, dims_value::Dimensions=𝐋^2*𝐓, density=true)
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, nothing)
end

function OmnidirectionalSpectrum(func::Function, freqs::AbstractVector{<:Number}; density=true)
    dims_frequency = dimension(eltype(freqs))
    dims_value = dimension(func(0*upreferred(eltype(freqs))))
    discrete = (
        frequency = freqs,
        value = func.(frequency)
    )
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, discrete)
end

function OmnidirectionalSpectrum(spectrum::OmnidirectionalSpectrum, freqs::AbstractVector{<:Number})
    Ts, Tf, density = typeof(spectrum).parameters
    discrete = (
            frequency = freqs,
            value = spectrum.func.(frequency)
        )
    return OmnidirectionalSpectrum{Ts, Tf, density}(spectrum.func, discrete)
end

function OmnidirectionalSpectrum(freqs::AbstractVector{<:Number}, vals::AbstractVector{<:Number}; density=true, interpolation::Function=linear_interpolation)
    dims_frequency = dimension(eltype(freqs))
    dims_value = dimension(eltype(vals))
    discrete = (
        frequency = freqs,
        value = vals
    )
    func = x->interpolation(freqs, vals; extrapolation_bc = 0*unit(eltype(vals)))(x)
    @assert func.(freqs) ≈ vals
    return OmnidirectionalSpectrum{typeof(dims_value), typeof(dims_frequency), density}(func, discrete)
end

# Interface
"""
    isdiscrete(spectrum::OmnidirectionalSpectrum)

Return whether the OmnidirectionalSpectrum uses discrete vectors.

# Examples
```jldoctest
a = 1
b = 2
a + b

# output

3
```
"""
isdiscrete(spectrum::OmnidirectionalSpectrum) = !isnothing(spectrum.discrete)
isdensity(spectrum::OmnidirectionalSpectrum) = typeof(spectrum) <: OmnidirectionalSpectrum{S, F, true} where {S,F}
isunitful(spectrum::OmnidirectionalSpectrum{Ts, Tf, D})  where {Ts, Tf, D} = isa(1*upreferred(Tf()), Frequency)
Unitful.unit(spectrum::OmnidirectionalSpectrum{Ts, Tf, D})  where {Ts, Tf, D} = isdiscrete(spectrum) ? unit(eltype(spectrum.discrete.value)) : unit(spectrum(0*upreferred(Tf())))
Unitful.dimension(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}) where {Ts, Tf, D} = Ts()
function quantity(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}) where {Ts, Tf, D}
    dims = Ts()*Tf()
    units = upreferred(dims)
    return dims, units
end

(spectrum::OmnidirectionalSpectrum)(freq::Number) =  spectrum.func(freq)

mutable struct DiscreteError <: Exception
    message::String
    DiscreteError() = new("continuous `OmnidirectionalSpectrum` instance was treated as discrete")
end
Base.showerror(io::IO, e::DiscreteError) = print(io, "DiscreteError: $(e.message)")

_firststate() = 1
Base.getindex(spectrum::OmnidirectionalSpectrum, i::Int) = isdiscrete(spectrum) ? spectrum.discrete.value[i] : throw(DiscreteError())
Base.getindex(spectrum::OmnidirectionalSpectrum, I) = isdiscrete(spectrum) ? [spectrum.discrete.value[i] for i in I] : throw(DiscreteError())
Base.firstindex(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? _firststate() : throw(DiscreteError())
Base.lastindex(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? length(spectrum.discrete.value) : throw(DiscreteError())
Base.length(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? length(spectrum.discrete.value) : throw(DiscreteError())
Base.size(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? size(spectrum.discrete.value) : throw(DiscreteError())
Base.eltype(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? eltype(spectrum.discrete.value) : throw(DiscreteError())
Base.iterate(spectrum::OmnidirectionalSpectrum) = isdiscrete(spectrum) ? (spectrum.discrete.value[firstindex(spectrum)], _firststate()+1) : throw(DiscreteError())
function Base.iterate(spectrum::OmnidirectionalSpectrum, state::Int)
    isdiscrete(spectrum) || throw(DiscreteError())
    return (state<=length(spectrum)) ? (spectrum.discrete.value[state], state+1) : nothing
end

# Plots recipes
function _labels(spectrum::OmnidirectionalSpectrum)
    x_label = "frequency"
    y_label = isdensity(spectrum) ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(spectrum::OmnidirectionalSpectrum{Ts, Tf, D}, args...) where {Ts, Tf, D}
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    if isdiscrete(spectrum)
        marker := :auto
        (spectrum.discrete.frequency, spectrum.discrete.value, args...)
    else
        isempty(args) && (args=(0*upreferred(Tf()), 10*upreferred(Tf())))
        (spectrum.func, args...)
>>>>>>> 08475a6c916900be849ac17d7e2b1d6f60e97cb5
    end
    return nothing
end

include("omnidirectional_spectra/continuous.jl")
