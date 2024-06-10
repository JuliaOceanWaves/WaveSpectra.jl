using Unitful: Quantity, Frequency
# Struct
struct DiscreteOmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity, D, N} <: AbstractArray{TS, N}
    value::AbstractVecOrMat{<:Quantity}
    frequency::AbstractVector{<:Quantity}
    function DiscreteOmnidirectionalSpectrum(value::AbstractVecOrMat{<:Quantity}, frequency::AbstractVector{<:Quantity}; density::Bool=true)
        # Parameters
        N=ndims(value)
        D=density
        TF=eltype(frequency)
        TS=eltype(value)
        # Check Parameters
        (dimension(TF) ∉ _frequency_dims) && error("invalid frequency dimensions")
        # Check Arguments
        (size(value)[1] ≠ length(frequency)) && error("'frequency' and first dimension of 'value' array must be same size")
        return new{TS,TF,D,N}(value, frequency)
    end
end

_firstState() = 1
#Interface
Base.size(spectrum::DiscreteOmnidirectionalSpectrum) = size(spectrum.value)
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, i::Int) = spectrum.value[i]
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, I::Vararg{Int, 2}) = spectrum.value[I...]
Base.getindex(spectrum::DiscreteOmnidirectionalSpectrum, I) = [spectrum.value[i] for i in I]
Base.length(spectrum::DiscreteOmnidirectionalSpectrum) = length(spectrum.value)
Base.iterate(spectrum::DiscreteOmnidirectionalSpectrum) = isempty(spectrum.value) ? nothing : (spectrum.value[_firstState()], _firstState() + 1)
Base.iterate(spectrum::DiscreteOmnidirectionalSpectrum, state::Int) = (state <= length(spectrum.value)) ? (spectrum.value[state], state+1) : nothing
Base.firstindex(spectrum::DiscreteOmnidirectionalSpectrum) = _firstState()
Base.lastindex(spectrum::DiscreteOmnidirectionalSpectrum) = length(spectrum.value)
Base.eltype(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, D, N}) where {TS, TF, D, N} = TS

#TODO
#convert from continuous

# Plots recipes
function _labels(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, D, N}) where {TS, TF, D, N}
    x_label = "frequency"
    y_label = D ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(spectrum::DiscreteOmnidirectionalSpectrum{TS, TF, D, N}, args...) where {TS, TF, D, N}
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (spectrum.frequency, spectrum.value, args...)
end
