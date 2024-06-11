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

# Call Methods
function DiscreteOmnidirectionalSpectrum(func::Function, frequency::AbstractVector{<:Quantity}; density::Bool=true)
    value = func.(frequency)
    return DiscreteOmnidirectionalSpectrum(value, frequency; density)
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
Base.firstindex(::DiscreteOmnidirectionalSpectrum) = _firstState()
Base.lastindex(spectrum::DiscreteOmnidirectionalSpectrum) = length(spectrum.value)
Base.eltype(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = TS

# Unitful Interface
Unitful.unit(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = unit(TS)
Unitful.dimension(::DiscreteOmnidirectionalSpectrum{TS}) where {TS} = dimension(TS)
frequency_unit(::DiscreteOmnidirectionalSpectrum{TS, TF}) where {TS, TF} = unit(TF)
frequency_dimension(::DiscreteOmnidirectionalSpectrum{TS, TF}) where {TS, TF} = dimension(TF)

function quantity(::DiscreteOmnidirectionalSpectrum{TS, TF, D}) where {TS, TF, D}
    dimensions = dimension(TS) * dimension(TF)
    units = unit(TS) * unit(TF)
    return dimensions, units
end

# Plots recipes
function _labels(::DiscreteOmnidirectionalSpectrum{TS, TF, D}) where {TS, TF, D}
    x_label = "frequency"
    y_label = D ? "spectral density" : "discrete (integral) spectrum"
    return (x_label, y_label)
end

@recipe function f(spectrum::DiscreteOmnidirectionalSpectrum, args...)
    _xlabel, _ylabel = _labels(spectrum)
    xlabel --> _xlabel
    ylabel --> _ylabel
    marker := :auto
    (spectrum.frequency, spectrum.value, args...)
end

function spectral_moment(spectrum::DiscreteOmnidirectionalSpectrum, n::Real=0; 
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    # There are no keyword arguments used to solve SampledIntegralProblems
    # https://docs.sciml.ai/Integrals/stable/basics/SampledIntegralProblem/
    sol = solve(SampledIntegralProblem(spectrum.value.*(spectrum.frequency .^float(n)), spectrum.frequency; dim=1), alg)
    (sol.retcode != ReturnCode.Success) && error("Solution unsuccessful with code $(sol.retcode)")
    return upreferred.(sol.u)
end

function integrate(spectrum::DiscreteOmnidirectionalSpectrum; 
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    return spectral_moment(spectrum, 0; alg)
end

function energy_period(spectrum::DiscreteOmnidirectionalSpectrum; 
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_n1 = spectral_moment(spectrum, -1; alg)
    m_0 = spectral_moment(spectrum, 0; alg)
    return m_n1 ./ m_0
end

function significant_waveheight(spectrum::DiscreteOmnidirectionalSpectrum; 
        alg::AbstractIntegralAlgorithm=TrapezoidalRule())
    m_0 = spectral_moment(spectrum, 0; alg)
    return @. 4*√m_0
end