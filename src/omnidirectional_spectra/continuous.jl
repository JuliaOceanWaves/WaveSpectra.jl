
# Struct
struct OmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity} <: Function
    func :: Function

    function OmnidirectionalSpectrum(func::Function, TS::DataType, TF::DataType)
        (dimension(TF) ∉ _frequency_dims) && error("invalid frequency dimensions")
        (typeof(func(ones(TF)[])) ≠ TS) && error("invalid spectrum dimensions")
        return new{TS, TF}(func)
    end
end

# Constructors
function OmnidirectionalSpectrum(func::Function, TF::DataType=typeof(1.0Hz))
    TS = typeof(func(ones(TF)[]))
    return OmnidirectionalSpectrum(func, TS, TF)
end

function OmnidirectionalSpectrum(
        value::AbstractVector{<:Quantity}, frequency::AbstractVector{<:Quantity};
        interpolation::Function=linear_interpolation
    )
    TS = eltype(value)
    TF = eltype(frequency)
    func = x -> interpolation(frequency, value; extrapolation_bc = 0*unit(TS))(x)
    return OmnidirectionalSpectrum(func, TS, TF)
end

# Call Methods
function (spectrum::OmnidirectionalSpectrum{TS,TF})(frequency::Quantity) where {TS,TF}
    (typeof(frequency) ≠ TF) && (frequency=convert(TF, frequency))
    return spectrum.func(frequency)
end

# Unitful interface
Unitful.unit(spectrum::OmnidirectionalSpectrum{TS, TF})  where {TS, TF} = unit(TS)
Unitful.dimension(spectrum::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TS)
frequency_unit(spectrum::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = unit(TF)
frequency_dimension(spectrum::OmnidirectionalSpectrum{TS,TF}) where {TS,TF} = dimension(TF)

function quantity(spectrum::OmnidirectionalSpectrum{TS, TF}) where {TS, TF}
    dimensions = dimension(TS) * dimension(TF)
    units = unit(TS) * unit(TF)
    return dimensions, units
end

# Plots recipe
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
    end
    return nothing
end
