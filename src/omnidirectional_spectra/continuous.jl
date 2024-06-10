
# Struct
struct OmnidirectionalSpectrum{TS<:Quantity, TF<:Quantity} <: Function
    func :: Function

    function OmnidirectionalSpectrum(func::Function, TS::DataType, TF::DataType)
        (dimension(TF) ∉ _frequency_dims) && error("invalid frequency dimensions")
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
function (spectrum::OmnidirectionalSpectrum{TS, TF})(frequency::TF) where TS
    return spectrum.func(frequency)
end

function (spectrum::OmnidirectionalSpectrum{TS, TF})(
        frequency::_frequency_quantity[dimension(TF)]
    ) where TS
    return spectrum.func(convert(TF, frequency))
end
