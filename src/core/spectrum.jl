"""
    Spectrum(
        data::AbstractMatrix,
        axis1::AbstractVector{<:Quantity},
        axis2::AbstractVector{<:Quantity}
    )

Directional wave spectrum with physical-unit axes and data.
"""
struct Spectrum{
    TDAT,
    TAX1 <: AbstractVector{<:Quantity},
    TAX2 <: AbstractVector{<:Quantity}
} <: AbstractSpectrum{TDAT}
    data::Matrix{TDAT}
    axis1::TAX1
    axis2::TAX2
    coordinates::Symbol
    axestypes::Tuple{Symbol, Symbol}
    axesnames::Tuple{Symbol, Symbol}

    function Spectrum(
            data::AbstractMatrix{},
            axis1::AbstractVector{<:Quantity},
            axis2::AbstractVector{<:Quantity}
    )
        validated = validate_superposition(data, axis1, axis2)

        return new{
            eltype(validated.data),
            typeof(validated.axis1),
            typeof(validated.axis2)
        }(
            validated.data,
            validated.axis1,
            validated.axis2,
            validated.coordinates,
            validated.axestypes,
            validated.axesnames
        )
    end
end

superposition_unit_aliases(::AbstractSpectrum) = (:superposition, :spectrum)

"""
    unit(x::AbstractSpectrum, quantity::Symbol)
    unit(x::AbstractSpectrum)

Return spectrum, axis, or integral units for a directional spectrum.
"""
unit(x::AbstractSpectrum) = unit(x, :spectrum)

function Spectrum(x::AxisArray)
    axes = axisvalues(x)
    if length(axes) != 2
        throw(DimensionMismatch("Exactly two axes are required for 'Spectrum'."))
    end
    return Spectrum(x.data, axes...)
end
