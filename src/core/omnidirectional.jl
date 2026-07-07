"""
    OmnidirectionalSpectrum(
        data::AbstractVector,
        axis::AbstractVector{<:Quantity}
    )

One-dimensional omnidirectional wave spectrum with a physical-unit axis.
"""
struct OmnidirectionalSpectrum{TDAT, TAX} <: AbstractOmnidirectionalSpectrum{TDAT}
    data::Vector{TDAT}
    axis::TAX
    axistype::Symbol
    axisname::Symbol

    function OmnidirectionalSpectrum(
            data::AbstractVector{},
            axis::AbstractVector{<:Quantity}
    )
        validated = validate_superposition1d(data, axis)
        return new{eltype(validated.data), typeof(validated.axis)}(
            validated.data,
            validated.axis,
            validated.axistype,
            validated.axisname
        )
    end
end

superposition_unit_aliases(::AbstractOmnidirectionalSpectrum) = (:superposition, :spectrum)

"""
    unit(x::AbstractOmnidirectionalSpectrum, quantity::Symbol)
    unit(x::AbstractOmnidirectionalSpectrum)

Return spectrum, axis, or integral units for an omnidirectional spectrum.
"""
unit(x::AbstractOmnidirectionalSpectrum) = unit(x, :spectrum)

function OmnidirectionalSpectrum(x::AxisArray)
    axis = axisvalues(x)
    if length(axis) != 1
        throw(DimensionMismatch(
            "Exactly one axis is required for 'OmnidirectionalSpectrum'.")
        )
    end
    return OmnidirectionalSpectrum(x.data, axis[1])
end
