"""
Functions for spectral shape and shape parameters.

Includes normalizing, non-dimensionalizing, or re-scaling spectra.
"""
module Shapes

using ..WaveSpectra: AbstractOmnidirectionalSpectrum, Dispersion, OmnidirectionalSpectrum,
                     uconvert, _check_typeconsistency, _ensure_increasing_axis
using ..WaveSpectra.Moments: energy_frequency, moment, significant_waveheight
using Unitful: Hz, Length, NoUnits, Quantity, unit

import ..WaveSpectra: axestypes
import ..WaveSpectra.Moments: energy_frequency

export OmnidirectionalSpectrumShape, scale

# non-dimensional normalized spectrum
"""
    OmnidirectionalSpectrumShape(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion()
    )
    OmnidirectionalSpectrumShape(
        data::AbstractVector,
        axis::AbstractVector
    )

Spectral shape of an omnidirectional wave spectrum.

The spectral shape is normalized to have a unit significant wave height and unit energy
period.
The data and its axis are dimensionless.

See also: [`OmnidirectionalSpectrum`](@ref), [`significant_waveheight`](@ref),
[`energy_frequency`](@ref).
"""
struct OmnidirectionalSpectrumShape{TDAT, TAX} <: AbstractOmnidirectionalSpectrum{TDAT}
    data::Vector{TDAT}
    axis::TAX
    axistype::Symbol
    axisname::Symbol

    function OmnidirectionalSpectrumShape(
            data::AbstractVector,
            axis::AbstractVector
    )
        # checks
        @assert length(data)==length(axis) "Data and axis lengths do not match!"
        _check_typeconsistency(data)
        _check_typeconsistency(axis)
        data, axis = _ensure_increasing_axis(data, axis)
        data = data .|> NoUnits
        axis = axis .|> NoUnits

        # assign axis type and name
        axistype = axisname = :ndfrequency

        return new{eltype(data), typeof(axis)}(data, axis, axistype, axisname)
    end
end

function OmnidirectionalSpectrumShape(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion())
    f = uconvert.(Hz, x.axis, dispersion)
    fₑ = energy_frequency(x; dispersion)
    Hₛ = significant_waveheight(x)
    axis = f / fₑ
    data = x.data * fₑ / Hₛ^2
    return OmnidirectionalSpectrumShape(data, axis)
end

# dimensional omnidirectional spectrum starting from a shape
"""
    OmnidirectionalSpectrum(
            x::OmnidirectionalSpectrumShape,
            waveheight::Length,
            frequency::Quantity;
            dispersion::Dispersion = Dispersion()
    )

Scale a spectral shape to have a specified significant wave height and energy frequency.

Returns an `OmnidirectionalSpectrum`.

See also: [`OmnidirectionalSpectrumShape`](@ref), [`significant_waveheight`](@ref),
[`energy_frequency`](@ref).
"""
function OmnidirectionalSpectrum(
        x::OmnidirectionalSpectrumShape,
        waveheight::Length,
        frequency::Quantity;
        dispersion::Dispersion = Dispersion()
)
    fₑ = uconvert(Hz, frequency, dispersion)
    Hₛ = waveheight
    data = x.data * Hₛ^2 / fₑ
    axis = uconvert.(unit(frequency), x.axis * fₑ, dispersion)
    return OmnidirectionalSpectrum(data, axis)
end

# scaling
"""
Scale a spectrum to have a specified significant wave height and energy frequency while preserving its shape.

Returns an `OmnidirectionalSpectrum`.

See also: [`significant_waveheight`](@ref), [`energy_frequency`](@ref).
"""
function scale(
        x::OmnidirectionalSpectrum,
        waveheight::Length,
        frequency::Quantity;
        dispersion::Dispersion = Dispersion()
)
    Hₛ₀ = significant_waveheight(x)
    fₑ₀ = energy_frequency(x; dispersion)
    Hₛ = waveheight
    fₑ = uconvert(Hz, frequency, dispersion)
    axis = x.axis * fₑ₀ / fₑ
    data = x.data * (Hₛ^2 / Hₛ₀^2) * (fₑ₀ / fₑ)
    return OmnidirectionalSpectrum(data, axis)
end

# other functions dispatch
energy_frequency(x::OmnidirectionalSpectrumShape) = moment(x, 0) / moment(x, -1)

axestypes(::Val{:shape}) = :ndfrequency
axestypes(::OmnidirectionalSpectrumShape) = axestypes(Val(:shape))
axestypes(x::AbstractVector{<:Real}) = axestypes(Val(:shape))

end
