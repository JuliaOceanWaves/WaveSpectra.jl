"""
Functions for spectral moments and moment-derived wave parameters.

This submodule provides helpers for integrating omnidirectional spectra into
commonly used bulk quantities such as energy frequency and significant wave height.
"""
module Moments

using ..WaveSpectra: AbstractSampledIntegralAlgorithm, AbstractOmnidirectionalSpectrum,
                     AbstractSpectrum, Dispersion, TrapezoidalRule,
                     integrate, uconvert, rebuild_superposition, °
using ..WaveSpectra.DispersionRelations: gravitywaves_deepwater
using Unitful: Hz, gn as g, m

export energy_frequency, mean_direction, mean_frequency, mean_wavelength, moment,
       significant_waveheight, steepness, zero_crossing_frequency

"""
    moment(
        x::AbstractOmnidirectionalSpectrum,
        n::Integer;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Spectral moment of order `n` for the omnidirectional spectrum `x`.
"""
function moment(x::AbstractOmnidirectionalSpectrum, n::Integer;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    include_zero && (n < 0) &&
        throw(ArgumentError(
            "Including the zero spectral value is not supported for negative moments."))
    return integrate(
        rebuild_superposition(x, x.data .* (x.axis .^ n), x.axis);
        method,
        include_zero
    )
end

"""
    significant_waveheight(
        x::AbstractOmnidirectionalSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Significant wave height of the omnidirectional spectrum `x`,
computed as ``4 \\sqrt{m_0}``.

See also: [`moment`](@ref).
"""
function significant_waveheight(x::AbstractOmnidirectionalSpectrum;
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    return 4√(moment(x, 0; method, include_zero))
end

"""
    energy_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Energy frequency of the omnidirectional spectrum `x`, computed as
``m_0 / m_{-1}`` for `S(f)`.

The energy period can be obtained as ``1/energy_frequency(x)``.
The provided spectrum is first converted to the temporal-linear-frequency representation
`S(f)`, using `dispersion` if the spectral variable is spatial.

See also: [`moment`](@ref).
"""
function energy_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    xf = uconvert(Hz, :axis, x, dispersion)
    return moment(xf, 0; method, include_zero) / moment(xf, -1; method, include_zero)
end

"""
    mean_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Mean frequency of the omnidirectional spectrum `x`, computed as
``m_1 / m_0`` for `S(f)`.

The mean period can be obtained as ``1/mean_frequency(x)``.
The provided spectrum is first converted to the temporal-linear-frequency representation
`S(f)`, using `dispersion` if the spectral variable is spatial.

See also: [`moment`](@ref).
"""
function mean_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    xf = uconvert(Hz, :axis, x, dispersion)
    return moment(xf, 1; method, include_zero) / moment(xf, 0; method, include_zero)
end

"""
    zero_crossing_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Zero-crossing frequency of the omnidirectional spectrum `x`,
computed as ``\\sqrt{m_2 / m_0}`` for `S(f)`.

The zero-crossing period can be obtained as ``1/zero_crossing_frequency(x)``.
The provided spectrum is first converted to the temporal-linear-frequency representation
`S(f)`, using `dispersion` if the spectral variable is spatial.

See also: [`moment`](@ref).
"""
function zero_crossing_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion(),
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    xf = uconvert(Hz, :axis, x, dispersion)
    return √(moment(xf, 2; method, include_zero) / moment(xf, 0; method, include_zero))
end

"""
    mean_wavelength(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater();
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Mean wavelength of the omnidirectional spectrum `x`.

By default, the characteristic spectral variable used is the energy frequency and the
dispersion relation is the deep-water gravity wave dispersion.

See also: [`energy_frequency`](@ref), [`mean_frequency`](@ref).
"""
function mean_wavelength(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater();
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    if include_zero
        f = mean_frequency_function(x; dispersion, method, include_zero)
    else
        f = mean_frequency_function(x; dispersion, method)
    end
    return uconvert(m, f, dispersion)
end

"""
    steepness(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater();
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule()
    )

Characteristic wave steepness of the omnidirectional spectrum `x`.

This is computed as ``H_s / \\lambda`` where `H_s` is the significant wave height and
`\\lambda` is a caharacteristic wavelength.
The characteristic wavelength is defined by the characteristic spectral variable and the
dispersion relation.
By default, the characteristic spectral variable is the energy frequency and the
dispersion relation is the deep-water gravity wave dispersion.

See also: [`significant_waveheight`](@ref), [`mean_wavelength`](@ref).
"""
function steepness(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater();
        include_zero::Bool = false,
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    Hₛ = significant_waveheight(x; method, include_zero)
    λ = mean_wavelength(x, mean_frequency_function, dispersion; method, include_zero)
    return Hₛ / λ
end

# Directional
function mean_direction(
        x::AbstractSpectrum;
        method::AbstractSampledIntegralAlgorithm = TrapezoidalRule())
    y = integrate(rebuild_superposition(x, x.data .* sin.(x.axis2)', x.axis1, x.axis2); method)
    z = integrate(rebuild_superposition(x, x.data .* cos.(x.axis2)', x.axis1, x.axis2); method)
    return atan(°, y, z)
end

end
