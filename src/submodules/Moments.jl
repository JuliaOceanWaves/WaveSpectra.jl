"""
Functions for spectral moments and moment-derived wave parameters.

This submodule provides helpers for integrating omnidirectional spectra into
commonly used bulk quantities such as energy frequency and significant wave height.
"""
module Moments

using ..WaveSpectra: AbstractOmnidirectionalSpectrum, Dispersion, integrate, uconvert,
                     _rebuild_spectrum
using ..WaveSpectra.DispersionRelations: gravitywaves_deepwater
using Unitful: Hz, gn as g, m

export moment

"""
    moment(x::AbstractOmnidirectionalSpectrum, n::Integer)

Spectral moment of order `n` for the omnidirectional spectrum `x`.
"""
function moment(x::AbstractOmnidirectionalSpectrum, n::Integer)
    integrate(_rebuild_spectrum(x, x.data .* (x.axis .^ n), x.axis))
end

"""
    significant_waveheight(x::AbstractOmnidirectionalSpectrum)

Significant wave height of the omnidirectional spectrum `x`,
computed as ``4 \\sqrt{m_0}``.

See also: [`moment`](@ref).
"""
significant_waveheight(x::AbstractOmnidirectionalSpectrum) = 4√(moment(x, 0))

"""
    energy_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion()
    )

Energy frequency of the omnidirectional spectrum `x`, computed as
``m_0 / m_{-1}`` for `S(f)`.

The energy period can be obtained as ``1/energy_frequency(x)``.
The provided spectrum is first converted to the frequency spectrum `S(f)`, using
`dispersion` if the spectrum is in terms of spatial frequency.

See also: [`moment`](@ref).
"""
function energy_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion())
    xf = uconvert(Hz, :axis, x, dispersion)
    return moment(xf, 0) / moment(xf, -1)
end

"""
    mean_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion()
    )

Mean frequency of the omnidirectional spectrum `x`, computed as
``m_1 / m_0`` for `S(f)`.

The mean period can be obtained as ``1/mean_frequency(x)``.
The provided spectrum is first converted to the frequency spectrum `S(f)`, using
`dispersion` if the spectrum is in terms of spatial frequency.

See also: [`moment`](@ref).
"""
function mean_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion())
    xf = uconvert(Hz, :axis, x, dispersion)
    return moment(xf, 1) / moment(xf, 0)
end

"""
    zero_crossing_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion()
    )

Zero-crossing frequency of the omnidirectional spectrum `x`,
computed as ``\\sqrt{m_2 / m_0}`` for `S(f)`.

The zero-crossing period can be obtained as ``1/zero_crossing_frequency(x)``.
The provided spectrum is first converted to the frequency spectrum `S(f)`, using
`dispersion` if the spectrum is in terms of spatial frequency.

See also: [`moment`](@ref).
"""
function zero_crossing_frequency(x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion())
    xf = uconvert(Hz, :axis, x, dispersion)
    return √(moment(xf, 2) / moment(xf, 0))
end

"""
    mean_wavelength(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater()
    )

Mean wavelength of the omnidirectional spectrum `x`.

By default, the characteristic frequency used is the energy frequency and the dispersion
relation is the deep-water gravity wave dispersion.

See also: [`energy_frequency`](@ref), [`mean_frequency`](@ref).
"""
function mean_wavelength(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater())
    f = mean_frequency_function(x; dispersion)
    return uconvert(m, f, dispersion)
end

"""
    steepness(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater()
    )

Characteristic wave steepness of the omnidirectional spectrum `x`.

This is computed as ``H_s / \\lambda`` where `H_s` is the significant wave height and
`\\lambda` is a caharacteristic wavelength.
The characteristic wavelength is defined by the characteristic frequency and the dispersion
relation.
By default, the representative frequency is the energy frequency and the
dispersion relation is the deep-water gravity wave dispersion.

See also: [`significant_waveheight`](@ref), [`mean_wavelength`](@ref).
"""
function steepness(
        x::AbstractOmnidirectionalSpectrum,
        mean_frequency_function::Function = energy_frequency,
        dispersion::Dispersion = gravitywaves_deepwater())
    Hₛ = significant_waveheight(x)
    λ = mean_wavelength(x, mean_frequency_function, dispersion)
    return Hₛ / λ
end

end
