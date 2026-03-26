"""
Functions for spectral moments and moment-derived wave parameters.

This submodule provides helpers for integrating omnidirectional spectra into
commonly used bulk quantities such as energy frequency and significant wave height.
"""
module Moments

using ..WaveSpectra: AbstractOmnidirectionalSpectrum, Dispersion, OmnidirectionalSpectrum,
                     integrate, uconvert
using Unitful: Hz

"""
    moment(x::AbstractOmnidirectionalSpectrum, n::Integer)

Returns the spectral moment of order `n` for the omnidirectional spectrum `x`.
"""
function moment(x::AbstractOmnidirectionalSpectrum, n::Integer)
    integrate(OmnidirectionalSpectrum(x.data .* (x.axis .^ n), x.axis))
end

"""
    significant_waveheight(x::AbstractOmnidirectionalSpectrum)

Returns the significant wave height of the omnidirectional spectrum `x`,
computed as ``4 \\sqrt{m_0}``.

See also: [`moment`](@ref).
"""
significant_waveheight(x::AbstractOmnidirectionalSpectrum) = 4√(moment(x, 0))

"""
    energy_frequency(
        x::AbstractOmnidirectionalSpectrum;
        dispersion::Dispersion = Dispersion()
    )

Returns the energy frequency of the omnidirectional spectrum `x`, computed as
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

Returns the mean frequency of the omnidirectional spectrum `x`, computed as
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

Returns the zero-crossing frequency of the omnidirectional spectrum `x`,
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

end
