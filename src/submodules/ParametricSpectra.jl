"""
Parametric omnidirectional spectra and directional spread functions.

Currently only a few parametric ocean wave spectra and spread functions are implemented.
The intention is to grow this library.
"""
module ParametricSpectra

using ..WaveSpectra: Dispersion, OmnidirectionalSpectrum, Spectrum, integrate, isdirection,
                     isspatial, istemporal, uconvert, unit, m, rad, s, θ₀
using Unitful: Hz, Length, Quantity
using DimensionfulAngles: Angle

export spectrum_jonswap, spectrum_pierson_moskowitz, spread_cartwright

# omnidirectional spectra
"""
    spectrum_pierson_moskowitz(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        energy_frequency::Quantity;
        dispersion::Dispersion=Dispersion(),
        peak_frequency::Bool=false
    )

Omnidirectional Pierson-Moskowitz spectrum for the provided frequency array.

The `frequencies` and `energy_frequency` can be of any frequency type
(e.g., period, wavenumber, etc) and are converted to frequency (Hz).
The spectrum is created using the frequency formulation and then converted to the frequency
type of the inputs.

You can provide the peak frequency instead of the energy frequency by passing
`peak_frequency=true`.
The peak frequency is defined as the peak of `S(f)` where `f` is the frequency (e.g. in Hz),
regardless of the dimensions of `frequencies`.
See documentation for more details.

Based on IEC TS 62600-2 ED2 Annex C.2 (2019).
Peak frequency and energy frequency relationship from ITTC Specialist Committee on Waves
(2002).
"""
function spectrum_pierson_moskowitz(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        energy_frequency::Quantity;
        dispersion::Dispersion = Dispersion(),
        peak_frequency::Bool = false
)
    n = length(frequencies)
    uaxis = unit(eltype(frequencies))
    frequencies = uconvert.(Hz, frequencies, dispersion)
    ind = (frequencies .≠ 0Hz)
    energy_frequency = uconvert(Hz, energy_frequency, dispersion)

    aₛ = significant_waveheight / 2  # significant wave amplitude
    fₚ = peak_frequency ? energy_frequency : (0.858 * energy_frequency)  # peak frequency
    f̅ = frequencies[ind]

    b = -(5 / 4) * (fₚ ./ f̅) .^ 4
    a = -b * aₛ^2 ./ f̅
    spec = zeros(n) * unit(eltype(a))
    spec[ind] = a .* exp.(b)

    spec = OmnidirectionalSpectrum(spec, frequencies)
    return uconvert(uaxis, :axis, spec, dispersion)
end

"""
    spectrum_jonswap(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        energy_frequency::Quantity,
        gamma::Union{Number, Nothing};
        dispersion::Dispersion=Dispersion(),
        peak_frequency::Bool=false
    )

Omnidirectional JONSWAP spectrum for the provided frequency axis.

The `frequencies` and `energy_frequency` can be of any frequency type
(e.g., period, wavenumber, etc) and are converted to frequency (Hz).
The spectrum is created using the frequency formulation and then converted to the frequency
type of the inputs.

You can provide the peak frequency instead of the energy frequency by passing
`peak_frequency=true`.
The peak frequency is defined as the peak of `S(f)` where `f` is the frequency (e.g. in Hz),
regardless of the dimensions of `frequencies`.
See documentation for more details.

When peak frequency is used, `gamma` can be set to `gamma=nothing` in which case it will be
calculated based on the standard formulation.

Based on IEC TS 62600-2 ED2 Annex C.2 (2019).
Peak frequency and energy frequency relationship from ITTC Specialist Committee on Waves
(2002).
"""
function spectrum_jonswap(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        energy_frequency::Quantity,
        gamma::Union{Number, Nothing};
        dispersion::Dispersion = Dispersion(),
        peak_frequency::Bool = false
)
    n = length(frequencies)
    uaxis = unit(eltype(frequencies))
    frequencies = uconvert.(Hz, frequencies, dispersion)
    ind = (frequencies .≠ 0Hz)
    energy_frequency = uconvert(Hz, energy_frequency, dispersion)

    if (!peak_frequency && isnothing(gamma))
        throw(ArgumentError("`gamma` must be specified when `peak_frequency` is `false`"))
    end
    γ = gamma

    aₛ = significant_waveheight / 2  # significant wave amplitude
    c = peak_frequency ? 1 : 0.8255 + 0.03852γ - 0.005537γ^2 + 0.0003154γ^3
    #  NOTE: for γ=1, `c=0.0.8587984` which is larger than the PM value of `c=0.858`
    fₚ = c * energy_frequency  # peak frequency
    f̅ = frequencies[ind]

    σₐ = 0.07
    σb = 0.09
    σ = [f ≤ fₚ ? σₐ : σb for f in f̅]

    if isnothing(γ)
        hₛ = 2aₛ
        tₚ = uconvert.(s, 1 / fₚ, dispersion)
        val = tₚ / √(hₛ)
        uval = s * m^(-1 // 2)
        γ = (val ≤ 3.6uval) ? 5 : (
            (val > 5uval) ? 1 : (
            exp(5.75 - 1.15val)
        ))
    end

    # PM spectrum
    b = -(5 / 4) * (fₚ ./ f̅) .^ 4
    a = -b * aₛ^2 ./ f̅
    spec = zeros(n) * unit(eltype(a))
    spec[ind] = a .* exp.(b)

    # JONSWAP
    spec[ind] = spec[ind] .* (1 - 0.287log(γ)) .*
                γ .^ exp.(-((f̅ .- fₚ) .^ 2) ./ (2 * σ .^ 2 * fₚ^2))

    spec = OmnidirectionalSpectrum(spec, frequencies)
    return uconvert(uaxis, :axis, spec, dispersion)
end

# spread functions
"""
    spread_cartwright(
        directions::AbstractVector{<:Angle},
        frequencies::AbstractVector{<:Quantity},
        mean_direction::Angle,
        spread::Angle;
        under90::Bool=false
    )

Cosine-squared directional spread function of Cartwright (1963).

If `under90=true`, values outside ±90° from `mean_direction` are set to zero before
normalization.
This spread function does not depend on the spectral variable.

Code adapted from the `wavespectra` Python package.
"""
function spread_cartwright(
        directions::AbstractVector{<:Angle},
        frequencies::AbstractVector{<:Quantity},
        mean_direction::Angle,
        spread::Angle,
        ;
        under90::Bool = false
)
    θ̅, θₘ, f = directions, mean_direction, frequencies
    Δθ̅ = mod2pi.(θ̅ .- θₘ)
    ind = Δθ̅ .> π * rad
    Δθ̅[ind] .= 2π * rad .- Δθ̅[ind]

    s̃ = (2(θ₀ / spread)^2) - 1
    spread_func = (cos.(0.5Δθ̅)) .^ 2s̃
    spread_func = repeat(spread_func', length(f), 1)

    # zero-out directions ±90° from mean direction
    under90 && (spread_func[:, Δθ̅ .≥ π / 4 * rad] .= 0)

    # normalize
    spread_func = Spectrum(spread_func, f, θ̅)
    norm_vec = integrate(spread_func, :direction)
    spread_func = Spectrum(spread_func ./ norm_vec.data, f, θ̅)
    return spread_func
end

end
