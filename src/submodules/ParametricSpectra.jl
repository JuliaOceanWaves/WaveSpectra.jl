"""
Parametric omnidirectional spectra and directional spreading functions.

Currently only a few parametric ocean wave spectra and spreading functions are implemented.
The intention is to grow this library.
"""
module ParametricSpectra

using ..WaveSpectra: Dispersion, OmnidirectionalSpectrum, Spectrum, integrate, isdirection,
                     isspatial, istemporal, uconvert, unit, m, rad, s, θ₀
using Unitful: Hz, Length, Quantity
using DimensionfulAngles: Angle

@inline function _check_frequency_axis(x, name::AbstractString)
    (istemporal(x) || isspatial(x)) && return nothing
    isdirection(x) &&
        throw(ArgumentError("`$name` must be spatial or temporal, not directional."))
    throw(ArgumentError("`$name` must be spatial or temporal."))
end

# omnidirectional spectra
"""
    spectrum_pierson_moskowitz(frequencies, significant_waveheight, peak_period;
        dispersion=Dispersion())

Returns an omnidirectional Pierson-Moskowitz spectrum for the provided frequency array.
The `frequencies` and `peak_period` can be spatial or temporal and are converted
using `dispersion`.

Based on IEC TS 62600-2 ED2 Annex C.2 (2019).
"""
function spectrum_pierson_moskowitz(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        peak_period::Quantity;
        dispersion::Dispersion = Dispersion()
)
    _check_frequency_axis(frequencies, "frequencies")
    _check_frequency_axis(peak_period, "peak_period")

    n = length(frequencies)
    ind = (frequencies .≠ 0)
    uaxis = unit(eltype(frequencies))

    aₛ = significant_waveheight / 2
    fₚ = uconvert.(Hz, peak_period, dispersion)
    f̅ = uconvert.(Hz, frequencies[ind], dispersion)

    b = -(5 / 4) * (fₚ ./ f̅) .^ 4
    a = -b * aₛ^2 ./ f̅
    spec = zeros(n) * unit(eltype(a))
    spec[ind] = a .* exp.(b)

    spec = OmnidirectionalSpectrum(spec, f̅)
    spec = uconvert(uaxis, :axis, spec, dispersion)
    return spec
end

"""
    spectrum_jonswap(frequencies, significant_waveheight, peak_period;
        dispersion=Dispersion(), gamma=nothing)

Returns an omnidirectional JONSWAP spectrum for the provided frequency axis.
The `frequencies` and `peak_period` can be spatial or temporal and are converted
using `dispersion`.

Based on IEC TS 62600-2 ED2 Annex C.2 (2019).
"""
function spectrum_jonswap(
        frequencies::AbstractVector{<:Quantity},
        significant_waveheight::Length,
        peak_period::Quantity;
        dispersion::Dispersion = Dispersion(),
        gamma::Union{Number, Nothing} = nothing
)
    _check_frequency_axis(frequencies, "frequencies")
    _check_frequency_axis(peak_period, "peak_period")

    n = length(frequencies)
    ind = (frequencies .≠ 0)
    uaxis = unit(eltype(frequencies))

    aₛ = significant_waveheight / 2
    fₚ = uconvert.(Hz, peak_period, dispersion)
    f̅ = uconvert.(Hz, frequencies[ind], dispersion)
    γ = gamma

    σₐ = 0.07
    σb = 0.09
    σ = [f ≤ fₚ ? σₐ : σb for f in f̅]

    if isnothing(γ)
        hₛ = 2aₛ
        tₚ = uconvert.(s, peak_period, dispersion)
        val = tₚ / √(hₛ)
        uval = s * m^(-1 // 2)
        γ = (val ≤ 3.6uval) ? 5 : (
            (val > 5uval) ? 1 : (
            exp(5.75 - 1.15val)
        ))
    end

    b = -(5 / 4) * (fₚ ./ f̅) .^ 4
    a = -b * aₛ^2 ./ f̅
    spec = zeros(n) * unit(eltype(a))
    spec[ind] = a .* exp.(b)

    spec[ind] = spec .* (1 - 0.287log(γ)) .*
                γ .^ exp.(-((f̅ .- fₚ) .^ 2) ./ (2 * σ .^ 2 * fₚ^2))

    spec = OmnidirectionalSpectrum(spec, f̅)
    spec = uconvert(uaxis, :axis, spec, dispersion)
    return spec
end

# spread functions
"""
    spread_cartwright(directions, frequencies, mean_direction, spread; under90=false)

Cosine-squared directional spreading of Cartwright (1963).
If `under90=true`, values outside ±90° from `mean_direction` are set to zero before
normalization.
This spread function does not have any frequency dependency.
"""
function spread_cartwright(
        directions::AbstractVector{<:Angle},
        frequencies::AbstractVector{<:Quantity},
        mean_direction::Angle,
        spread::Angle,
        ;
        under90::Bool = false
)
    _check_frequency_axis(frequencies, "frequencies")

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
    spread_func = Spectrum(spread_func ./ norm_vec, f, θ̅)
    return spread_func
end

end
