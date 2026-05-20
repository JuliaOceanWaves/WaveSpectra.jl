using Test
using DimensionfulAngles: °ᵃ as °, radᵃ as rad
using WaveSpectra
using WaveSpectra.Moments: energy_frequency, mean_direction, significant_waveheight
using WaveSpectra.ParametricSpectra
using Unitful: Hz, m, s

frequencies = collect(range(0.0, 1.0, length = 1024))[2:end] * Hz
directions = collect(range(0.0, 2π, length = 361)) * rad
Hₛ = 2.3m
Tₑ = 9.5s
fₑ = 1 / Tₑ
mean_dir = 60°
spread = 20°

pm = spectrum_pierson_moskowitz(frequencies, Hₛ, fₑ)
pm_peak = spectrum_pierson_moskowitz(frequencies, Hₛ, 0.858fₑ; peak_frequency = true)
jonswap = spectrum_jonswap(frequencies, Hₛ, fₑ, 3.3)
jonswap_peak = spectrum_jonswap(frequencies, Hₛ, fₑ, nothing; peak_frequency = true)
spread_func = spread_cartwright(directions, frequencies, mean_dir, spread)
spread_func_under90 = spread_cartwright(
    directions, frequencies, mean_dir, spread; under90 = true
)
Δθ = mod2pi.(directions .- mean_dir)
Δθ[Δθ .> π * rad] .= 2π * rad .- Δθ[Δθ .> π * rad]
under90_mask = Δθ .>= (π / 4) * rad

# Hs
@test significant_waveheight(pm)≈Hₛ rtol=5e-3
@test significant_waveheight(jonswap)≈Hₛ rtol=5e-3
@test significant_waveheight(pm_peak)≈Hₛ rtol=5e-3
@test significant_waveheight(jonswap_peak)≈Hₛ rtol=5e-3

# Te
@test (1 / energy_frequency(pm))≈Tₑ rtol=5e-3
@test (1 / energy_frequency(jonswap))≈Tₑ rtol=5e-3

# spread function
@test isspread(spread_func)
@test isspread(spread_func_under90)
@test all(spread_func_under90[:, under90_mask] .== 0.0rad^-1)

# JONSWAP
js1 = spectrum_jonswap(frequencies, Hₛ, fₑ, 1)
@test js1≈pm rtol=5e-3 atol=1e-6 * m^2 / Hz
@test_throws ArgumentError spectrum_jonswap(frequencies, Hₛ, fₑ, nothing)

# spectrum from omnidirectional spectrum and spread function
directional = Spectrum(pm, spread_func)
omni = OmnidirectionalSpectrum(directional)
spread = spread_function(directional)
@test directional isa Spectrum
@test mean_direction(directional) ≈ mean_dir atol = 0.5°
@test omni≈pm rtol=5e-3 atol=1e-6 * m^2 / Hz
idx = .!isnan.(spread)
@test spread[idx]≈spread_func[idx] rtol=5e-3 atol=1e-6/rad
