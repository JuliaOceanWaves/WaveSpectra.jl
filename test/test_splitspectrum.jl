using Test
using Unitful
using DimensionfulAngles.DefaultSymbols
using WaveSpectra

f = [0.1, 0.2, 0.3] .* Hz
θ = range(0, 2π, length = 17) .* rad
omni = OmnidirectionalSpectrum([1.0, 2.0, 3.0] .* (m^2 / Hz), f)
spread = WaveSpectra.ParametricSpectra.spread_cartwright(θ, f, π / 3 * rad, π / 8 * rad)
spec = Spectrum(omni, spread)

# Omnidirectional spectrum of a polar spectrum
@test OmnidirectionalSpectrum(spec) isa OmnidirectionalSpectrum
@test OmnidirectionalSpectrum(spec) ≈ omni

# Spread-function of a polar spectrum
spread2 = spread_function(spec)
@test spread2 isa Spectrum
@test isspread(spread2)
@test spread2 ≈ spread

# Split a spectrum
omni2, spread3 = split_spectrum(spec)
@test omni2 ≈ omni
@test spread3 ≈ spread

# Rebuild the original spectrum from omnidirectional spectrum and spread function
spec2 = Spectrum(omni2, spread3)
@test spec2 isa Spectrum
@test spec2 ≈ spec

# Test not a spread function if it does not integrate to 1
badspread = Spectrum(fill(1.0, length(f), length(θ)) .* (rad^-1), f, θ)
@test !isspread(badspread)

# Throw errors if cartesian
cart_axis = collect((-1:1) .* (rad / m))
cart = Spectrum(ones(3, 3) .* (m^4 / rad^2), cart_axis, cart_axis)
@test !isspread(cart)
@test_throws ArgumentError OmnidirectionalSpectrum(cart)
@test_throws ArgumentError Spectrum(omni, cart)

# Throw error if frequency axes do not match
omni_bad_axis = OmnidirectionalSpectrum(omni.data, f .+ 0.01Hz)
@test_throws ArgumentError Spectrum(omni_bad_axis, spread)
