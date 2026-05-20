using Test
using WaveSpectra
using WaveSpectra.Moments: energy_frequency, significant_waveheight
using WaveSpectra.Shapes
using Unitful: Hz, m

# direct construction
shape_direct = OmnidirectionalSpectrumShape(
    [3.0, 2.0, 1.0], [3.0, 2.0, 1.0])
@test shape_direct.data == [1.0, 2.0, 3.0]
@test collect(shape_direct.axis) == [1.0, 2.0, 3.0]
@test WaveSpectra.axestypes(shape_direct) == :ndfrequency
@test WaveSpectra.axestypes([1.0, 2.0, 3.0]) == :ndfrequency
@test WaveSpectra.axestypes(Val(:shape)) == :ndfrequency

# from an omnidirectional spectrum
spec = OmnidirectionalSpectrum([1.0, 2.0, 1.0] * (m^2 / Hz), [0.1, 0.2, 0.3] * Hz)
shape = OmnidirectionalSpectrumShape(spec)
@test significant_waveheight(shape) ≈ 1.0
@test energy_frequency(shape) ≈ 1.0

# round trip
spec_roundtrip = OmnidirectionalSpectrum(
    shape, significant_waveheight(spec), energy_frequency(spec))
@test spec_roundtrip ≈ spec

# scaling
Hₛ = 3.0m
fₑ = 0.25Hz
scaled = scale(spec, Hₛ, fₑ)
@test scaled.axis ≈ spec.axis * energy_frequency(spec) / (fₑ)
@test scaled.data ≈
      spec.data * ((Hₛ)^2 / significant_waveheight(spec)^2) *
      (energy_frequency(spec) / fₑ)
@test OmnidirectionalSpectrumShape(scaled) ≈ shape
