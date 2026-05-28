using Test
using Unitful: uconvert, m, Hz
using WaveSpectra
using WaveSpectra.DispersionRelations: gravitywaves_deepwater
using WaveSpectra.ParametricSpectra: spread_cartwright
using WaveSpectra.Moments

f = (1:3) * 1.0Hz
spec = OmnidirectionalSpectrum([2.0, 3.0, 4.0] * (m^2 / Hz), f)
dispersion = gravitywaves_deepwater()

# Trapezoidal integration
m₋₁ = (19 / 6) * (m^2 / Hz)  # (2/1/2 + 3/2 + 4/3/2)m^2/Hz^2 *1Hz = 19/6 m^2/Hz
m₀ = 6.0 * m^2  # (2/2 + 3 + 4/2)m^2/Hz * 1Hz = 6m^2
m₁ = 13.0 * m^2 * Hz  # (2*1/2 + 3*2 + 4*3/2)m^2 * 1Hz = 13m^2Hz
m₂ = 31.0 * m^2 * Hz^2  # (2*1*1/2 + 3*2*2 + 4*3*3/2)m^2Hz * 1Hz = (31)m^2Hz^2

@test moment(spec, -1) ≈ m₋₁
@test moment(spec, 0) ≈ m₀
@test moment(spec, 1) ≈ m₁
@test moment(spec, 2) ≈ m₂

@test significant_waveheight(spec) ≈ 4 * sqrt(m₀)
@test energy_frequency(spec) ≈ m₀ / m₋₁
@test mean_frequency(spec) ≈ m₁ / m₀
@test zero_crossing_frequency(spec) ≈ sqrt(m₂ / m₀)

@test mean_wavelength(spec) ≈ uconvert(m, energy_frequency(spec), dispersion)
@test mean_wavelength(spec, mean_frequency) ≈ uconvert(m, mean_frequency(spec), dispersion)

@test steepness(spec) ≈ significant_waveheight(spec) / mean_wavelength(spec)
@test steepness(spec, mean_frequency) ≈
      significant_waveheight(spec) / mean_wavelength(spec, mean_frequency)

m₀z = 7.0 * m^2
m₁z = 14.0 * m^2 * Hz
m₂z = 32.0 * m^2 * Hz^2

@test moment(spec, 0; include_zero = true) ≈ m₀z
@test moment(spec, 1; include_zero = true) ≈ m₁z
@test moment(spec, 2; include_zero = true) ≈ m₂z
@test_throws ArgumentError moment(spec, -1; include_zero = true)
@test significant_waveheight(spec; include_zero = true) ≈ 4 * sqrt(m₀z)
@test mean_frequency(spec; include_zero = true) ≈ m₁z / m₀z
@test zero_crossing_frequency(spec; include_zero = true) ≈ sqrt(m₂z / m₀z)
@test_throws ArgumentError energy_frequency(spec; include_zero = true)
@test mean_wavelength(spec, mean_frequency; include_zero = true) ≈
      uconvert(m, mean_frequency(spec; include_zero = true), dispersion)
@test steepness(spec, mean_frequency; include_zero = true) ≈
      significant_waveheight(spec; include_zero = true) /
      mean_wavelength(spec, mean_frequency; include_zero = true)

# Rectangular
method = WaveSpectra.RectangularRule()
m₋₁r = (29 / 6) * (m^2 / Hz)  # (2/1 + 3/2 + 4/3)m^2/Hz^2 *1Hz = 29/6 m^2/Hz
m₀r = 9.0 * m^2  # (2 + 3 + 4)m^2/Hz * 1Hz = 9m^2
m₁r = 20.0 * m^2 * Hz  # (2*1 + 3*2 + 4*3)m^2 * 1Hz = 20m^2Hz
m₂r = 50.0 * m^2 * Hz^2  # (2*1*1 + 3*2*2 + 4*3*3)m^2Hz * 1Hz = 50m^2Hz^2

@test moment(spec, -1; method) ≈ m₋₁r
@test moment(spec, 0; method) ≈ m₀r
@test moment(spec, 1; method) ≈ m₁r
@test moment(spec, 2; method) ≈ m₂r

@test significant_waveheight(spec; method) ≈ 4 * sqrt(m₀r)
@test energy_frequency(spec; method) ≈ m₀r / m₋₁r
@test mean_frequency(spec; method) ≈ m₁r / m₀r
@test zero_crossing_frequency(spec; method) ≈ sqrt(m₂r / m₀r)

@test mean_wavelength(spec; method) ≈
      uconvert(m, energy_frequency(spec; method), dispersion)
@test mean_wavelength(spec, mean_frequency; method) ≈
      uconvert(m, mean_frequency(spec; method), dispersion)

@test steepness(spec; method) ≈
      significant_waveheight(spec; method) / mean_wavelength(spec; method)
@test steepness(spec, mean_frequency; method) ≈
      significant_waveheight(spec; method) / mean_wavelength(spec, mean_frequency; method)

# Mean direction
mean_dir = 60°
spread = 20°
frequencies = collect(range(0.0, 1.0, length = 1024))[2:end] * Hz
directions = collect(range(0.0, 2π, length = 361)) * rad
spread_func = spread_cartwright(directions, frequencies, mean_dir, spread)
@test mean_direction(spread_func)≈mean_dir atol=0.5°
