# Quickstart

## Intro

Leveraging [AxisArrays.jl](https://juliaarrays.github.io/AxisArrays.jl/latest/)
we store additional information such as the 
[Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/)
and [DimensionfulAngles.jl](https://juliaoceanwaves.github.io/DimensionfulAngles.jl/stable/)
quantities. This helps preserve the units across transformations and ensure that all
operations respect the units of the data. DimensionfulAngles extends the functionality of
Unitful quantities to angles and this package leverages that functionality to ensure 
consistency when dealing with wave spectra. Below is an example of the former; normally 
Unitful would not correctly handle the conversion from _degrees_ to _radians_ as they are 
both considered unitless.

```julia
using WaveSpectra, AxisArrays, DimensionfulAngles

julia> f = (6:6:18) * Hz; Θ = (120:120:360) * °;

julia> A = AxisArray(ones(Float64, (3, 3)) * m^2/Hz/°, f, Θ);

julia> S1 = Spectrum(A)

3×3 Spectrum{m² °⁻¹ Hz⁻¹}{Hz}{°}
Spectral density of the quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m² °⁻¹ Hz⁻¹):
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0

julia> S2 = uconvert(m^2, Hz, rad, S1)

3×3 Spectrum{m² Hz⁻¹ rad⁻¹}{Hz}{rad}
Spectral density of the quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (rad)
and data(m² Hz⁻¹ rad⁻¹):
 57.29577951308232  57.29577951308232  57.29577951308232
 57.29577951308232  57.29577951308232  57.29577951308232
 57.29577951308232  57.29577951308232  57.29577951308232
```

## Conversions

In scenarios where the user is working with multiple spectra, this package will handle
conversions when appropriate:

```julia
julia> S1

3×3 Spectrum{m² °⁻¹ Hz⁻¹}{Hz}{°}
Spectral density of the quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m² °⁻¹ Hz⁻¹):
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0

julia> S2

3×3 Spectrum{m² Hz⁻¹ rad⁻¹}{Hz}{rad}
Spectral density of the quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (rad)
and data(m² Hz⁻¹ rad⁻¹):
 57.29577951308232  57.29577951308232  57.29577951308232
 57.29577951308232  57.29577951308232  57.29577951308232
 57.29577951308232  57.29577951308232  57.29577951308232

julia> S1 + S2

3×3 Spectrum{m² s rad⁻¹}{Hz}{°}
Spectral density of the quantity (° Hz m² s rad⁻¹) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m² s rad⁻¹):
 114.59155902616465  114.59155902616465  114.59155902616465
 114.59155902616465  114.59155902616465  114.59155902616465
 114.59155902616465  114.59155902616465  114.59155902616465
```

and will notify the user when otherwise incompatible.

```julia
julia> f = (6:6:18) * Hz; Θ = (120:120:360) * °;

julia> A3 = AxisArray(ones(Float64, (3, 3)) * m^3/Hz/°, f, Θ); # Cubic meters

julia> S3 = Spectrum(A3)

3×3 Spectrum{m³ °⁻¹ Hz⁻¹}{Hz}{°}
Spectral density of the quantity (m³) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m³ °⁻¹ Hz⁻¹):
 1.0  1.0  1.0
 1.0  1.0  1.0
 1.0  1.0  1.0

julia> S1 + S3

ERROR: DimensionError: 1.0 m² °⁻¹ Hz⁻¹ and 1.0 m³ °⁻¹ Hz⁻¹ are not dimensionally compatible.
```

The accepted spectral-variables types are temporal/spatial, frequency/period, and 
linear/angular combinations. Represented as a diagram [here](@ref spectral_var_cube).

## Characterization

This package also includes a few functions to characterize ocean wave spectra both for 
[directional spectra](@ref dir_spectra) and [omnidirectional spectra](@ref omnidir_spectra).

```julia
julia> f = (6:6:18) * Hz; Θ = (120:120:360) * °;

julia> A = AxisArray(Float64.([x+y for x in 0:2, y in 0:2]) * m^2/Hz/°, f, Θ);

julia> S = Spectrum(A)
3×3 Spectrum{m² °⁻¹ Hz⁻¹}{Hz}{°}
Spectral density of the quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m² °⁻¹ Hz⁻¹):
 0.0  1.0  2.0
 1.0  2.0  3.0
 2.0  3.0  4.0

julia> WaveSpectra.Moments.mean_direction(S)
-79.10660535086907°
```

```julia
julia> S_omni = OmnidirectionalSpectrum((1.0:3.0) .* m^2/Hz, (1.0:3.0) .* Hz)
3-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.0
 2.0
 3.0

julia> WaveSpectra.Moments.significant_waveheight(S_omni)
8.0 m

julia> WaveSpectra.Moments.energy_frequency(S_omni)
2.0 Hz
```

## Other Functions

This package also includes a few functions to create a parametric spectra.

```julia
julia> f = (1.0:10.0) .* Hz;

julia> WaveSpectra.ParametricSpectra.spectrum_pierson_moskowitz(f, 8.0 * m, 4.0 * Hz)
10-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.3423912370794646e-72
 0.001701611389553719
 1.3421276585162367
 1.376317426592531
 0.6727667924496734
 0.3121398175430826
 0.1535891690462104
 0.08116740160851418
 0.045764352079589225
 0.02727015323859412
```