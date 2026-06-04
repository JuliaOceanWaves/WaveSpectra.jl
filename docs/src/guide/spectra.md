# Spectra

Previously mentioned in the [Quickstart](@ref Quickstart) is the use of 
[AxisArrays.jl](https://juliaarrays.github.io/AxisArrays.jl/latest/)
for the structure of the data and both 
[Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/)
and [DimensionfulAngles.jl](https://juliaoceanwaves.github.io/DimensionfulAngles.jl/stable/)
to ensure the units of the data are respected in operations and conversions. The accepted 
spectral-variables types are temporal/spatial, frequency/period, and linear/angular 
combinations. Represented as a diagram [here](@ref spectral_var_cube).

## [Directional Spectra](@id directional_spectra)

You can construct a [Spectrum](@ref WaveSpectra.Spectrum) using an AxisArray matrix given 
that there are exactly two axes and both are spectral variables for a cartesian spectrum or 
only one is a direction for a polar spectrum. For the most control, one can pass in 
a 2-D Matrix and two separate axes with their respective units included.

```julia
julia> x1 = (1.0:4.0) .* m; x2 = (1.0:4.0) .* m;

julia> A = AxisArray(ones(Float64, (4,4)).*m, x1, x2);

julia> S = Spectrum(A)
4×4 Spectrum{m}{m}{m}
Spectral density of the quantity (m³) with cartesian coordinates:
  • Axis 1: Wavelength (m)
  • Axis 2: Wavelength (m)
and data(m):
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
```

```julia
julia> x = (1.0:4.0) .* m; Θ = (1.0:4.0) .* °;

julia> A = AxisArray(ones(Float64, (4,4)).*m, x, Θ);

julia> S = Spectrum(A)
4×4 Spectrum{m}{m}{°}
Spectral density of the quantity (° m²) with polar coordinates:
  • Axis 1: Wavelength (m)
  • Axis 2: Direction (°)
and data(m):
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
```

```julia
julia> x1 = (1.0:4.0) .* Hz; x2 = (1.0:4.0) .* °;

julia> S = Spectrum(ones(Float64, (4,4)) .* m^2, x1, x2)
4×4 Spectrum{m²}{Hz}{°}
Spectral density of the quantity (° Hz m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m²):
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
```

Please refer to the full syntax for each function [here](@ref directional_syntax).

## [Omnidirectional Spectra](@id omnidirectional_spectra)

You can construct an [OmnidirectionalSpectrum](@ref WaveSpectra.OmnidirectionalSpectrum)
using two vectors with units.

```julia
julia> S = (1.0:4.0) .* m^2/Hz; f = (1.0:4.0) .* Hz; 

julia> S = OmnidirectionalSpectrum(S, f)
4-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.0
 2.0
 3.0
 4.0
```

We also support the transformation of directional spectra into omnidirectional spectra. 

```julia
julia> x1 = (1.0:4.0) .* Hz; x2 = (1.0:4.0) .* °;

julia> S = Spectrum([x + y for x in 0:3, y in 1:4] .* m^2, x1, x2)
4×4 Spectrum{m²}{Hz}{°}
Spectral density of the quantity (° Hz m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m²):
 1  2  3  4
 2  3  4  5
 3  4  5  6
 4  5  6  7

julia> S_omni = OmnidirectionalSpectrum(S)
4-element OmnidirectionalSpectrum{° m²}{Hz}
Spectral density of the quantity (° Hz m²):
  • Axis: Frequency (Hz)
and data(° m²):
  7.5
 10.5
 13.5
 16.5

 julia> S_spread = spread_function(S)
4×4 Spectrum{°⁻¹}{Hz}{°}
Spectral density of the quantity (Hz) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(°⁻¹):
 0.13333333333333333  0.26666666666666666  0.4                  0.5333333333333333
 0.19047619047619047  0.2857142857142857   0.38095238095238093  0.47619047619047616
 0.2222222222222222   0.2962962962962963   0.37037037037037035  0.4444444444444444
 0.24242424242424243  0.30303030303030304  0.36363636363636365  0.42424242424242425
```

As well as the inverse through the use of the omnidirectional spectra and a directional 
spread function.

```julia
julia> S2 = Spectrum(S_omni, S_spread)
4×4 Spectrum{m²}{Hz}{°}
Spectral density of the quantity (° Hz m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m²):
 1.0  2.0  3.0  4.0
 2.0  3.0  4.0  5.0
 3.0  4.0  5.0  6.0
 4.0  5.0  6.0  7.0
```

Please refer to the full syntax for each function [here](@ref omnidirectional_syntax).

## Syntax

  - [Directional Spectra](@ref directional_syntax)
  - [Omnidirectional Spectra](@ref omnidirectional_syntax)
  - [Units](@ref unit_spectra_syntax)

### [Directional Spectra](@id directional_syntax)
```@autodocs; canonical=false
Modules = [WaveSpectra]
Filter = x -> (regex_match(r"Spectrum", x) && !regex_match(r"Omni", x))
```

### [Omnidirectional Spectra](@id omnidirectional_syntax)
```@autodocs; canonical=false
Modules = [WaveSpectra]
Filter = x -> (regex_match(r"Spectrum", x) && regex_match(r"Omni", x))
```
```@docs; canonical=false
WaveSpectra.spread_function
```

### [Unit Syntax](@id unit_spectra_syntax)

```@docs; canonical=false
WaveSpectra.uconvert
WaveSpectra.unit
```