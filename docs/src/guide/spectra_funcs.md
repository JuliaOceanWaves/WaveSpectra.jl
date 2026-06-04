# [Spectra Functions](@id spectra_functions)

## Spectrum Methods

There are a few other functions worth mentioning that work directly with the Spectrum 
object. You can query the information of the axes of a spectrum.

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

julia> axesinfo(S)
(((:temporal, :linear, :frequency), 𝐓⁻¹), ((:direction,), 𝐀))
```

Previously mentioned when transforming between 
[omnidirectional and directional spectra](@ref omnidirectional_spectra), but this package 
contains a method for capturing the spread function and whether a given Spectrum is a spread
 function.

```julia
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

julia> S_omni, S_spread = split_spectrum(S)
(4×1-OmnidirectionalSpectrum{° m²}{Hz}, 4×4-Spectrum{°⁻¹}{Hz}{°})

julia> isspread(S_spread)
true
```

Last functions worth mentioning transform directional spectra between polar and cartesian
coordinates.

```julia
julia> polar_to_cartesian(S)
2-dimensional AxisArray{Any,2,...} with axes:
    :point, 1:16
    :component, [:frequency_x, :frequency_y, :spectrum]
And data, a 16×3 Matrix{Any}:
 0.999848 Hz  0.0174524 Hz  1.0 m² rad Hz⁻¹
 1.9997 Hz    0.0349048 Hz  1.0 m² rad Hz⁻¹
 2.99954 Hz   0.0523572 Hz  1.0 m² rad Hz⁻¹
 3.99939 Hz   0.0698096 Hz  1.0 m² rad Hz⁻¹
 ⋮                          
 0.997564 Hz  0.0697565 Hz  4.0 m² rad Hz⁻¹
 1.99513 Hz   0.139513 Hz   2.5 m² rad Hz⁻¹
 2.99269 Hz   0.209269 Hz   2.0 m² rad Hz⁻¹
 3.99026 Hz   0.279026 Hz   1.75 m² rad Hz⁻¹

julia> x1 = (1.0:4.0) .* m; x2 = (1.0:4.0) .* m;

julia> S_cartesian = Spectrum([x + y for x in 0:3, y in 1:4] .* m^3, x1, x2)
4×4 Spectrum{m³}{m}{m}
Spectral density of the quantity (m⁵) with cartesian coordinates:
  • Axis 1: Wavelength (m)
  • Axis 2: Wavelength (m)
and data(m³):
 1  2  3  4
 2  3  4  5
 3  4  5  6
 4  5  6  7

julia> cartesian_to_polar(S_cartesian)
2-dimensional AxisArray{Any,2,...} with axes:
    :point, 1:16
    :component, [:wavelength, :direction, :spectrum]
And data, a 16×3 Matrix{Any}:
 1.41421 m  45.0°      1.41421 m⁴ rad⁻¹
 2.23607 m  26.5651°   4.47214 m⁴ rad⁻¹
 3.16228 m  18.4349°   9.48683 m⁴ rad⁻¹
 4.12311 m  14.0362°  16.4924 m⁴ rad⁻¹
 ⋮                    
 4.12311 m  75.9638°  16.4924 m⁴ rad⁻¹
 4.47214 m  63.4349°  22.3607 m⁴ rad⁻¹
 5.0 m      53.1301°  30.0 m⁴ rad⁻¹
 5.65685 m  45.0°     39.598 m⁴ rad⁻¹

```

Please refer to the full syntax for each function [here](@ref spectra_method_syntax).

## Dispersion Relations

The DispersionRelation.jl submodule contains the dispersion relations for linear interfacial
waves. Currently only gravity wave dispersion relations are implemented.

See also [https://en.wikipedia.org/wiki/Dispersion_(water_waves)](https://en.wikipedia.org/wiki/Dispersion_(water_waves)) 
for more information on dispersion.

```julia

julia> x = (1.0:3:10.0) * Hz;

julia> S = OmnidirectionalSpectrum(([1.0, 6.0, 3.0, 2.0])*m^2/Hz, x)
4-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density of the quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
 1.0
 6.0
 3.0
 2.0

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves_deepwater())
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.12420267319576646
 0.1863040097936497
 0.053229717083899904
 0.02484053463915329

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves_shallowwater(3m))
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.863258964143784
 5.179553784862704
 2.589776892431352
 1.726517928287568

julia> uconvert(rad/m, :axis, S, WaveSpectra.DispersionRelations.gravitywaves(3m))
4-element OmnidirectionalSpectrum{m³ rad⁻¹}{rad m⁻¹}
Spectral density of the quantity (m²):
  • Axis: Angular Wavenumber (rad m⁻¹)
and data(m³ rad⁻¹):
 0.12420267338189336
 0.1863040097936497
 0.053229717083899904
 0.02484053463915329

```

Please refer to the full syntax for each function [here](@ref dispersion_relation_syntax).

## Moments

The following examples are different functions used in literature for characterizing 
spectra. 

```julia
julia> x = (1.0:3:10.0) * Hz; S = OmnidirectionalSpectrum(([1.0, 6.0, 3.0, 2.0])*m^2, x);

julia> WaveSpectra.Moments.moment(S, 0)
31.5 m²

julia> WaveSpectra.Moments.energy_frequency(S)
4.152542372881356 Hz

julia> WaveSpectra.Moments.mean_frequency(S)
5.285714285714286 Hz

julia> WaveSpectra.Moments.mean_wavelength(S)
0.09051335476420994 m

julia> WaveSpectra.Moments.significant_waveheight(S)
22.44994432064365 m

julia> WaveSpectra.Moments.steepness(S)
248.02908232852977

julia> WaveSpectra.Moments.zero_crossing_frequency(S)
5.719640348333601 Hz

```

Please refer to the full syntax for each function [here](@ref moments_syntax).

## Utilities

Lastly these are functions that are either used in the above functions, but may be useful
for some users. 
```julia
julia> x = [1.0, 2.5, 3.1, 4.2, 5.22]*Hz; S = OmnidirectionalSpectrum([1.4, 5.2, 7.7, 3.1, 2].*m^2/Hz, x);

julia> evenspacing(S)
ERROR: ArgumentError: Vector `x` must be evenly spaced.

julia> S = OmnidirectionalSpectrum([1.4, 5.2, 7.7, 3.1, 2].*m^2/Hz, (1.0:5.0).*Hz);

julia> evenspacing(S)
(1.0 Hz, 1.0 Hz, 5)

julia> isevenlyspaced(S)
true
```

```julia
julia> integrate(S)
17.700000000000003 m²

julia> integrate(S, method=RectangularRule())
19.400000000000002 m²
```

Most notably is a function for plotting the polar spectrum.

```julia
julia> f = (1.0:10.0) * Hz; Θ = (0.0:20:360) * °;

julia> S = Spectrum([x + y for x in 0:9, y in 1:19], f, Θ)
10×19 Spectrum{1}{Hz}{°}
Spectral density of the quantity (° Hz) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(1):
  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
  2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20
  3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
  4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
  5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23
  6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
  7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25
  8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26
  9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27
 10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28

```

![polar_spectra](../assets/polar_spectra_plot.png)

Please refer to the full syntax for each function [here](@ref utilities_syntax).
 
## Syntax

  - [Spectrum Methods](@ref spectra_method_syntax)
  - [Dispersion Relation](@ref dispersion_relation_syntax)
  - [Moments](@ref moments_syntax)
  - [Utilities](@ref utilities_syntax)

### [Spectrum Methods](@id spectra_method_syntax)

```@docs; canonical=false
WaveSpectra.axesinfo
WaveSpectra.spread_function
WaveSpectra.split_spectrum
WaveSpectra.cartesian_to_polar
WaveSpectra.polar_to_cartesian
```

### [Dispersion Relation](@id dispersion_relation_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.DispersionRelations]
Order = [:function]
```

### [Moments](@id moments_syntax)

```@autodocs; canonical=false
Modules = [WaveSpectra.Moments]
Order = [:function]
```

### [Utilities](@id utilities_syntax)

```@docs; canonical=false
WaveSpectra.RectangularRule
WaveSpectra.evenspacing
WaveSpectra.integrate
WaveSpectra.isevenlyspaced
WaveSpectra.isspread
WaveSpectra.plot_spectrum!
WaveSpectra.plot_spectrum
```