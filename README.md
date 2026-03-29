[![Test](https://github.com/JuliaOceanWaves/WaveSpectra.jl/actions/workflows/Test.yaml/badge.svg)](https://github.com/JuliaOceanWaves/WaveSpectra.jl/actions/workflows/Test.yaml)
[![codecov](https://codecov.io/github/JuliaOceanWaves/WaveSpectra.jl/graph/badge.svg?token=YX6TR1YKBJ)](https://codecov.io/github/JuliaOceanWaves/WaveSpectra.jl)
[![deps](https://juliahub.com/docs/General/WaveSpectra/stable/deps.svg)](https://juliahub.com/ui/Packages/General/WaveSpectra?t=2)

# WaveSpectra.jl
A unit-aware spectral representation for ocean waves or other interfacial waves.

Spectra can be represented in polar or cartesian coordinates.
The representation keeps track of both the spectral density and the axes values.
At least on  axis must have units of frequency, or alternatively period, angular frequency, wavenumber, etc.
The second axis has units of either frequency (cartesian) or angle (polar).
Units are based on [`Unitful.jl`](https://github.com/JuliaPhysics/Unitful.jl) and [`DimensionfulAngles.jl`](https://github.com/JuliaOceanWaves/DimensionfulAngles.jl) with *angles* treated as a dimension.

## Basic Usage
We can manually create a Spectrum.

```julia-repl
julia> using DimensionfulAngles.DefaultSymbols
julia> using CairoMakie
julia> using WaveSpectra

julia> Nf = 100;
julia> Δf = 0.005Hz;
julia> f = (1:Nf) * Δf
(0.005:0.005:0.5) Hz

julia> Nθ = 36;
julia> Δθ = 360° / Nθ;
julia> θ = (0:Nθ-1) * Δθ
(0.0:10.0:350.0)°

julia> S = Spectrum(randn(Nf, Nθ)*m^2/Hz/°, f, θ)
100×36 Spectrum{m² °⁻¹ Hz⁻¹}{Hz}{°}
Spectral density for Quantity (m²) with polar coordinates:
  • Axis 1: Frequency (Hz)
  • Axis 2: Direction (°)
and data(m² °⁻¹ Hz⁻¹):
...
```

We can also define a spectrum using parametric represenations.
We start with a parametric omnidirectional spectrum with a significant wave height of 1.2m and a peak period of 8s.

```julia-repl
julia> Hₛ = 1.2m;
julia> Tₚ = 8s;

julia> S_omni = ParametricSpectra.spectrum_pierson_moskowitz(f, Hₛ, Tₚ)
100-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density for Quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
...

julia> plot_spectrum(S_omni)
```
<img height="400" alt="display" src="https://github.com/user-attachments/assets/15539198-a36a-4d3d-892e-00b496b37a69" />

Next we create a parametric directional spread function and combine it with the omnidirectional spectrum to get our spectrum.
We use a cosine squared model with mean direction 20°.

```julia-repl
julia> θₘ = 30°;

julia> spread = 40°;

julia> D = ParametricSpectra.spread_cartwright(θ, f, θₘ, spread);

julia> S = Spectrum(S_omni, D);

julia> plot_spectrum(S)
```
<img height="400" alt="display" src="https://github.com/user-attachments/assets/a8225890-9d04-4f51-9624-d34f3e76b271" />

### Integration
You can integrate a spectrum using hte `integrate` function.
If an axis is specified it integrates along that axis.
For example the omnidirectional spectrum can be obtained from integration along the direction axis.

```julia-repl
julia> integrate(S)
100-element OmnidirectionalSpectrum{m² Hz⁻¹}{Hz}
Spectral density for Quantity (m²):
  • Axis: Frequency (Hz)
and data(m² Hz⁻¹):
...

julia> all(integrate(S, :direction) .≈ OmnidirectionalSpectrum(S))
true
```

We can also integrate along both axes to get the total variance or energy per unit area.

```julia-repl
julia> integrate(S)
0.08956154541883508 m^2

julia> using Unitful: gn as g

julia> ρ = 1025kg/m^3
julia> 0.5*ρ*g*integrate(S) |> J/m^2
450.12809880807976 J m^-2
```

### Change of Variables
We can also convert between types of frequency axis, for example, between frequency and angular frequency or period.
This is a change of variable that requires some care.
```julia-repl
julia> using Unitful: uconvert

julia> S_omni_period = uconvert(s, :axis, S_omni);

julia> plot_spectrum(S_omni_period; ax_properties=Dict(:limits=>((0, 20), (0, nothing))))
```
<img height="400" alt="display" src="https://github.com/user-attachments/assets/7d899ed8-0664-421a-81f3-e0b893c669ff" />

When converting between temporal and a spatial frequency we need to specify a [dispersion relation](https://en.wikipedia.org/wiki/Dispersion_(water_waves)).

```julia-repl
julia> dispersion = DispersionRelations.gravitywaves_deepwater();

julia> S_omni_wavenumber = uconvert(rad/m, :axis, S_omni, dispersion);
```

We can do change of variables for a full spectrum too.

```julia-repl
julia> S2 = uconvert(rad/m, :frequency, S, dispersion);

julia> S2 = uconvert(rad, :direction, S2);

julia> plot_spectrum(S2)
```
<img height="400" alt="display" src="https://github.com/user-attachments/assets/fc56130b-d679-4feb-bb8c-1d5d369b237d" />


## Contributing

Contributions are welcome! 🎊 Please see the [contribution guidelines](https://github.com/JuliaOceanWaves/.github/blob/main/CONTRIBUTING.md) for ways to contribute to the project.
