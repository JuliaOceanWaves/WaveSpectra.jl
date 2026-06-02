# Spectra

Previously mentioned in the [Quickstart](@ref Quickstart) is the use of 
[AxisArrays.jl](https://juliaarrays.github.io/AxisArrays.jl/latest/)
for the structure of the data and both 
[Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/)
and [DimensionfulAngles.jl](https://juliaoceanwaves.github.io/DimensionfulAngles.jl/stable/)
to ensure the units of the data are respected in operations and conversions. The accepted 
spectral-variables types are temporal/spatial, frequency/period, and linear/angular 
combinations. Represented as a diagram [here](@ref spectral_var_cube).

## [Directional Spectra](@id dir_spectra)

You can construct a [Spectrum](@ref WaveSpectra.Spectrum) using an AxisArray matrix given 
that there are exactly two axes and both are spectral variables for a cartesian spectrum or 
only one is a direction for a polar spectrum.  It is possible to freely convert between 
directional spectrum and omnidirectional given the omnidirectional function with a 
directional spread function. Lastly, for the most control, one can pass in a 2-D Matrix and 
two separate axes with their respective units included.


```@autodocs; canonical=false
Modules = [WaveSpectra]
Filter = x -> (regex_match(r"Spectrum", x) && !regex_match(r"Omni", x))
```

```julia
Sample Cartesian Spectrum
```

```julia
Sample Polar Spectrum
```

## [Omnidirectional Spectra](@id omnidir_spectra)

Brief description of Omnidirectional spectrum

```julia
Sample Omnidirectional Spectrum
```