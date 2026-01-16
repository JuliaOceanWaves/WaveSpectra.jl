# WaveSpectra.jl

## Overview

WaveSpectra.jl provides a unit-aware tool for working with ocean wave spectra.
This package enables input, conversion, creation, and characterization 
for both polar, cartesian, and omnidirectional wave spectra.

## Quick start

```julia
using WaveSpectra

Nf = 20
Δf = 0.1Hz
f = (0:Nf-1) * Δf
Nθ = 36
Δθ = 360° / Nθ
θ = (0:Nθ-1) * Δθ

S = Spectrum(randn(Nf, Nθ)*m^2/Hz/°, f, θ)
```

## Detailed usage

### Creating and inputting spectra



### Spectrum conversion



### Spectrum characterization



## Reference

See the [API reference](@ref API) for all exported functions.

## Theory

