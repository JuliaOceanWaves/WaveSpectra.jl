# Quick Start

`WaveSpectra` stores ocean wave spectra with explicit frequency, direction, and
unit information. The common workflow is to build an omnidirectional or
directional spectrum, integrate it, and pass summary metrics to another model.

```julia
using DimensionfulAngles.DefaultSymbols
using Unitful
using WaveSpectra

f = range(0.05, 0.5; length = 32) .* u"Hz"
theta = range(0, 350; length = 36) .* °
S = abs.(randn(length(f), length(theta))) .* u"m^2/Hz" / °

spec = Spectrum(S, f, theta)
omni = OmnidirectionalSpectrum(spec)
energy = integrate(omni)
```

Plot recipes are available through `Plots`. Makie plotting support is loaded only
when Makie is present, so core simulation and optimization environments do not
need Makie as a hard dependency.
