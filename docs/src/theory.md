# Theory

Wave spectra represent sea-state energy density as a function of frequency and,
for directional spectra, propagation direction. `OmnidirectionalSpectrum` stores
`S(f)`, while `Spectrum` stores `S(f, theta)` with unit-aware axes.

Integrals are computed over sampled axes using internal rectangular or trapezoid
rules. This is deliberate: buoy records, parametric spectra, and simulation
outputs are all discrete arrays, and preserving the sampled grid avoids hidden
resampling during resource coupling. The exported `RectangularRule` can be used
when bin-centered data should be treated as constant over each interval.

SIRENOpt uses `WaveSpectra` as a resource layer. Device-specific wave-energy
conversion should consume spectra or derived power-flux metrics without changing
the spectral data model itself.
