
using DimensionfulAngles.DefaultSymbols
using Test
using WaveSpectra

Nf = 20
Δf = 0.1Hz
f = (0:(Nf - 1)) * Δf
Nθ = 36
Δθ = 360° / Nθ
θ = (0:(Nθ - 1)) * Δθ

S = Spectrum(randn(Nf, Nθ) * m^2 / Hz / °, f, θ)

@test size(S) == (Nf, Nθ)
