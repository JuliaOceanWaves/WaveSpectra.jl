#!/usr/bin/env julia
# Regular Julia script version of examples/ndbc.jl (Pluto-free), adapted for tests.
# Run from repo root: julia --project=test test/ndbc_script.jl

import Unitful
using DimensionfulAngles
using DimensionfulAngles.DefaultSymbols
using Test

include(joinpath(@__DIR__, "..", "src", "WaveSpectra.jl"))

# --- Settings ---
Nf = 20
Δf = 0.1Hz
f = (0:Nf-1) * Δf
Nθ = 36
Δθ = 360° / Nθ
θ = (0:Nθ-1) * Δθ

S = WaveSpectra.Spectrum(randn(Nf, Nθ)*m^2/Hz/°, f, θ)

@test size(S) == (Nf, Nθ)