module WaveSpectra

export DiscreteOmnidirectionalSpectrum, OmnidirectionalSpectrum
export integrate, frequency_dimension, frequency_unit, pm_spectrum, quantity, spectral_moment
export energy_period, significant_waveheight

# using Base: Base # extend: showerror
using Unitful: Unitful # extend: unit, dimension
using Interpolations: linear_interpolation
using Plots: @recipe
using Integrals: AbstractSampledIntegralAlgorithm, IntegralProblem, SampledIntegralProblem, TrapezoidalRule, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode


using Unitful: 𝐓, 𝐋, gn as g, Hz, m
using Unitful: dimension, unit, uconvert, upreferred
using Unitful: Frequency, Length, Time, Wavenumber
using Unitful: Dimensions, Frequency, NoDims, Quantity, Hz
using DimensionfulAngles: 𝐀
using DimensionfulAngles: AngularPeriod, AngularVelocity, AngularWavelength
using DimensionfulAngles: AngularWavenumber, Periodic

const _frequency_dims = [𝐓, 𝐓^-1, 𝐓*𝐀^-1, 𝐀*𝐓^-1, 𝐋, 𝐋^-1, 𝐋*𝐀^-1, 𝐀*𝐋^-1]

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
