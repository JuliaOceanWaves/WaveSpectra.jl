module WaveSpectra

export OmnidirectionalSpectrum, isdensity, isdiscrete, isunitful, quantity
export integrate, spectral_moment, energy_period, significant_waveheight

# using Base: Base # extend: showerror
using Unitful: Unitful # extend: unit, dimension
using Interpolations: linear_interpolation
using Plots: @recipe
using Unitful: Dimensions, Frequency, NoDims, Quantity, Hz, dimension, unit, upreferred, 𝐓, 𝐋
using Integrals: IntegralProblem, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
