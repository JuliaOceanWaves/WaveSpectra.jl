module WaveSpectra

export DiscreteOmnidirectionalSpectrum, OmnidirectionalSpectrum
export integrate, quantity, spectral_moment
# export energy_period, significant_waveheight

# using Base: Base # extend: showerror
using Unitful: Unitful # extend: unit, dimension
using Interpolations: linear_interpolation
using Plots: @recipe
using Integrals: AbstractSampledIntegralAlgorithm, IntegralProblem, SampledIntegralProblem, TrapezoidalRule, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode


using Unitful: 𝐓, 𝐋, Hz
using Unitful: dimension, unit, upreferred
using Unitful: Frequency, Length, Time, Wavenumber
using Unitful: Dimensions, Frequency, NoDims, Quantity, Hz
using DimensionfulAngles: 𝐀
using DimensionfulAngles: AngularPeriod, AngularVelocity, AngularWavelength
using DimensionfulAngles: AngularWavenumber, Periodic


# const _temporal_dims = [Time, Frequency, AngularPeriod, AngularVelocity]
# const _spatial_dims = [Length, Wavenumber, AngularWavelength, AngularWavenumber]
# const _frequency_dims = vcat(_temporal_dims, _spatial_dims)

const _frequency_quantity = Dict(
    𝐓 => Time, 𝐓^-1 => Frequency, 𝐓*𝐀^-1 => AngularPeriod, 𝐀*𝐓^-1 => AngularVelocity,
    𝐋 => Length, 𝐋^-1 => Wavenumber, 𝐋*𝐀^-1 => AngularWavelength,
    𝐀*𝐋^-1 => AngularWavenumber
)

const _frequency_dims = [𝐓, 𝐓^-1, 𝐓*𝐀^-1, 𝐀*𝐓^-1, 𝐋, 𝐋^-1, 𝐋*𝐀^-1, 𝐀*𝐋^-1]


# const _TemporalDims = Union{Time, Frequency, AngularPeriod, AngularVelocity}
# const _SpatialDims = Union{Length, Wavenumber, AngularWavelength, AngularWavenumber}
# const _FrequencyDims = Union{_TemporalDims, _SpatialDims}

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
