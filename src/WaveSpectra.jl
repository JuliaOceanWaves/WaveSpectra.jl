module WaveSpectra

module _NoUnit
    using Unitful: @unit, register
    export ∅
    @unit ∅ "∅" NoUnit 1 false
    function __init__()
        register(_NoUnit)
    end
end

export DiscreteOmnidirectionalSpectrum, OmnidirectionalSpectrum, WaveTimeSeries
export convert_frequency, frequency_dimension, frequency_unit, quantity, spectral_moment
export energy_period, significant_waveheight, ∅
export normalize, scale, slope_spectrum, deepwater, pierson_moskowitz_spectrum

using FFTW
using ._NoUnit: ∅
using Unitful: Unitful # extend: unit, dimension
using Interpolations: linear_interpolation
using Plots: @recipe
using Integrals: AbstractSampledIntegralAlgorithm, IntegralProblem, SampledIntegralProblem, TrapezoidalRule, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode
using Unitful: 𝐓, 𝐋, gn as g, Hz, m
using Unitful: dimension, register, unit, uconvert, upreferred, ustrip
using Unitful: Frequency, Length, Time, Wavenumber
using Unitful: Dimensions, DimensionlessQuantity, Frequency, NoDims, Quantity, Hz, s
using UnitfulEquivalences: UnitfulEquivalences # extend: edconvert
using UnitfulEquivalences: Equivalence, edconvert, dimtype
using DimensionfulAngles: 𝐀, radᵃ as rad, θ₀
using DimensionfulAngles: AngularPeriod, AngularVelocity, AngularWavelength
using DimensionfulAngles: AngularWavenumber, Dispersion

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
