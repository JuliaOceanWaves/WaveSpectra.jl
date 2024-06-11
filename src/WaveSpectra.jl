module WaveSpectra

module _NoUnit
    using Unitful: @unit, register
    export ∅
    @unit ∅ "∅" NoUnit 1 false
    function __init__()
        register(_NoUnit)
    end
end

export DiscreteOmnidirectionalSpectrum, OmnidirectionalSpectrum
export integrate, frequency_dimension, frequency_unit, pm_spectrum, quantity, spectral_moment
export energy_period, significant_waveheight

export DiscreteOmnidirectionalSpectrum, OmnidirectionalSpectrum
export convert_frequency, frequency_dimension, frequency_unit, quantity, spectral_moment
export energy_period, significant_waveheight, ∅
export normalize, scale, slope_spectrum, deepwater, pierson_moskowitz_spectrum

using ._NoUnit: ∅
using Unitful: Unitful # extend: unit, dimension
using Interpolations: linear_interpolation
using Plots: @recipe
using Integrals: AbstractSampledIntegralAlgorithm, IntegralProblem, SampledIntegralProblem, TrapezoidalRule, QuadGKJL, solve
using SciMLBase: AbstractIntegralAlgorithm, ReturnCode
using Unitful: 𝐓, 𝐋, gn as g, Hz, m
using Unitful: dimension, register, unit, uconvert, upreferred
using Unitful: Frequency, Length, Time, Wavenumber
using Unitful: Dimensions, DimensionlessQuantity, Frequency, NoDims, Quantity, Hz
using UnitfulEquivalences: Equivalence
using DimensionfulAngles: 𝐀, radᵃ as rad, θ₀
using DimensionfulAngles: AngularPeriod, AngularVelocity, AngularWavelength
using DimensionfulAngles: AngularWavenumber, Dispersion

const _frequency_dims = [𝐓, 𝐓^-1, 𝐓*𝐀^-1, 𝐀*𝐓^-1, 𝐋, 𝐋^-1, 𝐋*𝐀^-1, 𝐀*𝐋^-1, NoDims]
const _grad = Dict( #TODO
    (𝐓, 𝐓) => (x -> 1),
    (𝐓^-1, 𝐓^-1) => (x -> 1),
    (𝐀 * 𝐓^-1, 𝐀 * 𝐓^-1) => (x -> 1),
    (𝐓^-1, 𝐀 * 𝐓^-1) => (x -> 2π * rad),
    (𝐀 * 𝐓^-1, 𝐓^-1) => (x -> 1 / (2π * rad)),
    (𝐓^-1, 𝐓) => (x -> 1 / x^2),
    (𝐓, 𝐓^-1) => (x -> 1 / x^2),
)

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
