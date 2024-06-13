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
export frequency_dimension, frequency_unit, pm_spectrum, quantity, spectral_moment
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
    (𝐓, 𝐓) => (gf, gb -> (x -> 1)),
    (𝐓^-1, 𝐓^-1) => (gf, gb -> (x -> 1)),
    (𝐀*𝐓^-1, 𝐀*𝐓^-1) => (gf, gb -> (x -> 1)),
    (𝐀^-1*𝐓, 𝐀^-1*𝐓) => (gf, gb -> (x -> 1)),
    (𝐓^-1, 𝐀*𝐓^-1) => (gf, gb -> (x -> 2π * rad)),
    (𝐀*𝐓^-1, 𝐓^-1) => (gf, gb -> (x -> 1 / (2π * rad))),
    (𝐓^-1, 𝐓) => (gf, gb -> (x -> 1 / x^2)),
    (𝐓, 𝐓^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐓, 𝐓*𝐀^-1) => (gf, gb -> (x -> 1 / (2π * rad))),
    (𝐓*𝐀^-1, 𝐓) => (gf, gb -> (x -> 2π * rad)),
    (𝐀*𝐓^-1, 𝐓*𝐀^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐓*𝐀^-1, 𝐀*𝐓^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐓^-1, 𝐓*𝐀^-1) => (gf, gb -> (x -> 1/(2π*rad*x^2))),
    (𝐓*𝐀^-1, 𝐓^-1) => (gf, gb -> (x -> 1/(2π*rad*x^2))),
    (𝐓, 𝐀*𝐓^-1) => (gf, gb -> (x -> 2π*rad/(x^2))),
    (𝐀*𝐓^-1, 𝐓) => (gf, gb -> (x -> 2π*rad/(x^2))),
    (𝐋, 𝐋) => (gf, gb -> (x -> 1)),
    (𝐋^-1, 𝐋^-1) => (gf, gb -> (x -> 1)),
    (𝐀*𝐋^-1, 𝐀*𝐋^-1) => (gf, gb -> (x -> 1)),
    (𝐀^-1*𝐋, 𝐀^-1*𝐋) => (gf, gb -> (x -> 1)),
    (𝐋^-1, 𝐀*𝐋^-1) => (gf, gb -> (x -> 2π * rad)),
    (𝐀*𝐋^-1, 𝐋^-1) => (gf, gb -> (x -> 1 / (2π * rad))),
    (𝐋^-1, 𝐋) => (gf, gb -> (x -> 1 / x^2)),
    (𝐋, 𝐋^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐋, 𝐋*𝐀^-1) => (gf, gb -> (x -> 1 / (2π * rad))),
    (𝐋*𝐀^-1, 𝐋) => (gf, gb -> (x -> 2π * rad)),
    (𝐀*𝐋^-1, 𝐋*𝐀^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐋*𝐀^-1, 𝐀*𝐋^-1) => (gf, gb -> (x -> 1 / x^2)),
    (𝐋^-1, 𝐋*𝐀^-1) => (gf, gb -> (x -> 1 / (2π * rad * x^2))),
    (𝐋*𝐀^-1, 𝐋^-1) => (gf, gb -> (x -> 1 / (2π * rad * x^2))),
    (𝐋, 𝐀 * 𝐋^-1) => (gf, gb -> (x -> 2π * rad / (x^2))),
    (𝐀*𝐋^-1, 𝐋) => (gf, gb -> (x -> 2π * rad / (x^2))),

    (𝐀*𝐋^-1, 𝐀*𝐓^-1) =>(gf, gb -> (x -> gf(x))),
    (𝐀*𝐓^-1, 𝐀*𝐋^-1) =>(gf, gb -> (x -> gb(x))),

)

function _get_grad(dimension_from::Dimensions, dimension_to::Dimensions,
        dispersion_forward_gradient::Union{Function, Nothing},
        dispersion_backward_gradient::Union{Function, Nothing})
    f = _grad[dimension_from, dimension_to]
    return f(dispersion_forward_gradient, dispersion_backward_gradient)
end

include("omnidirectional_spectra.jl")
# include("spreading_functions.jl")
# include("directional_spectra.jl")

end
