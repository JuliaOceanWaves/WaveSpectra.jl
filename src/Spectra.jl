module Spectra

using Unitful: Dimensions, Frequency, Length, Quantity, Time, Units, Wavenumber
using Unitful: dimension, gn as g, m, s, 𝐋, 𝐓
using DimensionfulAngles: AngularPeriod, AngularVelocity, AngularWavelength
using DimensionfulAngles: AngularWavenumber, Dispersion, radᵃ as rad, θ₀, 𝐀
using UnitfulEquivalences: Equivalence, dimtype, edconvert
using Integrals: AbstractSampledIntegralAlgorithm, SampledIntegralProblem, TrapezoidalRule
using Integrals: UniformWeights, solve
using AxisArrays: Axis, AxisArray, ClosedInterval, axisvalues, (..)

import Base  # BroadcastStyle, copy, eltype, getindex, setindex!, show, similar, size
import Unitful: uconvert, unit
import UnitfulEquivalences: edconvert
import AxisArrays # axes
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future
import Integrals: find_weights

export OmnidirectionalSpectrum, RectangularRule, Spectrum
export axesinfo, deepwater, integrate, isevenlyspaced, omnidirectional_spectrum
export split_spectrum, spread_function
# export DispersionGradient, axesnames, axestypes, coordinates, iscartesian, isdirection
# export ispolar, isspatial, istemporal

include("core.jl")
include("utilities.jl")
include("integration.jl")
include("conversion.jl")
include("functions.jl")

end
