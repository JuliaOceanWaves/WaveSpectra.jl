module Spectra

using Unitful: Quantity, Dimensions, Units, Quantity, dimension, 𝐓, 𝐋, m, s, gn as g
using DimensionfulAngles: 𝐀, Dispersion, radᵃ as rad, θ₀
using UnitfulEquivalences: Equivalence
using Integrals: SampledIntegralProblem, AbstractSampledIntegralAlgorithm, TrapezoidalRule, solve
using AxisArrays: AxisArray, Axis, ClosedInterval, axisvalues, (..)

import Base  # size, getindex, setindex!, copy, similar, eltype, BroadcastStyle
import Unitful: unit, uconvert
import AxisArrays # axes
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future

export Spectrum, OmnidirectionalSpectrum
export axesinfo, integrate, spread_function, omnidirectional_spectrum, split_spectrum
export deepwater
# export ispolar, iscartesian, istemporal, isspatial, isdirection, axestypes, coordinates, axesnames

include("core.jl")
include("functions.jl")
include("integration.jl")
include("convert.jl")
include("dispersion.jl")

end
