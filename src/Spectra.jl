module Spectra

using Unitful: Quantity, Units, Dimensions, Frequency, dimension, uconvert, ustrip, unit, 𝐓, 𝐋, Hz, m, s
using DimensionfulAngles: radᵃ as rad, °ᵃ as °, 𝐀, AngularWavenumber, AngularVelocity
using UnitfulEquivalences: Equivalence, dimtype
using AxisArrays: AxisArray, Axis, AbstractInterval, axisnames, axisvalues, (..)
using NumericalIntegration: Trapezoidal, IntegrationMethod
using PhysicalConstants.CODATA2018: StandardAccelerationOfGravitation as g

import Base  # size, getindex, setindex!, copy, show, convert
import Unitful: unit
import AxisArrays: axisnames
import NumericalIntegration: integrate
# import ForwardDiff: derivative

export Spectrum, OmniSpectrum

include("core.jl")
include("show.jl")
include("functions.jl")

end
