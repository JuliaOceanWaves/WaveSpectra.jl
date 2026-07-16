"""
Package for working with wave spectral data.

This is part of the `JuliaOceanWaves` ecosystem.
`WaveSpectra.jl` uses units and is built on top of `DimensionfulAngles.jl`,
treating angles as a dimension.
This avoids common mistakes when working with spectral data.
"""
module WaveSpectra

# using
using Unitful: Dimensions, Quantity, Units, dimension, gn as g, Hz, kg, m, s, 𝐋, 𝐓
using DimensionfulAngles: Dispersion, θ₀, 𝐀
using DimensionfulAngles.DefaultSymbols: rad, °
using Integrals: AbstractSampledIntegralAlgorithm, SampledIntegralProblem, TrapezoidalRule,
                 UniformWeights, solve
using SpectralSuperpositions: AbstractSuperposition, AbstractSuperposition1D, axesinfo,
                              axesnames, axestypes, coordinates, isangular, iscartesian,
                              isdirection, isfrequency, islinear, isperiod, ispolar,
                              isspatial, isspectralvariable, istemporal, isevenlyspaced,
                              evenspacing, validate_superposition, validate_superposition1d,
                              rebuild_superposition
using AxisArrays: Axis, ClosedInterval, axisvalues, (..)

# imports to overload
import Base  # BroadcastStyle, copy, eltype, getindex, isapprox, setindex!, show, similar,
# size, (==), (!=), (<), (<=), (>), (>=)
import Unitful: uconvert, unit
import AxisArrays: AxisArrays, AxisArray # axes # in the future, do `import AxisArrays: axes as AAaxes`
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future
import Integrals: find_weights
import SpectralSuperpositions: superposition_unit_aliases

function plot_spectrum end
function plot_spectrum! end

# export
export DispersionRelations, Moments, ParametricSpectra, Shapes # modules
export OmnidirectionalSpectrum, RectangularRule, Spectrum # structs
export integrate, isspread, plot_spectrum, plot_spectrum!, polar_to_cartesian, # functions
       cartesian_to_polar, split_spectrum, spread_function
export g, Hz, kg, m, periodic, rad, s, uconvert, unit, # reexport from other packages
       °, (..)

# include
include("core.jl")
include("splitspectrum.jl")
include("integration.jl")
include("conversion.jl")
include("submodules/DispersionRelations.jl")
include("submodules/ParametricSpectra.jl")
include("submodules/Moments.jl")
include("submodules/Shapes.jl")

end
