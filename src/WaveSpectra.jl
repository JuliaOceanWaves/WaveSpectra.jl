"""
Package for working with wave spectral data.

This is part of the `JuliaOceanWaves` ecosystem.
`WaveSpectra.jl` uses units and is built on top of `DimensionfulAngles.jl`,
treating angles as a dimension.
This avoids common mistakes when working with spectral data.
"""
module WaveSpectra

# using
using Unitful: Dimensions, FreeUnits, NoDims, Quantity, Units, dimension, ustrip, gn as g,
    m, s, 𝐋, 𝐓
using DimensionfulAngles: Dispersion, radᵃ as rad, turnᵃ as τ, θ₀, 𝐀
using UnitfulEquivalences: Equivalence, dimtype, edconvert
using Integrals: AbstractSampledIntegralAlgorithm, SampledIntegralProblem, TrapezoidalRule,
    UniformWeights, solve
using AxisArrays: Axis, AxisArray, ClosedInterval, axisvalues, (..)
using Plots: mm as plots_mm, text, @recipe
using PrettyTables: HtmlTableStyle, pretty_table
using Makie: Axis as MAxis, Colorbar, Figure, PolarAxis, contourf!, lines!, plot!

# imports to overload
import Base  # BroadcastStyle, copy, eltype, getindex, setindex!, show, similar, size
import Unitful: uconvert, unit
import UnitfulEquivalences: edconvert
import AxisArrays # axes # in the future, do `import AxisArrays: axes as AAaxes`
const axes = Base.axes # name conflict will be fixed by AxisArrays in the future
import Integrals: find_weights

# export
export DispersionRelations, ParametricSpectra # modules
export OmnidirectionalSpectrum, RectangularRule, Spectrum # structs
export axesinfo, evenspacing, integrate, isevenlyspaced,  # functions
    plot_spectrum, plot_spectrum!, split_spectrum, spread_function
    # axesnames, axestypes, coordinates, iscartesian, isdirection, ispolar,
    # isspatial, istemporal

# include
include("core.jl")
include("splitspectrum.jl")
include("integration.jl")
include("conversion.jl")
include("plotting.jl")
include("submodules/DispersionRelations.jl")
include("submodules/ParametricSpectra.jl")

end
