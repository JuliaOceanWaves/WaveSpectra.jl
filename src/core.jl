# core functionality to represent directional and omnidirectional spectra

abstract type AbstractSpectrum{TDAT} <: AbstractSuperposition{TDAT} end
abstract type AbstractOmnidirectionalSpectrum{TDAT} <: AbstractSuperposition1D{TDAT} end

include("core/spectrum.jl")
include("core/omnidirectional.jl")
include("core/utilities.jl")
