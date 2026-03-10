# core functionality to represent directional and omnidirectional spectra

abstract type AbstractSpectrum{TDAT} <: AbstractMatrix{TDAT} end
abstract type AbstractOmnidirectionalSpectrum{TDAT} <: AbstractVector{TDAT} end

include("core/spectrum.jl")
include("core/omnidirectional.jl")
include("core/utilities.jl")
include("core/axes.jl")
include("core/show.jl")
