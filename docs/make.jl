using Documenter, Example

JULIA_DEBUG=Documenter
makedocs(sitename="WaveSpectra")

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
