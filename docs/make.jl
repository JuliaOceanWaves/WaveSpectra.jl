using Documenter, Example

ENV["JULIA_DEBUG"]=Documenter
makedocs(sitename="WaveSpectra")

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
