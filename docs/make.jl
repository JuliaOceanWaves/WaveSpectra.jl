using Documenter, Example
using WaveSpectra

#ENV["JULIA_DEBUG"]=Documenter
makedocs(sitename="WaveSpectra")

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
