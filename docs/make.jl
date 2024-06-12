using Documenter, Example
using WaveSpectra, Unitful

#ENV["JULIA_DEBUG"]=Documenter
makedocs(sitename="WaveSpectra")

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
