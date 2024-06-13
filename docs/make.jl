using Documenter, DocumenterInterLinks, Example
using WaveSpectra, Unitful, DimensionfulAngles

links = InterLinks(
    "DimensionfulAngles" => "https://cmichelenstrofer.github.io/DimensionfulAngles.jl/stable/",
    "Unitful" => "https://painterqubits.github.io/Unitful.jl/stable/",
)

#ENV["JULIA_DEBUG"]=Documenter
makedocs(sitename="WaveSpectra", plugins=[links])

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
