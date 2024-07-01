using Documenter, DocumenterInterLinks, Example
using WaveSpectra, Unitful, DimensionfulAngles

links = InterLinks(
    "DimensionfulAngles" => "https://cmichelenstrofer.github.io/DimensionfulAngles.jl/stable/",
    "Unitful" => "https://painterqubits.github.io/Unitful.jl/stable/",
)
ENV["UNITFUL_FANCY_EXPONENTS"] = true
#ENV["JULIA_DEBUG"]=Documenter
DocMeta.setdocmeta!(WaveSpectra, :DocTestSetup, :(using WaveSpectra); recursive=true)
makedocs(sitename="WaveSpectra",
         modules=[WaveSpectra],
         plugins=[links])

deploydocs(repo="github.com/JuliaOceanWaves/WaveSpectra.jl.git",)
