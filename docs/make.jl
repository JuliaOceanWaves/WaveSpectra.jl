using Documenter
using DocumenterCitations
using WaveSpectra
using Unitful
using DimensionfulAngles

ENV["UNITFUL_FANCY_EXPONENTS"] = true
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
    sitename = "WaveSpectra.jl",
    format = Documenter.HTML(assets = String["assets/citations.css"]),
    modules = [WaveSpectra],
    pages = [
        "Home" => "index.md",
        "Quick Start" => "quickstart.md",
        "Theory" => "theory.md",
        "API" => "api.md",
    ],
    plugins = [bib]
)

deploydocs(; repo = "github.com/JuliaOceanWaves/WaveSpectra.jl.git")
