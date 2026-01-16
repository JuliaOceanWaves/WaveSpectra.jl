using Documenter
using WaveSpectra

DocMeta.setdocmeta!(WaveSpectra, :DocTestSetup, :(using WaveSpectra); recursive = true)

makedocs(
    modules = [WaveSpectra],
    sitename = "WaveSpectra.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaOceanWaves/WaveSpectra.jl",
    devbranch = "main"
)