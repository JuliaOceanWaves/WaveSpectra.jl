using Documenter
using DocumenterCitations
using WaveSpectra
using Unitful
using DimensionfulAngles
using AxisArrays

ENV["UNITFUL_FANCY_EXPONENTS"] = true
bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

function regex_match(regex, x)
	return !isnothing(match(regex, string(x)))
end

makedocs(
	sitename = "WaveSpectra.jl",
	format = Documenter.HTML(
		prettyurls = get(ENV, "CI", "false") == "true",
		assets = String["assets/citations.css",
			"assets/images.css"],
	),
	modules = [WaveSpectra],
	pages = [
		"Home" => "index.md",
		"Package Guide" => ["Quickstart" => "guide/quickstart.md",
			"Spectra" => "guide/spectra.md",
			"Working with Spectra" => "guide/spectra_funcs.md",
			"Other Functions" => "guide/other.md"],
		"API" => "api.md",
	],
	plugins = [bib],
)

deploydocs(
	repo = "github.com/JuliaOceanWaves/WaveSpectra.jl.git",
	devbranch = "main",
	push_preview = true,
	versions = ["latest" => "v^", "v#.#"],
)
