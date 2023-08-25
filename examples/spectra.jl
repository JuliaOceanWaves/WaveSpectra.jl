### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 0b630718-c1fe-11ec-0fdc-99e430584e3e
let
	using Unitful: Quantity, Units, dimension, uconvert, unit, 𝐓, 𝐋, Hz, m, s, ustrip
	using DimensionfulAngles: radᵃ as rad, °ᵃ as °, 𝐀
	using UnitfulEquivalences: Equivalence
	using AxisArrays: AxisArray, Axis, axisnames, axisvalues, (..), ClosedInterval
	using NumericalIntegration
end;

# ╔═╡ e2646a89-0080-4f2f-a74c-76019514b6d0
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
		Expr(:toplevel,
			 :(eval(x) = $(Expr(:core, :eval))($name, x)),
			 :(include(x) = $(Expr(:top, :include))($name, x)),
			 :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
			 :(include($path))))
	m
end;

# ╔═╡ 3fe4c14b-7974-4c9c-9fe5-d2e74edcd888
Spectra = ingredients("../src/Spectra.jl").Spectra

# ╔═╡ f6292c5f-80f7-448d-8436-07b646d78b78
md"""
## 1. Spectrum
"""

# ╔═╡ 0cfbad4f-fb82-470d-8b0f-d9dbf362b16f
begin
	Nf = 20
	Δf = 0.1Hz
	f = (0:Nf-1)*Δf
end

# ╔═╡ a63202c5-ff50-4a7f-bc60-228bdc48455f
begin
	Nθ = 36
	Δθ = 360°/Nθ
	θ = (0:Nθ-1)*Δθ
end

# ╔═╡ f7889a0f-f846-402d-95b9-1b9abe05f3a5
S = Spectra.Spectrum(randn(Nf, Nθ)*m^2/Hz/°, f, θ)

# ╔═╡ 97522a2c-e43a-415f-8de3-7a13ad5572a0
md"### indexing"

# ╔═╡ 1833e4b8-ea80-4e04-a6f3-b58e614166d0
axesnames(S)

# ╔═╡ f3dfafdb-56b6-445d-8653-d1692d7a234b
S[direction=40° .. 80°]

# ╔═╡ dccbbfdb-aece-43fb-bf95-c35c8fdf721a
S[frequency=0.2Hz .. 0.5Hz]

# ╔═╡ 06750eb9-794e-42fa-ae52-d7121b9b9b2f
S[direction=40° .. 80°, frequency=0.2Hz .. 0.5Hz]

# ╔═╡ 194b0c5f-aa6b-4dfb-a26c-49bdcdcdbd2d
S[direction=10°]

# ╔═╡ 8340e2b2-caa6-4806-a404-01bc81e74335
S[direction=10°, frequency=0.2Hz]

# ╔═╡ 6bd5fe2f-12db-4851-8e7f-f38b995f5eb0
md"### assignment"

# ╔═╡ cfe1e2e3-23d0-45b0-9525-f061bba17c3f
S

# ╔═╡ 44a3569c-21a1-4fa6-962e-d7d8da88bf36
begin
	S2 = copy(S)
	S2[frequency=0.2Hz .. 0.3Hz, direction=10.0°] *= 10
	S2[frequency=0.2Hz, direction=10.0°] *= 100
	S2[frequency=0.2Hz] *= 0
	S2[frequency=0.2Hz, direction=10.0°] = 100m^2/Hz/°
end

# ╔═╡ c8c2dfa9-c771-4c7b-ae52-92c6a63f0ddc
S2

# ╔═╡ bb90dbfd-a5ca-42e3-9d60-018960034584
md"""
#### AxisArray bug?
Can't do similar assignment with `Spectrum`.
"""

# ╔═╡ 43e98356-5ec0-4d96-a614-934942f5aa88
B = AxisArray(randn(3,2))

# ╔═╡ 1f6590b0-7ae3-488b-9b97-e7732690ec07
B[row=1..2, col=2] = [0.2, 0.3]

# ╔═╡ ee09bc43-5014-4977-bdfb-05102d036e00
md"### units"

# ╔═╡ 666dc7e3-8e36-472d-aa08-bba54087b518
unit(S), unit(S, :frequency), unit(S, :direction), unit(S, :integral)

# ╔═╡ bad687c8-e0c1-4a18-87d1-4fb2ccfb82e2
md"### arrays of Spectrum"

# ╔═╡ 254082a1-4446-47ad-92db-adaf5fec84b2
[S,S]

# ╔═╡ 7c33dd32-a7f7-436c-99c8-a4b35e18c71b
[[S] [S]]

# ╔═╡ 6ac0f90b-eb42-4f9f-a7c9-fdec45641075
[[S] [S];;;[S] [S]]

# ╔═╡ ec4ffbcc-8408-4358-82cd-1aad39a3249d
AxisArray([[S] [S]; [S] [S]], location=[:Hawaii, :PuertoRico], season=[:winter, :summer])

# ╔═╡ f631f0c0-0f1c-477e-b36d-fca554cd93a0
AxisArray(S)

# ╔═╡ 136f4823-f6eb-4f4d-a4a3-cf891a8b5506
md"""
## 2. Omni-directional spectrum
"""

# ╔═╡ a0b35a84-a8e2-4ca7-9790-d908f2cae323
begin
	Nω = 10
	Δω = 0.1rad/s
	ω = (0:Nω-1)*Δω
end

# ╔═╡ 03126da0-7416-457a-b081-e676e0808042
Sω = Spectra.OmniSpectrum(randn(Nω) * m^2/(rad/s), ω)

# ╔═╡ 7c0b52ec-db20-4e15-939a-686967020046
md"### indexing"

# ╔═╡ 132a12f9-32c9-4458-a078-94bce25b5983
Spectra.axesname(Sω)

# ╔═╡ 81508fa5-c8ee-4a7b-aa9d-5807385a864d
Sω[angular_frequency=0.1rad/s]

# ╔═╡ 4bf04ae1-8f9c-4a32-99d1-c40356a649cd
Sω[angular_frequency=0.1rad/s .. 0.5rad/s]

# ╔═╡ c7b33e7e-f8ce-4029-b42d-54619dec6099
Sω[0.2rad/s]

# ╔═╡ 79fff198-88a4-4dbb-a367-5a6027436302
Sω[0.2rad/s .. 0.4rad/s]

# ╔═╡ ef01bfc1-f134-43b0-a3ff-b74007f55701
md"### assignment"

# ╔═╡ 8017e36c-e305-4693-b895-1c4e4333e6ae
Sω

# ╔═╡ abc73aeb-5316-4798-950a-c25e8ad6898b
begin
	Sω2 = copy(Sω)
	Sω2[angular_frequency=0.2rad/s] = 100m^2/(rad/s)
	Sω2[0.3rad/s] = 200m^2/(rad/s)
end;

# ╔═╡ 938a93f6-7272-411a-b1a4-d20943720292
Sω2

# ╔═╡ 8dce359a-5111-4eb3-ab8a-57b8a548d2ca
md"#### AxisArray bug?"

# ╔═╡ df655a49-a852-4fed-805f-abb7121de775
Sω2[0.8rad/s .. 0.9rad/s] = [0m^2/(rad/s), 0m^2/(rad/s)]

# ╔═╡ 70c80fab-62a6-4eb0-93a6-ef6692cbd1b9
md"### units"

# ╔═╡ 443ecb7d-c04d-4999-9414-73f4d5e9601f
unit(Sω), unit(Sω, :angular_frequency), unit(Sω, :integral)

# ╔═╡ 63b90a22-4d65-4952-a58a-fe8151cba3b0
md"### arrays of OmniSpectrum"

# ╔═╡ 087c300e-73eb-49d0-bc7c-c19882ef8192
[Sω, Sω]

# ╔═╡ 62cb67ce-5a29-4976-ba55-cf3084de3323
[[Sω] [Sω]; [Sω] [Sω]]

# ╔═╡ c02e45da-4360-4732-833b-405b4cd151b0
[[Sω] [Sω];;; [Sω] [Sω]]

# ╔═╡ 3635ca3f-ef64-41e8-a179-a98b473cef01
AxisArray([[Sω] [Sω]; [Sω] [Sω]], location=[:Hawaii, :PuertoRico], season=[:winter, :summer])

# ╔═╡ cb4b004f-f606-4bae-bb21-7a3a01524300
AxisArray(Sω)

# ╔═╡ 74f23746-8a72-4133-8162-e903df905b27
[S, Sω]

# ╔═╡ 57b415de-a173-41f7-95ce-e3250b272a96
(S, Sω)

# ╔═╡ 688fbba9-4950-4158-bcc9-6c81c9a530fa
md"## 3. Integrals"

# ╔═╡ 04865b34-240b-4891-93b2-1c7abe52864a
integrate(S), integrate(Sω)

# ╔═╡ f102117e-d4ee-44b5-b6ad-b6e98724762e
integrate(Sω), integrate(Sω, TrapezoidalEvenFast()), integrate(Sω, SimpsonEven())

# ╔═╡ 1bd83ee2-a95f-443f-a138-6ff903b745a7
integrate(S, :frequency)

# ╔═╡ dae4dbfa-cc82-4c66-b8b3-f4affec8cb67
integrate(S, :direction, TrapezoidalEvenFast())

# ╔═╡ 010a824e-0671-4748-ab92-5b416e8bd0ad
md"## 4. Change of variables"

# ╔═╡ 3629ad15-b992-475a-a94a-ea3270311a70


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AxisArrays = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
DimensionfulAngles = "2d4b8d7a-02d9-40f9-9abe-9c695b77de0d"
NumericalIntegration = "e7bfaba1-d571-5449-8927-abc22e82249b"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulEquivalences = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"

[compat]
AxisArrays = "~0.4.6"
DimensionfulAngles = "~0.2.0"
NumericalIntegration = "~0.3.3"
Unitful = "~1.15.0"
UnitfulEquivalences = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "6a726d212360dfc94757feecca73d17368890b53"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "4e88377ae7ebeaf29a047aa1ee40826e0b708a5d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.7.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DimensionfulAngles]]
deps = ["Unitful", "UnitfulEquivalences"]
git-tree-sha1 = "2a485014726b128d8faa59064ea2aa68c8094a5d"
uuid = "2d4b8d7a-02d9-40f9-9abe-9c695b77de0d"
version = "0.2.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NumericalIntegration]]
deps = ["Interpolations", "LinearAlgebra", "Logging"]
git-tree-sha1 = "2a4ef5fc235053f9747d59cfdee19bcb8ba1e833"
uuid = "e7bfaba1-d571-5449-8927-abc22e82249b"
version = "0.3.3"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

    [deps.Ratios.weakdeps]
    FixedPointNumbers = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore"]
git-tree-sha1 = "fffc14c695c17bfdbfa92a2a01836cdc542a1e46"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.6.1"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "1d5708d926c76a505052d0d24a846d5da08bc3a4"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c4d2a349259c8eba66a00a540d550f122a3ab228"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.15.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulEquivalences]]
deps = ["Unitful"]
git-tree-sha1 = "76fc2f7fdc87531a1018eb7d647df7c29daf36b7"
uuid = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"
version = "0.2.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═0b630718-c1fe-11ec-0fdc-99e430584e3e
# ╠═e2646a89-0080-4f2f-a74c-76019514b6d0
# ╠═3fe4c14b-7974-4c9c-9fe5-d2e74edcd888
# ╟─f6292c5f-80f7-448d-8436-07b646d78b78
# ╠═0cfbad4f-fb82-470d-8b0f-d9dbf362b16f
# ╠═a63202c5-ff50-4a7f-bc60-228bdc48455f
# ╠═f7889a0f-f846-402d-95b9-1b9abe05f3a5
# ╟─97522a2c-e43a-415f-8de3-7a13ad5572a0
# ╠═1833e4b8-ea80-4e04-a6f3-b58e614166d0
# ╠═f3dfafdb-56b6-445d-8653-d1692d7a234b
# ╠═dccbbfdb-aece-43fb-bf95-c35c8fdf721a
# ╠═06750eb9-794e-42fa-ae52-d7121b9b9b2f
# ╠═194b0c5f-aa6b-4dfb-a26c-49bdcdcdbd2d
# ╠═8340e2b2-caa6-4806-a404-01bc81e74335
# ╟─6bd5fe2f-12db-4851-8e7f-f38b995f5eb0
# ╠═cfe1e2e3-23d0-45b0-9525-f061bba17c3f
# ╠═44a3569c-21a1-4fa6-962e-d7d8da88bf36
# ╠═c8c2dfa9-c771-4c7b-ae52-92c6a63f0ddc
# ╟─bb90dbfd-a5ca-42e3-9d60-018960034584
# ╠═43e98356-5ec0-4d96-a614-934942f5aa88
# ╠═1f6590b0-7ae3-488b-9b97-e7732690ec07
# ╟─ee09bc43-5014-4977-bdfb-05102d036e00
# ╠═666dc7e3-8e36-472d-aa08-bba54087b518
# ╟─bad687c8-e0c1-4a18-87d1-4fb2ccfb82e2
# ╠═254082a1-4446-47ad-92db-adaf5fec84b2
# ╠═7c33dd32-a7f7-436c-99c8-a4b35e18c71b
# ╠═6ac0f90b-eb42-4f9f-a7c9-fdec45641075
# ╠═ec4ffbcc-8408-4358-82cd-1aad39a3249d
# ╠═f631f0c0-0f1c-477e-b36d-fca554cd93a0
# ╟─136f4823-f6eb-4f4d-a4a3-cf891a8b5506
# ╠═a0b35a84-a8e2-4ca7-9790-d908f2cae323
# ╠═03126da0-7416-457a-b081-e676e0808042
# ╟─7c0b52ec-db20-4e15-939a-686967020046
# ╠═132a12f9-32c9-4458-a078-94bce25b5983
# ╠═81508fa5-c8ee-4a7b-aa9d-5807385a864d
# ╠═4bf04ae1-8f9c-4a32-99d1-c40356a649cd
# ╠═c7b33e7e-f8ce-4029-b42d-54619dec6099
# ╠═79fff198-88a4-4dbb-a367-5a6027436302
# ╟─ef01bfc1-f134-43b0-a3ff-b74007f55701
# ╠═8017e36c-e305-4693-b895-1c4e4333e6ae
# ╠═abc73aeb-5316-4798-950a-c25e8ad6898b
# ╠═938a93f6-7272-411a-b1a4-d20943720292
# ╟─8dce359a-5111-4eb3-ab8a-57b8a548d2ca
# ╠═df655a49-a852-4fed-805f-abb7121de775
# ╟─70c80fab-62a6-4eb0-93a6-ef6692cbd1b9
# ╠═443ecb7d-c04d-4999-9414-73f4d5e9601f
# ╟─63b90a22-4d65-4952-a58a-fe8151cba3b0
# ╠═087c300e-73eb-49d0-bc7c-c19882ef8192
# ╠═62cb67ce-5a29-4976-ba55-cf3084de3323
# ╠═c02e45da-4360-4732-833b-405b4cd151b0
# ╠═3635ca3f-ef64-41e8-a179-a98b473cef01
# ╠═cb4b004f-f606-4bae-bb21-7a3a01524300
# ╠═74f23746-8a72-4133-8162-e903df905b27
# ╠═57b415de-a173-41f7-95ce-e3250b272a96
# ╟─688fbba9-4950-4158-bcc9-6c81c9a530fa
# ╠═04865b34-240b-4891-93b2-1c7abe52864a
# ╠═f102117e-d4ee-44b5-b6ad-b6e98724762e
# ╠═1bd83ee2-a95f-443f-a138-6ff903b745a7
# ╠═dae4dbfa-cc82-4c66-b8b3-f4affec8cb67
# ╟─010a824e-0671-4748-ab92-5b416e8bd0ad
# ╠═3629ad15-b992-475a-a94a-ea3270311a70
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
