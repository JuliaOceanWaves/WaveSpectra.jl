using Test, Documenter, WaveSpectra

ENV["UNITFUL_FANCY_EXPONENTS"] = true
# ENV["JULIA_DEBUG"]=Documenter

doctest(WaveSpectra)