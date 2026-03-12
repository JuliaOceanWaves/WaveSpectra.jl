using Test, SafeTestsets, Documenter

@time @testset verbose=true "WaveSpectra.jl" begin
    @time @safetestset "Test Core Capabilities" begin
        include("test_core.jl")
    end
    @time @safetestset "Test Split Spectrum" begin
        include("test_splitspectrum.jl")
    end
    @time @safetestset "Test Integration" begin
        include("test_integration.jl")
    end
    @time @safetestset "Test Conversion" begin
        include("test_conversion.jl")
    end
    @time @safetestset "Test Dispersion Relations" begin
        include("test_dispersion.jl")
    end
    @time @safetestset "Test Parametric Spectra" begin
        include("test_parametric.jl")
    end
end
