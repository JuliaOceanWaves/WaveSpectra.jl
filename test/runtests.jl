using Test, SafeTestsets

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
    @time @safetestset "Test Plotting" begin
        include("test_plotting.jl")
    end

    # submodules
    @time @safetestset "Test Module: DispersionRelations" begin
        include("test_module_dispersion.jl")
    end
    @time @safetestset "Test Module: ParametricSpectra" begin
        include("test_module_parametric.jl")
    end
    @time @safetestset "Test Module: Moments" begin
        include("test_module_moments.jl")
    end
    @time @safetestset "Test Module: Shapes" begin
        include("test_module_shapes.jl")
    end

    # documentation
    @time @safetestset "Doc Tests" begin
        include("test_doctest.jl")
    end
end
