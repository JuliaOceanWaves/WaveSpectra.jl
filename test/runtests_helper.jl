using Test, Documenter


@time @testset verbose=true "WaveSpectra.jl" begin
    # Write your tests here.
    @time @testset "Doc Tests" begin include("test_doctest.jl") end
    @time @testset "Function Tests" begin include("func_test.jl") end
end
